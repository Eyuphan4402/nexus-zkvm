#![deny(unsafe_code)]

use ark_crypto_primitives::sponge::constraints::{CryptographicSpongeVar, SpongeWithGadget};
use ark_ec::{
    short_weierstrass::{Projective, SWCurveConfig},
    AdditiveGroup,
};
use ark_ff::{Field, PrimeField};
use ark_r1cs_std::{
    boolean::Boolean,
    eq::EqGadget,
    fields::{fp::FpVar, FieldVar},
    R1CSVar,
};
use ark_relations::r1cs::SynthesisError;
use ark_spartan::polycommitments::PolyCommitmentScheme;
use ark_std::ops::Neg;

pub(crate) mod primary;

use crate::{
    commitment::CommitmentScheme,
    folding::hypernova::{
        cyclefold::nimfs::SQUEEZE_ELEMENTS_BIT_SIZE,
        ml_sumcheck::protocol::verifier::SQUEEZE_NATIVE_ELEMENTS_NUM,
    },
    gadgets::{
        cyclefold::secondary,
        emulated::{cast_field_element_unique, short_weierstrass::EmulatedFpAffineVar},
    },
};

pub fn multifold<G1, G2, C1, C2, RO>(
    config: &<RO::Var as CryptographicSpongeVar<G1::ScalarField, RO>>::Parameters,
    vk: &FpVar<G1::ScalarField>,
    sumcheck_rounds: usize,
    U: &primary::LCCSInstanceFromR1CSVar<G1, C1>,
    U_secondary: &secondary::RelaxedR1CSInstanceVar<G2, C2>,
    u: &primary::CCSInstanceFromR1CSVar<G1, C1>,
    commitment_W_proof: &secondary::ProofVar<G2, C2>,
    hypernova_proof: &primary::ProofFromR1CSVar<G1, RO>,
    should_enforce: &Boolean<G1::ScalarField>,
) -> Result<
    (
        primary::LCCSInstanceFromR1CSVar<G1, C1>,
        secondary::RelaxedR1CSInstanceVar<G2, C2>,
    ),
    SynthesisError,
>
where
    G1: SWCurveConfig<BaseField = G2::ScalarField, ScalarField = G2::BaseField>,
    G2: SWCurveConfig,
    C1: PolyCommitmentScheme<Projective<G1>>,
    C2: CommitmentScheme<Projective<G2>>,
    G1::BaseField: PrimeField,
    G2::BaseField: PrimeField,
    RO: SpongeWithGadget<G1::ScalarField>,
{
    let cs = U.cs();
    let mut random_oracle = RO::Var::new(cs.clone(), config);

    // Absorbing the necessary values into the random oracle
    random_oracle.absorb(&U_secondary)?;
    random_oracle.absorb(&vk)?;
    random_oracle.absorb(&U.var())?;
    random_oracle.absorb(&u.var())?;

    let (rho, rho_bits) = random_oracle
        .squeeze_emulated_field_elements_with_sizes::<G1::BaseField>(&[SQUEEZE_ELEMENTS_BIT_SIZE])?;
    let rho = &rho[0];
    let rho_bits = &rho_bits[0];
    let rho_scalar = Boolean::le_bits_to_fp(rho_bits)?;

    // HyperNova Verification Circuit - implementation specific to R1CS origin for constraints
    const NUM_MATRICES: usize = 3;
    const NUM_MULTISETS: usize = 2;
    const MAX_CARDINALITY: usize = 2;

    // Generating gamma and beta values from the oracle
    let gamma: FpVar<G1::ScalarField> = random_oracle.squeeze_field_elements(1)?[0].clone();
    let beta: Vec<FpVar<G1::ScalarField>> =
        random_oracle.squeeze_field_elements(sumcheck_rounds)?;

    // Generate gamma powers
    let gamma_powers: Vec<FpVar<G1::ScalarField>> = (1..=NUM_MATRICES)
        .map(|j| gamma.pow_le(&Boolean::constant_vec_from_bytes(&j.to_le_bytes())))
        .collect::<Result<Vec<FpVar<G1::ScalarField>>, SynthesisError>>()?;

    // Initial expected value based on gamma powers and instance
    let mut expected: FpVar<G1::ScalarField> = gamma_powers.iter().zip(U.var().vs.iter()).fold(
        FpVar::<G1::ScalarField>::Constant(G1::ScalarField::ZERO),
        |acc, (a, b)| acc + (a * b),
    );

    // Interpolation constants for Lagrange interpolation
    let interpolation_constants = [
        (G1::ScalarField::from(0), G1::ScalarField::from(-6)),
        (G1::ScalarField::from(1), G1::ScalarField::from(2)),
        (G1::ScalarField::from(2), G1::ScalarField::from(-2)),
        (G1::ScalarField::from(3), G1::ScalarField::from(6)),
    ];

    // Absorb the polynomial information
    random_oracle.absorb(&hypernova_proof.var().poly_info.var())?;

    let mut rs_p: Vec<FpVar<G1::ScalarField>> = vec![];
    for round in 0..sumcheck_rounds {
        // Absorb the sumcheck proof for the current round
        random_oracle.absorb(&hypernova_proof.var().sumcheck_proof[round])?;
        let r = random_oracle.squeeze_field_elements(SQUEEZE_NATIVE_ELEMENTS_NUM)?[0].clone();
        random_oracle.absorb(&r)?;

        // Evaluate and enforce expected value consistency
        let evals = &hypernova_proof.var().sumcheck_proof[round];
        expected.conditional_enforce_equal(&(&evals[0] + &evals[1]), should_enforce)?;

        // Lagrange interpolation to evaluate the polynomial
        let prod: FpVar<G1::ScalarField> = (0..(MAX_CARDINALITY + 2)).fold(
            FpVar::<G1::ScalarField>::Constant(G1::ScalarField::ONE),
            |acc, i| acc * (&r - interpolation_constants[i].0),
        );

        expected = (0..(MAX_CARDINALITY + 2))
            .map(|i| {
                let num = &prod * &evals[i];
                let denom = (&r - interpolation_constants[i].0) * interpolation_constants[i].1;
                num.mul_by_inverse(&denom)
            })
            .collect::<Result<Vec<FpVar<G1::ScalarField>>, SynthesisError>>()?
            .iter()
            .sum();

        rs_p.push(r);
    }

    // Calculate the final values for e1 and e2
    let e1 = (0..U.var().rs.len())
        .map(|i| {
            &U.var().rs[i] * &rs_p[i]
                + (FpVar::<G1::ScalarField>::Constant(G1::ScalarField::ONE) - &U.var().rs[i])
                    * (FpVar::<G1::ScalarField>::Constant(G1::ScalarField::ONE) - &rs_p[i])
        })
        .fold(
            FpVar::<G1::ScalarField>::Constant(G1::ScalarField::ONE),
            |acc, x| acc * x,
        );

    let e2 = (0..beta.len())
        .map(|i| {
            &beta[i] * &rs_p[i]
                + (FpVar::<G1::ScalarField>::Constant(G1::ScalarField::ONE) - &beta[i])
                    * (FpVar::<G1::ScalarField>::Constant(G1::ScalarField::ONE) - &rs_p[i])
        })
        .fold(
            FpVar::<G1::ScalarField>::Constant(G1::ScalarField::ONE),
            |acc, x| acc * x,
        );

    // Combine results to compute cl and final constraints
    let cl = gamma_powers
        .iter()
        .zip(hypernova_proof.var().sigmas.iter())
        .fold(
            FpVar::<G1::ScalarField>::Constant(G1::ScalarField::ZERO),
            |acc, (a, b)| acc + (a * b),
        )
        * e1;

    let cSs = [
        (G1::ScalarField::ONE, vec![0, 1]),
        (G1::ScalarField::from(2), vec![1, 3]),
    ];

    // Return the result with necessary constraints
    
