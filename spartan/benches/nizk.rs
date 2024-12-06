#![allow(clippy::assertions_on_result_states)]

extern crate core;
extern crate criterion;
extern crate digest;
extern crate libspartan;
extern crate merlin;

use ark_bls12_381::G1Projective;
use ark_ec::CurveGroup;
use libspartan::{Instance, NIZKGens, NIZK};
use merlin::Transcript;
use criterion::*;

// Helper function to run the NIZK prove benchmark
fn run_nizk_prove_benchmark<G: CurveGroup>(c: &mut Criterion, s: usize) {
    let num_vars = (2_usize).pow(s as u32);
    let num_cons = num_vars;
    let num_inputs = 10;

    let (inst, vars, inputs) =
        Instance::<G::ScalarField>::produce_synthetic_r1cs(num_cons, num_vars, num_inputs);

    let gens = NIZKGens::<G>::new(num_cons, num_vars, num_inputs);

    let mut prover_transcript = Transcript::new(b"example");
    c.bench_function(&format!("NIZK_prove_{}", num_vars), move |b| {
        b.iter(|| {
            NIZK::prove(
                black_box(&inst),
                black_box(vars.clone()),
                black_box(&inputs),
                black_box(&gens),
                black_box(&mut prover_transcript),
            );
        });
    });
}

// Helper function to run the NIZK verify benchmark
fn run_nizk_verify_benchmark<G: CurveGroup>(c: &mut Criterion, s: usize) {
    let num_vars = (2_usize).pow(s as u32);
    let num_cons = num_vars;
    let num_inputs = 10;

    let (inst, vars, inputs) =
        Instance::<G::ScalarField>::produce_synthetic_r1cs(num_cons, num_vars, num_inputs);

    let gens = NIZKGens::<G>::new(num_cons, num_vars, num_inputs);

    // Produce a proof of satisfiability
    let mut prover_transcript = Transcript::new(b"example");
    let proof = NIZK::prove(&inst, vars.clone(), &inputs, &gens, &mut prover_transcript);

    let mut verifier_transcript = Transcript::new(b"example");
    c.bench_function(&format!("NIZK_verify_{}", num_cons), move |b| {
        b.iter(|| {
            assert!(proof
                .verify(
                    black_box(&inst),
                    black_box(&inputs),
                    black_box(&mut verifier_transcript),
                    black_box(&gens)
                )
                .is_ok());
        });
    });
}

// Benchmark for NIZK proving
fn nizk_prove_benchmark<G: CurveGroup>(c: &mut Criterion) {
    for &s in [10, 12, 16].iter() {
        let plot_config = PlotConfiguration::default().summary_scale(AxisScale::Logarithmic);
        let mut group = c.benchmark_group("NIZK_prove_benchmark");
        group.plot_config(plot_config);

        run_nizk_prove_benchmark::<G>(c, s);

        group.finish();
    }
}

// Benchmark for NIZK verification
fn nizk_verify_benchmark<G: CurveGroup>(c: &mut Criterion) {
    for &s in [10, 12, 16].iter() {
        let plot_config = PlotConfiguration::default().summary_scale(AxisScale::Logarithmic);
        let mut group = c.benchmark_group("NIZK_verify_benchmark");
        group.plot_config(plot_config);

        run_nizk_verify_benchmark::<G>(c, s);

        group.finish();
    }
}

// Set the criterion configuration to control the benchmark sample size
fn set_duration() -> Criterion {
    Criterion::default().sample_size(10)
}

// Criterion benchmark group configuration
criterion_group! {
    name = benches_nizk;
    config = set_duration();
    targets = nizk_prove_benchmark::<G1Projective>, nizk_verify_benchmark::<G1Projective>
}

// Criterion main entry point
criterion_main!(benches_nizk);
