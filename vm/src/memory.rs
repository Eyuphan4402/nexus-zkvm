//! Virtual Machine Memory

pub mod cacheline;
pub mod jolt;
pub mod paged;
pub mod path;
pub mod trie;

use ark_r1cs_std::fields::fp::FpVar;
use ark_relations::r1cs::{ConstraintSystemRef, SynthesisError};
use ark_serialize::{CanonicalDeserialize, CanonicalSerialize};

use crate::circuit::F;
use crate::error::Result;
use crate::rv32::{LOP, SOP};
use cacheline::CacheLine;

/// A `MemoryProof` implementation provides the memory commitments and
/// in-circuit proofs of the memory commitments. Each memory controller
/// must provide proof values that implement this trait.
pub trait MemoryProof: Default + Clone + CanonicalSerialize + CanonicalDeserialize {
    /// Type of in-circuit parameters needed, if any. For instance, hash parameters.
    type Params;

    /// Generate required parameters within the given constraint system.
    fn params(cs: ConstraintSystemRef<F>) -> Result<Self::Params, SynthesisError>;

    /// Generate in-circuit verification of the memory proof.
    fn circuit(
        &self,
        cs: ConstraintSystemRef<F>,
        params: &Self::Params,
        root: &FpVar<F>,
        data: &[FpVar<F>],
    ) -> Result<(), SynthesisError>;

    /// Return the memory commitment related to this proof.
    fn commit(&self) -> F;

    /// Return the `CacheLine` data related to this proof.
    fn data(&self) -> [F; 2];
}

/// A `Memory` implementation is responsible for managing the machine's
/// memory, and providing access to `CacheLine`s. Each implementation
/// has associated commitment values and circuits for checking the
/// commitments to the memory contents. A `Memory` implementation and
/// the associated `MemoryProof`s are referred to as a "memory controller".
pub trait Memory: Default {
    /// Type of the memory commitment proof values generated by
    /// this memory controller.
    type Proof: MemoryProof;

    /// Query the cacheline at `addr`
    fn query(&self, addr: u32) -> (CacheLine, Self::Proof);

    /// Update the cacheline at `addr` using the function `f`.
    fn update<F>(&mut self, addr: u32, f: F) -> Result<Self::Proof>
    where
        F: Fn(&mut CacheLine) -> Result<()>;

    /// read instruction at address
    fn read_inst(&self, addr: u32) -> Result<(u32, Self::Proof)> {
        let (cl, path) = self.query(addr);
        Ok((cl.lw(addr)?, path))
    }

    /// write instruction at address
    fn write_inst(&mut self, addr: u32, val: u32) -> Result<()> {
        let _ = self.update(addr, |cl| cl.sw(addr, val))?;
        Ok(())
    }

    /// perform load according to `lop`
    fn load(&self, lop: LOP, addr: u32) -> Result<(u32, Self::Proof)> {
        let (cl, path) = self.query(addr);
        Ok((cl.load(lop, addr)?, path))
    }

    /// Load n bytes from an address
    fn load_n(&self, address: u32, len: u32) -> Result<Vec<u8>> {
        (address..address + len)
            .map(|addr| self.load(LOP::LBU, addr).map(|b| b.0 as u8))
            .collect()
    }

    /// perform store according to `sop`
    fn store(&mut self, sop: SOP, addr: u32, val: u32) -> Result<Self::Proof> {
        self.update(addr, |cl| cl.store(sop, addr, val))
    }
}

#[cfg(test)]
mod test {
    use super::{paged::Paged, trie::MerkleTrie, *};
    use crate::rv32::{LOP::*, SOP::*};

    #[test]
    fn test_mem_merkle() {
        test_mem(MerkleTrie::default());
    }

    #[test]
    fn test_mem_paged() {
        test_mem(Paged::default());
    }

    fn test_mem(mut mem: impl Memory) {
        // read before write
        assert_eq!(mem.load(LW, 0x1000).unwrap().0, 0);

        mem.store(SW, 0x1000, 1).unwrap();
        mem.store(SB, 0x1100, 1).unwrap();
        mem.store(SB, 0x1101, 2).unwrap();
        mem.store(SB, 0x1103, 3).unwrap();
        mem.store(SB, 0x1104, 4).unwrap();
        mem.store(SB, 0x11000, 1).unwrap();

        assert_eq!(mem.load(LBU, 0x10ff).unwrap().0, 0);
        assert_eq!(mem.load(LBU, 0x1100).unwrap().0, 1);
        assert_eq!(mem.load(LBU, 0x1101).unwrap().0, 2);
        assert_eq!(mem.load(LBU, 0x1103).unwrap().0, 3);
        assert_eq!(mem.load(LBU, 0x1104).unwrap().0, 4);
        assert_eq!(mem.load(LBU, 0x1105).unwrap().0, 0);
        assert_eq!(mem.load(LBU, 0x11000).unwrap().0, 1);
        assert_eq!(mem.load(LBU, 0x11001).unwrap().0, 0);

        mem.store(SH, 0x1100, 0x708).unwrap();
        assert_eq!(mem.load(LBU, 0x1100).unwrap().0, 8);
        assert_eq!(mem.load(LBU, 0x1101).unwrap().0, 7);
        assert_eq!(mem.load(LHU, 0x1100).unwrap().0, 0x708);
        assert_eq!(mem.load(LHU, 0x1200).unwrap().0, 0);

        mem.store(SW, 0x1200, 0x10203040).unwrap();
        assert_eq!(mem.load(LBU, 0x1200).unwrap().0, 0x40);
        assert_eq!(mem.load(LBU, 0x1201).unwrap().0, 0x30);
        assert_eq!(mem.load(LBU, 0x1202).unwrap().0, 0x20);
        assert_eq!(mem.load(LBU, 0x1203).unwrap().0, 0x10);
        assert_eq!(mem.load(LHU, 0x1200).unwrap().0, 0x3040);
        assert_eq!(mem.load(LHU, 0x1202).unwrap().0, 0x1020);
        assert_eq!(mem.load(LW, 0x1200).unwrap().0, 0x10203040);

        mem.store(SH, 0x1300, 0x81).unwrap();
        assert_eq!(mem.load(LB, 0x1300).unwrap().0, 0xffffff81);

        mem.store(SH, 0x1300, 0x8321).unwrap();
        assert_eq!(mem.load(LH, 0x1300).unwrap().0, 0xffff8321);
    }
}
