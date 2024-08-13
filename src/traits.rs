//! Traits used in the BDFG21 and KZG Schemes
use ark_ec::pairing::{Pairing, PairingOutput};
use ark_serialize::{CanonicalDeserialize, CanonicalSerialize};
use ark_std::vec::Vec;
use merlin::Transcript;

use crate::{Commitment, Error};

/// A curve-agnostic trait for a KZG commitment scheme
pub trait Committer<E: Pairing> {
    /// Commit to the given polynomial
    fn commit(&self, poly: impl AsRef<[E::ScalarField]>) -> Result<Commitment<E>, Error>;
}

/// A curve-agnostic trait for making KZG opening proofs
pub trait KZGProof<E: Pairing>: Sized {
    /// The output proof type
    type Proof: Clone;

    /// Compute the witness polynomial for a given point. This can be passed into `open` to create
    /// a proof.
    fn compute_witness_polynomial(
        &self,
        poly: Vec<E::ScalarField>,
        point: E::ScalarField,
    ) -> Result<Vec<E::ScalarField>, Error>;

    /// Creates a proof from a witness polynomial. This essentially commits to the witness
    /// polynomial
    fn open(&self, witness_poly: Vec<E::ScalarField>) -> Result<Self::Proof, Error>;

    /// Verifies a proof against a commitment
    fn verify<M: MSMEngine<E = E>>(
        &self,
        commit: &crate::Commitment<E>,
        point: E::ScalarField,
        value: E::ScalarField,
        proof: &Self::Proof,
    ) -> Result<bool, crate::Error>;
}

/// A curve-agnostic trait for a BDFG commitment scheme *with precomputation*
pub trait PolyMultiProof<E: Pairing>: Sized {
    /// The output proof type
    type Proof: Clone;

    /// Creates a of the given polynomials at the given point set index
    fn open(
        &self,
        transcript: &mut Transcript,
        evals: &[impl AsRef<[E::ScalarField]>],
        polys: &[impl AsRef<[E::ScalarField]>],
        point_set_index: usize,
    ) -> Result<Self::Proof, Error>;

    /// Verifies a proof against the given set of commitments and points
    fn verify(
        &self,
        transcript: &mut Transcript,
        commits: &[Commitment<E>],
        point_set_index: usize,
        evals: &[impl AsRef<[E::ScalarField]>],
        proof: &Self::Proof,
    ) -> Result<bool, Error>;
}

/// A curve-agnostic trait for a BDFG commitment scheme *without precomputation*
pub trait PolyMultiProofNoPrecomp<E: Pairing>: Sized {
    /// The output proof type
    type Proof: Clone;

    /// Creates a proof of the given polynomials and evals at the given points
    fn open(
        &self,
        transcript: &mut Transcript,
        evals: &[impl AsRef<[E::ScalarField]>],
        polys: &[impl AsRef<[E::ScalarField]>],
        points: &[E::ScalarField],
    ) -> Result<Self::Proof, Error>;

    /// Verifies a proof against the given set of commitments and points
    fn verify(
        &self,
        transcript: &mut Transcript,
        commits: &[Commitment<E>],
        points: &[E::ScalarField],
        evals: &[impl AsRef<[E::ScalarField]>],
        proof: &Self::Proof,
    ) -> Result<bool, Error>;
}

/// A curve-agnostic trait for fast multi-scalar multiplication
pub trait MSMEngine: Clone + Copy {
    /// The curve type implemented
    type E: Pairing;
    /// The prepared G1 Scalars
    type G1Prepared: Clone;
    /// The prepared G2 Scalars
    type G2Prepared: Clone;

    /// Prepare the given points for multi-scalar multiplication
    fn prepare_g1(g: Vec<<Self::E as Pairing>::G1Affine>) -> Self::G1Prepared;

    /// Prepare the given points for multi-scalar multiplication
    fn prepare_g2(g: Vec<<Self::E as Pairing>::G2Affine>) -> Self::G2Prepared;

    /// Perform a multi-scalar multiplication on the given G1 elements
    fn multi_scalar_mul_g1(
        g: &Self::G1Prepared,
        s: impl AsRef<[<Self::E as Pairing>::ScalarField]>,
    ) -> Result<<Self::E as Pairing>::G1, Error>;

    /// Perform a multi-scalar multiplication on the given G1 elements
    fn multi_scalar_mul_g2(
        g: &Self::G2Prepared,
        s: impl AsRef<[<Self::E as Pairing>::ScalarField]>,
    ) -> Result<<Self::E as Pairing>::G2, Error>;

    /// Checks that e(p1, q1) == e(p2, q2)
    fn pairing_eq_check(
        p1: <Self::E as Pairing>::G1Affine,
        q1: <Self::E as Pairing>::G2Affine,
        p2: <Self::E as Pairing>::G1Affine,
        q2: <Self::E as Pairing>::G2Affine,
    ) -> bool;

    /// Computes e(p1, q1)
    fn pairing(
        p1: <Self::E as Pairing>::G1Affine,
        q1: <Self::E as Pairing>::G2Affine,
    ) -> PairingOutput<Self::E>;
}

/// Utility trait for serialization and deserialization
pub trait AsBytes<const N: usize>: Sized {
    /// Convert to bytes
    fn to_bytes(&self) -> Result<[u8; N], Error>;
    /// Convert from bytes
    fn from_bytes(bytes: &[u8; N]) -> Result<Self, Error>;
}

#[cfg(feature = "ark-bls12-381")]
impl AsBytes<48> for Commitment<ark_bls12_381::Bls12_381> {
    fn to_bytes(&self) -> Result<[u8; 48], Error> {
        let mut out = [0u8; 48];
        self.0.serialize_compressed(&mut out[..])?;
        Ok(out)
    }

    fn from_bytes(bytes: &[u8; 48]) -> Result<Self, Error> {
        Ok(Self(ark_bls12_381::G1Affine::deserialize_compressed(
            &bytes[..],
        )?))
    }
}

#[cfg(feature = "ark-bls12-381")]
impl AsBytes<48> for crate::method1::Proof<ark_bls12_381::Bls12_381> {
    fn to_bytes(&self) -> Result<[u8; 48], Error> {
        let mut out = [0u8; 48];
        self.0.serialize_compressed(&mut out[..])?;
        Ok(out)
    }

    fn from_bytes(bytes: &[u8; 48]) -> Result<Self, Error> {
        Ok(Self(ark_bls12_381::G1Affine::deserialize_compressed(
            &bytes[..],
        )?))
    }
}

#[cfg(feature = "ark-bls12-381")]
impl AsBytes<96> for crate::method2::Proof<ark_bls12_381::Bls12_381> {
    fn to_bytes(&self) -> Result<[u8; 96], Error> {
        let mut out = [0u8; 96];
        self.0.serialize_compressed(&mut out[..48])?;
        self.1.serialize_compressed(&mut out[48..])?;
        Ok(out)
    }

    fn from_bytes(bytes: &[u8; 96]) -> Result<Self, Error> {
        Ok(Self(
            ark_bls12_381::G1Affine::deserialize_compressed(&bytes[..48])?,
            ark_bls12_381::G1Affine::deserialize_compressed(&bytes[48..])?,
        ))
    }
}

#[cfg(feature = "ark-bls12-381")]
impl AsBytes<32> for ark_bls12_381::Fr {
    fn to_bytes(&self) -> Result<[u8; 32], Error> {
        let mut out = [0u8; 32];
        self.serialize_compressed(&mut out[..])?;
        Ok(out)
    }

    fn from_bytes(bytes: &[u8; 32]) -> Result<Self, Error> {
        Ok(Self::deserialize_compressed(&bytes[..])?)
    }
}

#[cfg(not(feature = "ark-bls12-381"))]
impl<T: CanonicalDeserialize + CanonicalSerialize, const N: usize> AsBytes<N> for T {
    fn to_bytes(&self) -> Result<[u8; N], Error> {
        let mut out = [0u8; N];
        self.serialize_compressed(&mut out[..])?;
        Ok(out)
    }

    fn from_bytes(bytes: &[u8; N]) -> Result<Self, Error> {
        Ok(Self::deserialize_compressed(&bytes[..])?)
    }
}
