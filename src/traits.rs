use ark_ec::pairing::Pairing;
use ark_serialize::{CanonicalDeserialize, CanonicalSerialize};
use ark_std::vec::Vec;
use merlin::Transcript;

use crate::{Commitment, Error};

pub trait Committer<E: Pairing> {
    fn commit(&self, poly: impl AsRef<[E::ScalarField]>) -> Result<Commitment<E>, Error>;
}

pub trait KZGProof<E: Pairing>: Sized {
    type Proof: Clone;

    fn compute_witness_polynomial(
        &self,
        poly: Vec<E::ScalarField>,
        point: E::ScalarField,
    ) -> Result<Vec<E::ScalarField>, Error>;

    fn open(&self, witness_poly: Vec<E::ScalarField>) -> Result<Self::Proof, Error>;

    fn verify(
        &self,
        commit: &crate::Commitment<E>,
        point: E::ScalarField,
        value: E::ScalarField,
        proof: &Self::Proof,
    ) -> Result<bool, crate::Error>;
}

pub trait PolyMultiProof<E: Pairing>: Sized {
    type Proof: Clone;

    fn open(
        &self,
        transcript: &mut Transcript,
        evals: &[impl AsRef<[E::ScalarField]>],
        polys: &[impl AsRef<[E::ScalarField]>],
        point_set_index: usize,
    ) -> Result<Self::Proof, Error>;

    fn verify(
        &self,
        transcript: &mut Transcript,
        commits: &[Commitment<E>],
        point_set_index: usize,
        evals: &[impl AsRef<[E::ScalarField]>],
        proof: &Self::Proof,
    ) -> Result<bool, Error>;
}

pub trait PolyMultiProofNoPrecomp<E: Pairing>: Sized {
    type Proof: Clone;

    fn open(
        &self,
        transcript: &mut Transcript,
        evals: &[impl AsRef<[E::ScalarField]>],
        polys: &[impl AsRef<[E::ScalarField]>],
        points: &[E::ScalarField],
    ) -> Result<Self::Proof, Error>;

    fn verify(
        &self,
        transcript: &mut Transcript,
        commits: &[Commitment<E>],
        points: &[E::ScalarField],
        evals: &[impl AsRef<[E::ScalarField]>],
        proof: &Self::Proof,
    ) -> Result<bool, Error>;
}

pub trait AsBytes<const N: usize>: Sized {
    fn to_bytes(&self) -> Result<[u8; N], Error>;
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
