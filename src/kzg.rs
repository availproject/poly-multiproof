//! KZG10 proof system for each instance of the multiproof scheme
use core::ops::Mul;

use ark_ec::{pairing::Pairing, AffineRepr};
use ark_poly::{univariate::DensePolynomial, DenseUVPolynomial};

use crate::{
    m1_cycl::M1CyclPrecomp,
    method1::{precompute::M1Precomp, M1NoPrecomp},
    poly_div_q_r,
    traits::{Committer, KZGProof, MSMEngine},
    Error,
};
use ark_std::vec;
use ark_std::vec::Vec;

impl<E: Pairing, A: Committer<E> + WithSrs<E>> KZGProof<E> for A {
    type Proof = crate::method1::Proof<E>;

    fn compute_witness_polynomial(
        &self,
        poly: Vec<<E as Pairing>::ScalarField>,
        point: <E as Pairing>::ScalarField,
    ) -> Result<Vec<<E as Pairing>::ScalarField>, Error> {
        use ark_ff::One;
        let poly = DensePolynomial::from_coefficients_vec(poly);
        let divisor = DensePolynomial::from_coefficients_vec(vec![-point, E::ScalarField::one()]);
        let (q, _r) = poly_div_q_r(poly.into(), divisor.into())?;
        Ok(q)
    }

    fn open(&self, witness_poly: Vec<<E as Pairing>::ScalarField>) -> Result<Self::Proof, Error> {
        self.commit(witness_poly)
            .map(|commit| Self::Proof { 0: commit.0 })
    }

    fn verify<M: MSMEngine<E = E>>(
        &self,
        commit: &crate::Commitment<E>,
        point: <E as Pairing>::ScalarField,
        value: <E as Pairing>::ScalarField,
        proof: &Self::Proof,
    ) -> Result<bool, crate::Error> {
        let g1 = self.g1s()[0];
        let g2 = self.g2s()[0];
        let lhsg1 = commit.0.into_group() - g1.mul(value);
        let lhsg2 = g2;

        let rhsg1 = proof.0;
        let rhsg2 = self.g2s()[1].into_group() - g2.mul(point);

        Ok(M::pairing_eq_check(
            lhsg1.into(),
            lhsg2,
            rhsg1,
            rhsg2.into(),
        ))
    }
}

// Impls
trait WithSrs<E: Pairing> {
    fn g1s(&self) -> &[E::G1Affine];
    fn g2s(&self) -> &[E::G2Affine];
}

impl<E: Pairing, M: MSMEngine<E = E>> WithSrs<E> for M1NoPrecomp<E, M> {
    fn g1s(&self) -> &[E::G1Affine] {
        &self.powers_of_g1
    }

    fn g2s(&self) -> &[E::G2Affine] {
        &self.powers_of_g2
    }
}

impl<E: Pairing, M: MSMEngine<E = E>> WithSrs<E> for M1Precomp<E, M> {
    fn g1s(&self) -> &[E::G1Affine] {
        &self.inner.powers_of_g1
    }

    fn g2s(&self) -> &[E::G2Affine] {
        &self.inner.powers_of_g2
    }
}

impl<E: Pairing, M: MSMEngine<E = E>> WithSrs<E> for M1CyclPrecomp<E, M> {
    fn g1s(&self) -> &[E::G1Affine] {
        &self.inner.powers_of_g1
    }

    fn g2s(&self) -> &[E::G2Affine] {
        &self.inner.powers_of_g2
    }
}
