//! # BDFG Method 2 with precomputation
use ark_ec::pairing::Pairing;
use ark_poly::univariate::DensePolynomial;
use ark_std::vec::Vec;
use merlin::Transcript;

#[cfg(feature = "parallel")]
use rayon::prelude::*;

use super::{vanishing_polynomial, Error, Proof};
use crate::lagrange::LagrangeInterpContext;
use crate::traits::{Committer, PolyMultiProof};
use crate::{cfg_iter, Commitment};

/// Method 2 with precomputation
#[derive(Clone, Debug)]
pub struct M2Precomp<E: Pairing> {
    /// The inner method 2 object without precomputation
    pub inner: super::M2NoPrecomp<E>,
    point_sets: Vec<Vec<E::ScalarField>>,
    vanishing_polys: Vec<DensePolynomial<E::ScalarField>>,
    lagrange_ctxs: Vec<LagrangeInterpContext<E::ScalarField>>,
}

impl<E: Pairing> M2Precomp<E> {
    /// Make a precompute-optimized version of a method 2 object for the given sets of points
    pub fn from_inner(
        inner: super::M2NoPrecomp<E>,
        point_sets: Vec<Vec<E::ScalarField>>,
    ) -> Result<Self, Error> {
        let vanishing_polys = cfg_iter!(point_sets)
            .map(|(_, ps)| vanishing_polynomial(ps))
            .collect();
        let lagrange_ctxs = cfg_iter!(point_sets)
            .map(|(_, ps)| LagrangeInterpContext::new_from_points(ps.as_ref()))
            .collect::<Result<Vec<_>, Error>>()?;

        Ok(M2Precomp {
            inner,
            point_sets,
            vanishing_polys,
            lagrange_ctxs,
        })
    }
}

impl<E: Pairing> Committer<E> for M2Precomp<E> {
    fn commit(
        &self,
        poly: impl AsRef<[<E as Pairing>::ScalarField]>,
    ) -> Result<Commitment<E>, Error> {
        self.inner.commit(poly)
    }
}

impl<E: Pairing> PolyMultiProof<E> for M2Precomp<E> {
    type Proof = Proof<E>;

    fn open(
        &self,
        transcript: &mut Transcript,
        evals: &[impl AsRef<[E::ScalarField]>],
        polys: &[impl AsRef<[E::ScalarField]>],
        point_set_index: usize,
    ) -> Result<Proof<E>, Error> {
        self.inner.open_with_vanishing_poly(
            transcript,
            evals,
            polys,
            &self.point_sets[point_set_index],
            &self.vanishing_polys[point_set_index],
        )
    }

    fn verify(
        &self,
        transcript: &mut Transcript,
        commits: &[Commitment<E>],
        point_set_index: usize,
        evals: &[impl AsRef<[E::ScalarField]>],
        proof: &Proof<E>,
    ) -> Result<bool, Error> {
        self.inner.verify_with_lag_ctx_vanishing_poly(
            transcript,
            commits,
            &self.point_sets[point_set_index],
            evals,
            proof,
            &self.lagrange_ctxs[point_set_index],
            &self.vanishing_polys[point_set_index],
        )
    }
}

#[cfg(test)]
mod tests {
    use super::M2Precomp;
    use crate::{method2::M2NoPrecomp, test_rng, testing::test_basic_precomp};
    use ark_bls12_381::{Bls12_381, Fr};
    use ark_std::{vec, vec::Vec, UniformRand};

    #[test]
    fn test_basic_open_works() {
        let points = (0..30)
            .map(|_| Fr::rand(&mut test_rng()))
            .collect::<Vec<_>>();
        let s = M2NoPrecomp::<Bls12_381>::new(256, &mut test_rng());
        let s = M2Precomp::<Bls12_381>::from_inner(s, vec![points.clone()])
            .expect("Failed to construct");
        test_basic_precomp(&s, &points)
    }
}
