use ark_ec::{pairing::Pairing, scalar_mul::fixed_base::FixedBase, CurveGroup, ScalarMul};
use ark_ff::{FftField, Field, PrimeField};
use ark_poly::{
    univariate::{DenseOrSparsePolynomial, DensePolynomial},
    DenseUVPolynomial,
};
use ark_std::ops::{Add, Mul};
use ark_std::rand::RngCore;
#[cfg(test)]
use ark_std::test_rng;

pub mod method1;
pub mod method2;

#[derive(Debug, Eq, PartialEq)]
pub enum Error {
    PolynomialTooLarge {
        n_coeffs: usize,
        expected_max: usize,
    },
    DivisorIsZero,
    NoPolynomialsGiven,
}

pub trait MultiOpenKzg<E: Pairing> {
    type Setup;
    type Commitment;
    type Proof;
    fn new(max_degree: usize, max_pts: usize, rng: &mut impl RngCore) -> Self::Setup;
    fn commit(
        setup: &Self::Setup,
        poly: impl AsRef<[E::ScalarField]>,
    ) -> Result<Self::Commitment, Error>;
    fn open(
        &self,
        polys: &[impl AsRef<[E::ScalarField]>],
        points: &[E::ScalarField],
    ) -> Result<Self::Proof, Error>;
    fn verify(
        &self,
        commits: &[Self::Commitment],
        pts: &[E::ScalarField],
        evals: &[impl AsRef<[E::ScalarField]>],
        proof: &Self::Proof,
    ) -> Result<bool, Error>;
}

pub(crate) fn gen_powers<F: Field>(element: F, len: usize) -> Vec<F> {
    let mut powers = vec![F::one(); len];
    for i in 1..len {
        powers[i] = element * powers[i - 1];
    }
    powers
}

#[inline]
pub(crate) fn curve_msm<G: ScalarMul + CurveGroup>(
    bases: &[G::Affine],
    scalars: &[G::ScalarField],
) -> Result<G, Error> {
    if scalars.len() > bases.len() {
        return Err(Error::PolynomialTooLarge {
            n_coeffs: scalars.len(),
            expected_max: bases.len(),
        });
    }
    let scalars = scalars.iter().map(|x| x.into_bigint()).collect::<Vec<_>>();
    let sp = G::msm_bigint(&bases[..scalars.len()], &scalars);
    Ok(sp)
}

pub(crate) fn vanishing_polynomial<F: Field>(points: impl AsRef<[F]>) -> DensePolynomial<F> {
    let one = DensePolynomial::from_coefficients_vec(vec![F::one()]);
    points
        .as_ref()
        .iter()
        .map(|&point| DensePolynomial::from_coefficients_vec(vec![-point, F::one()]))
        .fold(one, |x, y| x.naive_mul(&y))
}

/// Does polynomial division, returning q, r
pub(crate) fn poly_div_q_r<F: Field>(
    num: DenseOrSparsePolynomial<F>,
    denom: DenseOrSparsePolynomial<F>,
) -> Result<(Vec<F>, Vec<F>), Error> {
    if denom.is_zero() {
        return Err(Error::DivisorIsZero);
    }
    let (q, r) = num.divide_with_q_and_r(&denom).expect("Cannot return none");
    Ok((q.coeffs, r.coeffs))
}

pub(crate) fn linear_combination<F: Field>(
    polynomials: &[impl AsRef<[F]>],
    challenges: &[F],
) -> Option<Vec<F>> {
    polynomials
        .as_ref()
        .iter()
        .zip(challenges.iter())
        .map(|(p, &c)| &DensePolynomial::from_coefficients_slice(p.as_ref()) * c)
        .reduce(|x, y| x + y)?
        .coeffs
        .into()
}

pub(crate) fn gen_curve_powers<G: ScalarMul + CurveGroup>(
    powers: &[G::ScalarField],
    rng: &mut impl RngCore,
) -> Vec<G::Affine> {
    let g = G::rand(rng);
    let window_size = FixedBase::get_mul_window_size(powers.len());
    let scalar_size = G::ScalarField::MODULUS_BIT_SIZE as usize;
    let g_table = FixedBase::get_window_table::<G>(scalar_size, window_size, g);
    let powers_of_g_proj = FixedBase::msm::<G>(scalar_size, window_size, &g_table, powers);
    G::normalize_batch(&powers_of_g_proj)
}

/// This computes the inverse of each `j`-th lagrange polynomial,
/// constructed from `points`, evaluated at `points[j]`
///
/// Ex: for `x_2`, it would be
///     `(x_2 - x_1) * (x_2 - x_3) * (x_2 - x_4) * .. * (x_2 - x_n)
///
/// This is used so that when the 2nd lagrange polynomial is construction,
/// we can divide it by this factor so that it evaluates to 1 on `x_2` and
/// zero at all other `x_i`
pub(crate) fn lagrange_poly_inverses<F: Field>(points: &[F]) -> Vec<F> {
    let mut invs = Vec::with_capacity(points.len());
    for (j, x_j) in points.iter().enumerate() {
        let mut prod = F::one();
        for (k, x_k) in points.iter().enumerate() {
            if j == k {
                continue;
            }
            prod *= *x_j - *x_k;
        }
        // TODO: Enforce this constraint
        invs.push(prod.inverse().expect("Point cannot be zero"));
    }
    invs
}

pub(crate) fn gen_lagrange_polynomials<F: FftField>(points: &[F]) -> Vec<DensePolynomial<F>> {
    let mut lang = Vec::new();
    for (j, _x_j) in points.iter().enumerate() {
        let mut l_poly: DensePolynomial<F> = DensePolynomial::from_coefficients_vec(vec![F::one()]);
        for (k, x_k) in points.iter().enumerate() {
            if j == k {
                continue;
            }
            let tmp_poly: DensePolynomial<F> =
                DensePolynomial::from_coefficients_vec(vec![-(*x_k), F::one()]);
            // This does fft mul... not sure if it's actually faster
            l_poly = l_poly.mul(&tmp_poly);
        }
        lang.push(l_poly);
    }
    lang
}

pub(crate) fn do_lagrange_interpolation<F: FftField>(
    evals: &[impl AsRef<[F]>],
    points: &[F],
    inverses: &[F],
    polys: &[DensePolynomial<F>],
) -> Vec<DensePolynomial<F>> {
    evals
        .iter()
        .map(|evals_i| {
            let mut res = DensePolynomial::from_coefficients_vec(vec![F::zero()]);
            for (j, (_x_j, y_j)) in points.iter().zip(evals_i.as_ref().iter()).enumerate() {
                let l_poly = polys[j].mul(inverses[j] * y_j);
                res = (&res).add(&l_poly);
            }
            res
        })
        .collect()
}

pub(crate) fn lagrange_interp<F: FftField>(
    evals: &[impl AsRef<[F]>],
    points: &[F],
) -> Vec<DensePolynomial<F>> {
    let inverses = lagrange_poly_inverses(points);
    let polys = gen_lagrange_polynomials(points);
    do_lagrange_interpolation(evals, points, &inverses, &polys)
}
