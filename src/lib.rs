#![cfg_attr(not(feature = "std"), no_std)]
#![warn(future_incompatible, nonstandard_style, missing_docs)]
#![cfg_attr(not(feature = "blst"), deny(unsafe_code))]
//! ## `poly-multiproof`
//! `poly-multiproof` is a library for generating BDFG21-esque proofs against standard KZG
//! commitments. It supports both the method 1 and method 2 proofs from the paper, and contains
//! an assembly-optimized implementation of method 1 for the BLS12-381 curve. When the points that
//! will be committed to are known beforehand, separate `precompute` modules can be used which
//! pre-compute lagrange polynomials and vanishing polynomials for the points, which can speed up
//! proof generation by a significant amount, especially for larger proof sizes.
//!
//! ### Features
//! * `blst` enables a specific `bls12-381` implementation which uses `blst` for curve msm.
//! * `parallel` enables parallel computation for
//!   * PMP setup generation
//!   * operations in the `data_availability_grid` example
//! * `print-trace` enables some tracing that shows the time certain things take to execute
//!
//! See [the `poly-multiproof` documentation](https://docs.rs/poly-multiproof) for more details.
//!
//! ### Examples
//!
//! An example of using pmp for a grid data availability scheme with 1d erasure encoding is in `examples/data_availability_grid.rs`. To run it with a nice timer, do
//! ```bash
//! cargo run --example data_availability_grid --release --features print-trace,blst,parallel
//! ```
//!
//! which will print out something like this example for a 256x256 grid
//! ```bash
//! Start:   create pmp
//! End:     create pmp ........................................330.880ms
//! Start:   create grid
//! ··Start:   erasure encoding columns
//! ··End:     erasure encoding columns ........................32.143ms
//! ··Start:   computing polynomials from evals
//! ··End:     computing polynomials from evals ................22.562ms
//! ··Start:   computing commitments
//! ··End:     computing commitments ...........................707.127ms
//! End:     create grid .......................................771.776ms
//! Start:   opening to grid
//! End:     opening to grid ...................................244.446ms
//! Start:   verifying grid
//! End:     verifying grid ....................................223.557ms
//! ```
//!
//! There are nice constants in the top of the file to play with.
//!
//! ### Benchmarks
//!
//! To run benchmarks with `arkworks-rs` asm optimizations on x86 machines, run
//! ```bash
//! RUSTFLAGS="-C target-feature=+bmi2,+adx" cargo +nightly criterion --features asm
//! ```
//! or to run with the goal of plotting, run
//! ```bash
//! RUSTFLAGS="-C target-feature=+bmi2,+adx" cargo +nightly criterion --features asm --plotting-backend disabled -- --quick --quiet &> bench_out.txt
//! ```
//! The logs in `bench_out.txt` can then be parsed and plotted in `Plot Benches.ipynb`.
//! Using `--quick` is nice since there are many many inputs benchmarked and it will still take an hour or so to run with `--quick`.
//!
use ark_ec::{scalar_mul::fixed_base::FixedBase, CurveGroup, ScalarMul};
use ark_ff::{Field, PrimeField};
use ark_poly::{
    univariate::{DenseOrSparsePolynomial, DensePolynomial},
    DenseUVPolynomial, EvaluationDomain, GeneralEvaluationDomain,
};
use ark_serialize::{CanonicalDeserialize, CanonicalSerialize, Compress, SerializationError};
use ark_std::{vec, vec::Vec};
use merlin::Transcript;
#[cfg(test)]
use rand::thread_rng as test_rng;

// Public uses
pub use ark_ff;
pub use ark_poly;
pub use ark_serialize;
pub use ark_ec::pairing::Pairing;
pub use merlin;

pub mod m1_cycl;
pub mod method1;
pub mod method2;

pub mod kzg;

pub(crate) mod lagrange;

pub mod traits;

pub mod poly_ops;
pub mod utils;

pub mod msm;

#[cfg(test)]
pub mod testing;

#[cfg(feature = "ark-bls12-381")]
pub use ark_bls12_381;

/// Crate error type
#[derive(Debug, Eq, PartialEq)]
#[cfg_attr(feature = "std", derive(thiserror::Error))]
pub enum Error {
    /// Too many scalars were given. This normally happens when you initialize a `PMP` with too few
    /// points for the degree of the polyomial you want to commit/open to.
    #[cfg_attr(feature = "std", error("Polynomial given is too large"))]
    TooManyScalars {
        /// The number of scalars given
        n_coeffs: usize,
        /// The maximum number of scalars expected
        expected_max: usize,
    },
    /// Found a zero divisor when computing a polynomial quotient
    #[cfg_attr(feature = "std", error("A divisor was zero"))]
    DivisorIsZero,
    /// No polynomials were given
    #[cfg_attr(feature = "std", error("Expected polynomials, none were given"))]
    NoPolynomialsGiven,
    /// The evaluations given to a method did not match with the expected number for the size of
    /// the polynomial.
    #[cfg_attr(feature = "std", error("Given evaluations were the incorrect size"))]
    EvalsIncorrectSize {
        /// The index of the polynomial that had incorrect evals
        poly: usize,
        /// The number of evals
        n_evals: usize,
        /// The expected number of evals
        expected: usize,
    },
    /// Error serializing
    #[cfg_attr(feature = "std", error("Serialization error"))]
    SerializationError,
    /// No points were given
    #[cfg_attr(feature = "std", error("Not given any points"))]
    NoPointsGiven,
    /// Evals and polynomials had different sizes
    #[cfg_attr(
        feature = "std",
        error("Given {n_eval_rows} evaluations, but {n_polys} polynomials")
    )]
    EvalsAndPolysDifferentSizes {
        /// The number of rows of evals
        n_eval_rows: usize,
        /// The number of polynomials
        n_polys: usize,
    },
    /// Evals and points had different sizes
    #[cfg_attr(feature = "std", error("Given {n_points} points, but {n_evals} evals"))]
    EvalsAndPointsDifferentSizes {
        /// The number of points
        n_points: usize,
        /// The number of eval rows
        n_evals: usize,
    },
    /// Evals and commits had different sizes
    #[cfg_attr(
        feature = "std",
        error("Given {n_commits} commits, but {n_evals} evals")
    )]
    EvalsAndCommitsDifferentSizes {
        /// The number of eval rows
        n_evals: usize,
        /// The number of commits
        n_commits: usize,
    },
    /// Failed to construct a domain of the given size
    #[cfg_attr(feature = "std", error("Unable to construct a domain of size {0}"))]
    DomainConstructionFailed(usize),
    /// Subgroup index was invalid
    #[cfg_attr(
        feature = "std",
        error("Invalid subgroup index {idx} for {n_splits} splits")
    )]
    InvalidSubgroupIndex {
        /// Index
        idx: usize,
        /// Number of splits
        n_splits: usize,
    },
    /// Invalid input length
    #[cfg_attr(
        feature = "std",
        error("Invalid input length: {got}, expected {expected}")
    )]
    InvalidInputLength {
        /// Expected length
        expected: usize,
        /// Actual length
        got: usize,
    },
}

impl From<SerializationError> for Error {
    fn from(_: SerializationError) -> Self {
        Self::SerializationError
    }
}

/// A KZG commitment, consisting of a single G1 group element
#[derive(Debug, Clone, CanonicalSerialize, CanonicalDeserialize)]
pub struct Commitment<E: Pairing>(pub E::G1Affine);

impl<E: Pairing> Commitment<E> {
    /// Given a set of commitments and a target output size that is a power of 2,
    /// extend the commitments to the target size using FFTs.
    pub fn extend_commitments(
        commits: impl AsRef<[Commitment<E>]>,
        output_size: usize,
    ) -> Result<Vec<Self>, Error> {
        let mut vals: Vec<E::G1> = commits
            .as_ref()
            .iter()
            .map(|x| x.0.into())
            .collect::<Vec<_>>();
        let domain = GeneralEvaluationDomain::<E::ScalarField>::new(vals.len())
            .ok_or(Error::DomainConstructionFailed(vals.len()))?;
        let domain_ext = GeneralEvaluationDomain::<E::ScalarField>::new(output_size)
            .ok_or(Error::DomainConstructionFailed(output_size))?;
        domain.ifft_in_place(&mut vals);
        domain_ext.fft_in_place(&mut vals);
        Ok(vals
            .into_iter()
            .map(|x| Commitment(x.into()))
            .collect::<Vec<_>>())
    }
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
        return Err(Error::TooManyScalars {
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

pub(crate) fn gen_curve_powers_proj<G: ScalarMul + CurveGroup>(
    powers: &[G::ScalarField],
    base: G,
) -> Vec<G> {
    let window_size = FixedBase::get_mul_window_size(powers.len());
    let scalar_size = G::ScalarField::MODULUS_BIT_SIZE as usize;
    let g_table = FixedBase::get_window_table::<G>(scalar_size, window_size, base);
    FixedBase::msm::<G>(scalar_size, window_size, &g_table, powers)
}

pub(crate) fn gen_curve_powers<G: ScalarMul + CurveGroup>(
    powers: &[G::ScalarField],
    base: G,
) -> Vec<G::Affine> {
    G::normalize_batch(&gen_curve_powers_proj(powers, base))
}

pub(crate) fn get_field_size<F: Field + CanonicalSerialize>() -> usize {
    F::zero().serialized_size(Compress::Yes)
}

pub(crate) fn transcribe_points_and_evals<F: CanonicalSerialize>(
    transcript: &mut Transcript,
    points: &[F],
    evals: &[impl AsRef<[F]>],
    field_size_bytes: usize,
) -> Result<(), Error> {
    let n_points = points.len();
    let mut eval_bytes = vec![0u8; field_size_bytes * n_points * evals.len()];
    for (i, e) in evals.iter().enumerate() {
        if e.as_ref().len() != n_points {
            return Err(Error::EvalsIncorrectSize {
                poly: i,
                n_evals: e.as_ref().len(),
                expected: n_points,
            });
        }
        for (j, p) in e.as_ref().iter().enumerate() {
            let start = (i * n_points + j) * field_size_bytes;
            p.serialize_compressed(&mut eval_bytes[start..start + field_size_bytes])?;
        }
    }
    transcript.append_message(b"open evals", &eval_bytes);
    let mut point_bytes = vec![0u8; field_size_bytes * n_points];
    for (i, p) in points.iter().enumerate() {
        p.serialize_compressed(&mut point_bytes[i * field_size_bytes..(i + 1) * field_size_bytes])?;
    }
    transcript.append_message(b"open points", &point_bytes);
    Ok(())
}

pub(crate) fn transcribe_generic<F: CanonicalSerialize>(
    transcript: &mut Transcript,
    label: &'static [u8],
    f: &F,
) -> Result<(), Error> {
    let elt_size = f.serialized_size(Compress::Yes);
    let mut buf = vec![0u8; elt_size];
    f.serialize_compressed(&mut buf)?;
    transcript.append_message(label, &buf);
    Ok(())
}

pub(crate) fn get_challenge<F: PrimeField>(
    transcript: &mut Transcript,
    label: &'static [u8],
    field_size_bytes: usize,
) -> F {
    let mut challenge_bytes = vec![0u8; field_size_bytes];
    transcript.challenge_bytes(label, &mut challenge_bytes);
    F::from_be_bytes_mod_order(&challenge_bytes)
}

pub(crate) fn check_opening_sizes<F>(
    evals: &[impl AsRef<[F]>],
    polys: &[impl AsRef<[F]>],
    n_points: usize,
) -> Result<(), Error> {
    if evals.len() != polys.len() {
        return Err(Error::EvalsAndPolysDifferentSizes {
            n_eval_rows: evals.len(),
            n_polys: polys.len(),
        });
    }
    for e in evals {
        if e.as_ref().len() != n_points {
            return Err(Error::EvalsAndPointsDifferentSizes {
                n_evals: e.as_ref().len(),
                n_points,
            });
        }
    }
    Ok(())
}

pub(crate) fn check_verify_sizes<F, C>(
    commits: &[C],
    evals: &[impl AsRef<[F]>],
    n_points: usize,
) -> Result<(), Error> {
    if evals.len() != commits.len() {
        return Err(Error::EvalsAndCommitsDifferentSizes {
            n_evals: evals.len(),
            n_commits: commits.len(),
        });
    }
    for e in evals {
        if e.as_ref().len() != n_points {
            return Err(Error::EvalsAndPointsDifferentSizes {
                n_evals: e.as_ref().len(),
                n_points,
            });
        }
    }
    Ok(())
}

/// This macro is used to iterate over a slice in parallel if the `parallel` feature is enabled.
#[macro_export]
macro_rules! cfg_iter {
    ($e: expr) => {{
        #[cfg(feature = "parallel")]
        let result = $e.par_iter().enumerate();

        #[cfg(not(feature = "parallel"))]
        let result = $e.iter().enumerate();

        result
    }};
}
