use crate::{msm::blst::BlstMSMEngine, test_rng, traits::KZGProof, vec, Error, Vec};
use ark_bls12_381::Bls12_381;
use ark_ec::pairing::Pairing;
use ark_ff::{One, UniformRand, Zero};
use ark_poly::{univariate::DensePolynomial, DenseUVPolynomial, Polynomial};
use merlin::Transcript;

use crate::traits::{Committer, PolyMultiProof, PolyMultiProofNoPrecomp};

pub fn test_basic_no_precomp<E: Pairing, P: PolyMultiProofNoPrecomp<E> + Committer<E>>(s: &P) {
    let points = (0..30)
        .map(|_| E::ScalarField::rand(&mut test_rng()))
        .collect::<Vec<_>>();
    let polys = (0..20)
        .map(|_| DensePolynomial::<E::ScalarField>::rand(50, &mut test_rng()))
        .collect::<Vec<_>>();
    let evals: Vec<Vec<_>> = polys
        .iter()
        .map(|p| points.iter().map(|x| p.evaluate(x)).collect())
        .collect();
    let coeffs = polys.iter().map(|p| p.coeffs.clone()).collect::<Vec<_>>();
    let commits = coeffs
        .iter()
        .map(|p| s.commit(p).expect("Commit failed"))
        .collect::<Vec<_>>();
    let mut transcript = Transcript::new(b"testing");
    let open = s
        .open(&mut transcript, &evals, &coeffs, &points)
        .expect("Open failed");
    let mut transcript = Transcript::new(b"testing");
    assert_eq!(
        Ok(true),
        s.verify(&mut transcript, &commits, &points, &evals, &open)
    );
}

/// Basic test for a precomp. Assumes `points` are the zero-th pointset.
pub fn test_basic_precomp<E: Pairing, P: PolyMultiProof<E> + Committer<E>>(
    s: &P,
    points: &Vec<E::ScalarField>,
) {
    let polys = (0..20)
        .map(|_| DensePolynomial::<E::ScalarField>::rand(50, &mut test_rng()))
        .collect::<Vec<_>>();
    let evals: Vec<Vec<_>> = polys
        .iter()
        .map(|p| points.iter().map(|x| p.evaluate(x)).collect())
        .collect();
    let coeffs = polys.iter().map(|p| p.coeffs.clone()).collect::<Vec<_>>();
    let commits = coeffs
        .iter()
        .map(|p| s.commit(p).expect("Commit failed"))
        .collect::<Vec<_>>();
    let mut transcript = Transcript::new(b"testing");
    let open = s
        .open(&mut transcript, &evals, &coeffs, 0)
        .expect("Open failed");
    let mut transcript = Transcript::new(b"testing");
    assert_eq!(
        Ok(true),
        s.verify(&mut transcript, &commits, 0, &evals, &open)
    );
}

pub fn test_size_errors<E: Pairing, P: PolyMultiProofNoPrecomp<E> + Committer<E>>(s: &P) {
    let points = (0..20)
        .map(|_| E::ScalarField::rand(&mut test_rng()))
        .collect::<Vec<_>>();
    let polys = (0..20)
        .map(|_| DensePolynomial::<E::ScalarField>::rand(50, &mut test_rng()))
        .collect::<Vec<_>>();
    let evals: Vec<Vec<_>> = polys
        .iter()
        .map(|p| points.iter().map(|x| p.evaluate(x)).collect())
        .collect();
    let coeffs = polys.iter().map(|p| p.coeffs.clone()).collect::<Vec<_>>();
    let commits = coeffs
        .iter()
        .map(|p| s.commit(p).expect("Commit failed"))
        .collect::<Vec<_>>();
    // Eval row commit mismatch
    assert_eq!(
        Err(Error::EvalsAndPolysDifferentSizes {
            n_eval_rows: 19,
            n_polys: 20,
        }),
        s.open(
            &mut Transcript::new(b"testing"),
            &evals[..19],
            &coeffs,
            &points
        )
        .map(|_| ())
    );
    // Eval poly mismatch
    assert_eq!(
        Err(Error::EvalsAndPolysDifferentSizes {
            n_eval_rows: 20,
            n_polys: 19,
        }),
        s.open(
            &mut Transcript::new(b"testing"),
            &evals,
            &coeffs[..19],
            &points
        )
        .map(|_| ())
    );
    // Evals point mismatch
    assert_eq!(
        Err(Error::EvalsAndPointsDifferentSizes {
            n_evals: 19,
            n_points: 20,
        }),
        s.open(
            &mut Transcript::new(b"testing"),
            &evals.iter().map(|s| &s[..19]).collect::<Vec<_>>(),
            &coeffs,
            &points,
        )
        .map(|_| ())
    );

    let open = s
        .open(&mut Transcript::new(b"testing"), &evals, &coeffs, &points)
        .expect("Open failed");
    // Eval row commit mismatch
    assert_eq!(
        Err(Error::EvalsAndCommitsDifferentSizes {
            n_evals: 20,
            n_commits: 19,
        }),
        s.verify(
            &mut Transcript::new(b"testing"),
            &commits[..19],
            &points,
            &evals,
            &open
        )
    );
    // Eval point mismatch
    assert_eq!(
        Err(Error::EvalsAndPointsDifferentSizes {
            n_evals: 19,
            n_points: 20,
        }),
        s.verify(
            &mut Transcript::new(b"testing"),
            &commits,
            &points,
            &evals.iter().map(|a| &a[..19]).collect::<Vec<_>>(),
            &open
        )
    );
}

pub fn test_kzg(srs: &(impl KZGProof<Bls12_381> + Committer<Bls12_381>)) {
    use ark_bls12_381::Fr;

    fn run(srs: &(impl KZGProof<Bls12_381> + Committer<Bls12_381>), coeffs: Vec<Fr>) {
        let pt = Fr::rand(&mut test_rng());
        let value = DensePolynomial::from_coefficients_vec(coeffs.clone()).evaluate(&pt);

        let commit = srs.commit(&coeffs).unwrap();
        let witness = srs.compute_witness_polynomial(coeffs, pt).unwrap();
        let proof = srs.open(witness).unwrap();

        let veri = srs
            .verify::<BlstMSMEngine>(&commit, pt, value, &proof)
            .unwrap();
        assert!(veri);
    }

    // Test random poly works
    run(srs, DensePolynomial::<Fr>::rand(50, &mut test_rng()).coeffs);
    // Test zero poly works
    run(srs, vec![Fr::zero(); 50]);
    // Test constant poly works
    run(srs, vec![Fr::one(); 50]);
    // Test unit poly works
    let mut unit = vec![Fr::zero(); 50];
    unit[0] = Fr::one();
    run(srs, unit);
}
