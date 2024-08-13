//! Multi-scalar multiplication engines
use std::marker::PhantomData;

use ark_ec::pairing::{Pairing, PairingOutput};

use crate::{traits::MSMEngine, Error};

/// The MSM engine for the Arkworks library
#[derive(Clone, Copy)]
pub struct ArkMSMEngine<E: Pairing>(PhantomData<E>);

impl<E: Pairing> MSMEngine for ArkMSMEngine<E> {
    type E = E;

    type G1Prepared = Vec<<Self::E as Pairing>::G1Affine>;
    type G2Prepared = Vec<<Self::E as Pairing>::G2Affine>;

    fn prepare_g1(g: Vec<<Self::E as Pairing>::G1Affine>) -> Self::G1Prepared {
        g
    }

    fn prepare_g2(g: Vec<<Self::E as Pairing>::G2Affine>) -> Self::G2Prepared {
        g
    }

    fn multi_scalar_mul_g1(
        g: &Self::G1Prepared,
        s: impl AsRef<[<Self::E as Pairing>::ScalarField]>,
    ) -> Result<<Self::E as Pairing>::G1, Error> {
        crate::curve_msm(g, s.as_ref())
    }

    fn multi_scalar_mul_g2(
        g: &Self::G2Prepared,
        s: impl AsRef<[<Self::E as Pairing>::ScalarField]>,
    ) -> Result<<Self::E as Pairing>::G2, Error> {
        crate::curve_msm(g, s.as_ref())
    }

    fn pairing_eq_check(
        a: <Self::E as Pairing>::G1Affine,
        b: <Self::E as Pairing>::G2Affine,
        c: <Self::E as Pairing>::G1Affine,
        d: <Self::E as Pairing>::G2Affine,
    ) -> bool {
        Self::E::pairing(a, b) == Self::E::pairing(c, d)
    }

    fn pairing(
        p1: <Self::E as Pairing>::G1Affine,
        q1: <Self::E as Pairing>::G2Affine,
    ) -> PairingOutput<Self::E> {
        Self::E::pairing(p1, q1)
    }
}

#[cfg(feature = "blst")]
/// The MSM engine for the blst library
pub mod blst {
    use ::blst::{
        blst_final_exp, blst_fp, blst_fp12, blst_fp12_finalverify, blst_fp2, blst_miller_loop,
        blst_p1, blst_p1_affine, blst_p1_mult, blst_p2, blst_p2_affine, blst_p2_mult, MultiPoint,
    };
    use ark_ec::{pairing::PairingOutput, AffineRepr};
    use ark_ff::{BigInt, Zero};
    use ark_serialize::CanonicalSerialize;

    use super::*;

    /// The MSM engine for the blst library
    #[derive(Clone, Copy)]
    pub struct BlstMSMEngine;

    impl MSMEngine for BlstMSMEngine {
        type E = ark_bls12_381::Bls12_381;

        type G1Prepared = PreparedG1s;
        type G2Prepared = PreparedG2s;

        fn prepare_g1(g: Vec<<Self::E as Pairing>::G1Affine>) -> Self::G1Prepared {
            PreparedG1s::from_affines(g)
        }

        fn prepare_g2(g: Vec<<Self::E as Pairing>::G2Affine>) -> Self::G2Prepared {
            PreparedG2s::from_affines(g)
        }

        fn multi_scalar_mul_g1(
            g: &Self::G1Prepared,
            s: impl AsRef<[<Self::E as Pairing>::ScalarField]>,
        ) -> Result<<Self::E as Pairing>::G1, Error> {
            g.msm(s.as_ref())
        }

        fn multi_scalar_mul_g2(
            g: &Self::G2Prepared,
            s: impl AsRef<[<Self::E as Pairing>::ScalarField]>,
        ) -> Result<<Self::E as Pairing>::G2, Error> {
            g.msm(s.as_ref())
        }

        fn pairing_eq_check(
            p1: <Self::E as Pairing>::G1Affine,
            q1: <Self::E as Pairing>::G2Affine,
            p2: <Self::E as Pairing>::G1Affine,
            q2: <Self::E as Pairing>::G2Affine,
        ) -> bool {
            check_pairings_equal(p1, q1, p2, q2)
        }

        fn pairing(
            p1: <Self::E as Pairing>::G1Affine,
            q1: <Self::E as Pairing>::G2Affine,
        ) -> PairingOutput<Self::E> {
            PairingOutput(pairing(p1, q1))
        }
    }

    /// Check that two pairings are equal by doing two miller loops and a single final exponentiation
    // TODO: this should really be checked by someone who understands the blst lib
    pub fn check_pairings_equal(
        p1: ark_bls12_381::G1Affine,
        q1: ark_bls12_381::G2Affine,
        p2: ark_bls12_381::G1Affine,
        q2: ark_bls12_381::G2Affine,
    ) -> bool {
        let bp1 = convert_g1_affine(p1);
        let bq1 = convert_g2_affine(q1);
        let bp2 = convert_g1_affine(p2);
        let bq2 = convert_g2_affine(q2);

        unsafe {
            let mut ret1 = blst::blst_fp12::default();
            blst::blst_miller_loop(
                &mut ret1 as *mut blst::blst_fp12,
                &bq1 as *const blst::blst_p2_affine,
                &bp1 as *const blst::blst_p1_affine,
            );
            let mut ret2 = blst::blst_fp12::default();
            blst::blst_miller_loop(
                &mut ret2 as *mut blst::blst_fp12,
                &bq2 as *const blst::blst_p2_affine,
                &bp2 as *const blst::blst_p1_affine,
            );
            blst::blst_fp12_finalverify(
                &ret1 as *const blst::blst_fp12,
                &ret2 as *const blst::blst_fp12,
            )
        }
    }

    /// Compute a pairing
    pub fn pairing(p: ark_bls12_381::G1Affine, q: ark_bls12_381::G2Affine) -> ark_bls12_381::Fq12 {
        let bp = convert_g1_affine(p);
        let bq = convert_g2_affine(q);
        let ret = unsafe {
            let mut ret1 = blst::blst_fp12::default();
            blst::blst_miller_loop(
                &mut ret1 as *mut blst::blst_fp12,
                &bq as *const blst::blst_p2_affine,
                &bp as *const blst::blst_p1_affine,
            );
            let mut ret2 = blst::blst_fp12::default();
            blst_final_exp(
                &mut ret2 as *mut blst::blst_fp12,
                &mut ret1 as *const blst::blst_fp12,
            );
            ret2
        };
        ark_bls12_381::Fq12 {
            c0: ark_ff::CubicExtField {
                c0: ark_ff::QuadExtField {
                    c0: ark_ff::Fp(BigInt(ret.fp6[0].fp2[0].fp[0].l), PhantomData),
                    c1: ark_ff::Fp(BigInt(ret.fp6[0].fp2[0].fp[1].l), PhantomData),
                },
                c1: ark_ff::QuadExtField {
                    c0: ark_ff::Fp(BigInt(ret.fp6[0].fp2[1].fp[0].l), PhantomData),
                    c1: ark_ff::Fp(BigInt(ret.fp6[0].fp2[1].fp[1].l), PhantomData),
                },
                c2: ark_ff::QuadExtField {
                    c0: ark_ff::Fp(BigInt(ret.fp6[0].fp2[2].fp[0].l), PhantomData),
                    c1: ark_ff::Fp(BigInt(ret.fp6[0].fp2[2].fp[1].l), PhantomData),
                },
            },
            c1: ark_ff::CubicExtField {
                c0: ark_ff::QuadExtField {
                    c0: ark_ff::Fp(BigInt(ret.fp6[1].fp2[0].fp[0].l), PhantomData),
                    c1: ark_ff::Fp(BigInt(ret.fp6[1].fp2[0].fp[1].l), PhantomData),
                },
                c1: ark_ff::QuadExtField {
                    c0: ark_ff::Fp(BigInt(ret.fp6[1].fp2[1].fp[0].l), PhantomData),
                    c1: ark_ff::Fp(BigInt(ret.fp6[1].fp2[1].fp[1].l), PhantomData),
                },
                c2: ark_ff::QuadExtField {
                    c0: ark_ff::Fp(BigInt(ret.fp6[1].fp2[2].fp[0].l), PhantomData),
                    c1: ark_ff::Fp(BigInt(ret.fp6[1].fp2[2].fp[1].l), PhantomData),
                },
            },
        }
    }

    /// Prepared G1 affines for fast MSM
    #[derive(Clone)]
    pub struct PreparedG1s {
        first: Option<blst_p1>,
        all: Vec<blst_p1_affine>,
        len: usize,
    }

    impl PreparedG1s {
        /// Construct from a set of affine points
        fn from_affines(value: Vec<ark_bls12_381::G1Affine>) -> Self {
            let len = value.len();
            let first = value.get(0).map(|p1| convert_g1(p1.into_group()));
            let all: Vec<_> = value.into_iter().map(convert_g1_affine).collect();
            Self { first, all, len }
        }

        /// Perform the MSM
        fn msm(
            &self,
            mut scalars: &[ark_bls12_381::Fr],
        ) -> Result<ark_bls12_381::G1Projective, Error> {
            scalars = trim_zeros(scalars);
            check_scalars(scalars, self.len)?;
            if scalars.len() == 0 || self.len == 0 {
                return Ok(Zero::zero());
            }
            let scalars_le = prep_scalars(scalars);
            let res_p1 = if scalars.len() == 1 {
                let mut out = blst_p1::default();
                unsafe {
                    blst_p1_mult(
                        &mut out as *mut blst_p1,
                        &self.first.expect("len != 0"),
                        scalars_le.as_ptr(),
                        255,
                    )
                }
                out
            } else {
                let a: &[blst_p1_affine] = &self.all[..scalars.len()];
                a.mult(&scalars_le, 255)
            };
            Ok(ark_bls12_381::G1Projective {
                x: ark_ff::Fp(BigInt(res_p1.x.l), PhantomData),
                y: ark_ff::Fp(BigInt(res_p1.y.l), PhantomData),
                z: ark_ff::Fp(BigInt(res_p1.z.l), PhantomData),
            })
        }
    }

    /// Prepared G2 points for fast MSM
    #[derive(Clone)]
    pub struct PreparedG2s {
        first: Option<blst_p2>,
        all: Vec<blst_p2_affine>,
        len: usize,
    }

    impl PreparedG2s {
        fn from_affines(value: Vec<ark_bls12_381::G2Affine>) -> Self {
            let len = value.len();
            let first = value.get(0).map(|p2| convert_g2(p2.into_group()));
            let all: Vec<_> = value.into_iter().map(convert_g2_affine).collect();
            Self { first, all, len }
        }

        fn msm(
            &self,
            mut scalars: &[ark_bls12_381::Fr],
        ) -> Result<ark_bls12_381::G2Projective, Error> {
            scalars = trim_zeros(scalars);
            check_scalars(scalars, self.len)?;
            if scalars.len() == 0 || self.len == 0 {
                return Ok(Zero::zero());
            }
            let scalars_le = prep_scalars(scalars);
            let res_p2 = if scalars.len() == 1 {
                let mut out = blst_p2::default();
                unsafe {
                    blst_p2_mult(
                        &mut out as *mut blst_p2,
                        &self.first.expect("len != 0"),
                        scalars_le.as_ptr(),
                        255,
                    )
                }
                out
            } else {
                let a: &[blst_p2_affine] = &self.all[..scalars.len()];
                a.mult(&scalars_le, 255)
            };
            Ok(ark_bls12_381::G2Projective {
                x: ark_ff::QuadExtField {
                    c0: ark_ff::Fp(BigInt(res_p2.x.fp[0].l), PhantomData),
                    c1: ark_ff::Fp(BigInt(res_p2.x.fp[1].l), PhantomData),
                },
                y: ark_ff::QuadExtField {
                    c0: ark_ff::Fp(BigInt(res_p2.y.fp[0].l), PhantomData),
                    c1: ark_ff::Fp(BigInt(res_p2.y.fp[1].l), PhantomData),
                },
                z: ark_ff::QuadExtField {
                    c0: ark_ff::Fp(BigInt(res_p2.z.fp[0].l), PhantomData),
                    c1: ark_ff::Fp(BigInt(res_p2.z.fp[1].l), PhantomData),
                },
            })
        }
    }

    fn convert_g1(p: ark_bls12_381::G1Projective) -> blst_p1 {
        let x = blst_fp { l: p.x.0 .0 };
        let y = blst_fp { l: p.y.0 .0 };
        let z = blst_fp { l: p.z.0 .0 };
        blst_p1 { x, y, z }
    }

    fn convert_g2(p: ark_bls12_381::G2Projective) -> blst_p2 {
        let x = blst_fp2 {
            fp: [blst_fp { l: p.x.c0.0 .0 }, blst_fp { l: p.x.c1.0 .0 }],
        };
        let y = blst_fp2 {
            fp: [blst_fp { l: p.y.c0.0 .0 }, blst_fp { l: p.y.c1.0 .0 }],
        };
        let z = blst_fp2 {
            fp: [blst_fp { l: p.z.c0.0 .0 }, blst_fp { l: p.z.c1.0 .0 }],
        };
        blst_p2 { x, y, z }
    }

    fn convert_g1_affine(p: ark_bls12_381::G1Affine) -> blst_p1_affine {
        blst_p1_affine {
            x: blst_fp { l: p.x.0 .0 },
            y: blst_fp { l: p.y.0 .0 },
        }
    }

    fn convert_g2_affine(q: ark_bls12_381::G2Affine) -> blst_p2_affine {
        blst_p2_affine {
            x: blst_fp2 {
                fp: [blst_fp { l: q.x.c0.0 .0 }, blst_fp { l: q.x.c1.0 .0 }],
            },
            y: blst_fp2 {
                fp: [blst_fp { l: q.y.c0.0 .0 }, blst_fp { l: q.y.c1.0 .0 }],
            },
        }
    }

    fn prep_scalars(scalars: &[ark_bls12_381::Fr]) -> Vec<u8> {
        let mut scalars_le = vec![0u8; 32 * scalars.len()];
        for (i, s) in scalars.iter().enumerate() {
            // This _must_ be little endian bytes for this to work
            s.serialize_compressed(&mut scalars_le[i * 32..(i + 1) * 32])
                .unwrap();
        }
        scalars_le
    }

    fn check_scalars(scalars: &[ark_bls12_381::Fr], point_len: usize) -> Result<(), Error> {
        if point_len < scalars.len() {
            return Err(Error::TooManyScalars {
                n_coeffs: scalars.len(),
                expected_max: point_len,
            });
        }
        Ok(())
    }

    fn trim_zeros(mut scalars: &[ark_bls12_381::Fr]) -> &[ark_bls12_381::Fr] {
        while scalars.last().map(|s| s.is_zero()).unwrap_or(false) {
            scalars = &scalars[..scalars.len() - 1];
        }
        scalars
    }
    #[cfg(test)]
    mod tests {
        use ark_bls12_381::{Fr, G1Affine, G1Projective, G2Affine, G2Projective};
        use ark_ec::{pairing::Pairing, CurveGroup};
        use ark_ff::UniformRand;
        use ark_std::vec::Vec;
        use rand::thread_rng;

        use crate::curve_msm;

        use super::*;

        #[test]
        fn test_msm_errors() {
            let g1s = (0..512)
                .map(|_| ark_bls12_381::G1Projective::rand(&mut thread_rng()).into_affine())
                .collect::<Vec<_>>();

            let scalars = (0..512)
                .map(|_| ark_bls12_381::Fr::rand(&mut thread_rng()))
                .collect::<Vec<_>>();

            fn run(err: bool, g1s: &[ark_bls12_381::G1Affine], scalars: &[ark_bls12_381::Fr]) {
                let pg1 = PreparedG1s::from_affines(g1s.to_vec());
                let res = pg1.msm(&scalars);
                assert_eq!(res.is_err(), err);
            }

            run(true, &[], scalars[0..2].as_ref());
            run(true, g1s[0..1].as_ref(), scalars[0..2].as_ref());
            run(true, g1s[0..10].as_ref(), scalars[0..20].as_ref());
            run(true, g1s[0..19].as_ref(), scalars[0..20].as_ref());
            run(false, g1s[0..20].as_ref(), scalars[0..20].as_ref());
        }

        #[test]
        fn test_msm_works() {
            let g1s = (0..512)
                .map(|_| ark_bls12_381::G1Projective::rand(&mut thread_rng()).into_affine())
                .collect::<Vec<_>>();
            let g2s = (0..512)
                .map(|_| ark_bls12_381::G2Projective::rand(&mut thread_rng()).into_affine())
                .collect::<Vec<_>>();
            let scalars = (0..512)
                .map(|_| ark_bls12_381::Fr::rand(&mut thread_rng()))
                .collect::<Vec<_>>();

            let pg1 = PreparedG1s::from_affines(g1s.clone());
            let pg2 = PreparedG2s::from_affines(g2s.clone());

            let res1 = pg1.msm(&scalars).unwrap();
            let res2 = pg2.msm(&scalars).unwrap();

            let alt_res1 = curve_msm::<ark_bls12_381::G1Projective>(&g1s, &scalars).unwrap();
            let alt_res2 = curve_msm::<ark_bls12_381::G2Projective>(&g2s, &scalars).unwrap();

            assert_eq!(res1, alt_res1);
            assert_eq!(res2, alt_res2);
        }

        #[test]
        fn test_single_works() {
            let g1s = vec![ark_bls12_381::G1Affine::rand(&mut thread_rng())];
            let g2s = vec![ark_bls12_381::G2Affine::rand(&mut thread_rng())];
            let scalars = vec![ark_bls12_381::Fr::rand(&mut thread_rng())];

            let pg1 = PreparedG1s::from_affines(g1s.to_vec());
            let pg2 = PreparedG2s::from_affines(g2s.to_vec());

            let res1 = pg1.msm(&scalars).unwrap();
            let res2 = pg2.msm(&scalars).unwrap();

            let alt_res1 = curve_msm::<ark_bls12_381::G1Projective>(&g1s, &scalars).unwrap();
            let alt_res2 = curve_msm::<ark_bls12_381::G2Projective>(&g2s, &scalars).unwrap();

            assert_eq!(res1, alt_res1);
            assert_eq!(res2, alt_res2);
        }

        #[test]
        fn test_no_scalars_returns_zero() {
            let g1s = vec![G1Affine::rand(&mut thread_rng())];
            let g2s = vec![G2Affine::rand(&mut thread_rng())];

            let pg1 = PreparedG1s::from_affines(g1s);
            let pg2 = PreparedG2s::from_affines(g2s);

            assert_eq!(Ok(G1Projective::zero()), pg1.msm(&[]));
            assert_eq!(Ok(G2Projective::zero()), pg2.msm(&[]));
        }

        #[test]
        fn test_all_zeros_returns_zero() {
            let g1s = (0..10)
                .map(|_| G1Projective::rand(&mut thread_rng()).into_affine())
                .collect::<Vec<_>>();
            let g2s = (0..10)
                .map(|_| G2Projective::rand(&mut thread_rng()).into_affine())
                .collect::<Vec<_>>();

            let pg1 = PreparedG1s::from_affines(g1s);
            let pg2 = PreparedG2s::from_affines(g2s);

            let scalars = vec![Zero::zero(); 10];
            let res1 = pg1.msm(&scalars);
            let res2 = pg2.msm(&scalars);

            assert_eq!(Ok(G1Projective::zero()), res1);
            assert_eq!(Ok(G2Projective::zero()), res2);
        }

        #[test]
        fn test_pairings() {
            use ark_ff::One;
            let z = Fr::rand(&mut thread_rng());
            let p1 = G1Affine::rand(&mut thread_rng());
            let q1 = G2Affine::rand(&mut thread_rng());

            let p2 = (p1 * z).into_affine();
            let q2 = (q1 * (Fr::one() / z)).into_affine();

            let ref1 = ark_bls12_381::Bls12_381::pairing(p1, q1).0;
            let ref2 = ark_bls12_381::Bls12_381::pairing(p2, q2).0;
            let pair1 = pairing(p1, q1);
            let pair2 = pairing(p2, q2);
            assert_eq!(ref1, pair1);
            assert_eq!(ref2, pair2);
            assert_eq!(ref1, ref2);
            assert!(check_pairings_equal(p1, q1, p2, q2));
        }
    }
}
