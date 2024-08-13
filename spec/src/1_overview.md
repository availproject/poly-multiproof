# Overview

The Polynomial Multiproof (PMP) scheme is a polynomial commitment scheme that allows for efficiently creating/verifying opening proofs for multiple polynomials at multiple points. 

Here's why this project is cool: **Opening gets faster the more points we open a polynomial at.**

This is a huge deal for data availability systems, which are bottlenecked by opening time.

For some sample numbers, here are times opening/verifying degree 32768 - 1 polynomials using BLS12-381[^1] on a single core of an M1 Macbook Pro:

![Benchmarks](./m1cycl_32768.png)

This specification outlines implementation requirements for the polynomial multiproof (PMP) scheme.
This builds on previous methods such as KZG10[^2], and is derived entirely from BDFG21[^3]. 
It should be seen as choosing a special case of BDFG21 which allows for significant optimizations that make the protocol more viable for use in applications like Data Availability.

The scheme consists of four methods:

1. `Setup` which sets up a structured reference string for the curve and does some 
2. `Commit` which commits to a polynomial
3. `M1Open/M2Open` which computes a single opening proof that a set of polynomials are equal to specific values at a set of points
4. `M1Verify/M2Verify` which verify an opening proof against commitments and evaluations

### API

The API is designed around two traits, the first looks roughly like
```rust
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
```

This takes a Merlin transcript for the Fiat-Shamir transform, some points $x_j$, some polynomials $f_i$ and the evaluations of each polynomial at each point $f_i(x_j)$.

This API is the most flexible, but not the most efficient.

The second trait is more efficient, but requires knowledge of which points are being opened to beforehand.
```rust
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
```

When creating this object, there will be some way to specify the points, or they will be predetermined by the implementation. The API is largely the same.

### Applications to Data Availability
Data availability systems can benefit greatly from PMP, since polynomial commitment-based DA systems are tremendously constrained by opening time.
Even with very fast KZG implementations, opening is too slow, and constrains the amount of data which can be processed through the system

> Opening to a degree 255 BLS12-381 scalar field polynomial takes at best ~3-6ms, which means opening to every cell of 256 of those polynomials will take >200s.

Because PMP allows many cells to be opened to at once, the data availability grid can be chunked into a smaller grid, allowing for a faster, more secure system, with far higher througput [^4].
This allows polynomial commitment-based DA systems scale up in size without taking more compute.

### Curve selection

The protocol depends on the selection of a pairing based curve. 
Attributes of an ideal curve for this protocol are, in order of importance,
1. Security
2. Fast G1 scalar multiplication
3. Large scalar field size
4. Fast G2 scalar multiplication
5. Small compressed G1 scalar size

BLS12-381 satisfies the first two critera well, and while other curves likely could satisfy more requirements, they are currently not as well audited. Hence BLS12-381 is a reasonable choice to start with.

The scheme given in BDFG21 is interactive, so Merlin transcripts are used to do the Fiat-Shamir transform.

Optimizations made from BDFG21 are outlined in [Optimizations](./a_optimizations.md)

[^1]: All times are with a Xeon E5-2676 v3 @ 2.40GHz vCPU on AWS

[^2]: [KZG10 Paper](https://www.iacr.org/archive/asiacrypt2010/6477178/6477178.pdf)

[^3]: [BDFG21 Paper](https://eprint.iacr.org/2020/081.pdf)

[^4]: More FFTs are required as data size increases, but these are dwarfed by opening times. Additionally, exactly how the chunking is done has nuanced security implications in different scenarios. 
