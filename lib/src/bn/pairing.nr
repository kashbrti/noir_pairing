use crate::bn::fp2::Fp2;
use crate::bn::g1::G1Affine;

/// Encapsulates sparse non-zero coefficients in Fp12.
///
/// There are multiple operations over elements of Fp12 where the majority of coefficients are zero,
/// which results in unnecessary operations.
///
/// We define a data structure representing a more succinct form of the following.
///
/// ```
/// Fp12 {
///     c0: Fp6 {
///        c0: Fp2 { .. },
///        c1: Fp2::zero(),
///        c2: Fp2::zero(),
///     },
///     c1: Fp6 {
///        c0: Fp2 { .. },
///        c1: Fp2 { .. },
///        c2: Fp2::zero(),
///     },
/// }
/// ```
///
/// In particular, we define this for the point doubling/addition and line evaluation algorithms.
///
/// Reference Document: https://eprint.iacr.org/2010/354.pdf (Algorithm 26, 27).
pub(crate) struct Ell {
    pub(crate) c0_c0: Fp2,
    pub(crate) c1_c0: Fp2,
    pub(crate) c1_c1: Fp2,
}

impl Ell {
    /// Constructs an Ell instance with coefficients in Fp2.
    pub(crate) fn new(c0_c0: Fp2, c1_c0: Fp2, c1_c1: Fp2) -> Self {
        Self { c0_c0, c1_c0, c1_c1 }
    }

    /// Scales the coefficients of an Ell instance by an affine point.
    pub(crate) fn scale_by_affine(self, p: G1Affine) -> Self {
        Self { c0_c0: self.c0_c0, c1_c0: self.c1_c0.scale(p.x), c1_c1: self.c1_c1.scale(p.y) }
    }
}
