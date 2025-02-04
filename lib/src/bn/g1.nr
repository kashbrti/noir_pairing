use crate::bn::fp2::Fp;
use bignum::BigNum;

/// BN254 G1 Curve in Affine Space
pub struct G1Affine {
    /// Element of `Fp`
    pub x: Fp,
    /// Element of `Fp`
    pub y: Fp,
    /// Whether the point is the point at infinity (identity)
    pub infinity: bool,
}

impl G1Affine {
    /// Returns the generator point in `G1`.
    pub fn generator() -> Self {
        Self {
            x: BigNum { limbs: [0x01, 0x00, 0x00] },
            y: BigNum { limbs: [0x02, 0x00, 0x00] },
            infinity: false,
        }
    }

    /// Performs Y-axis negation for an element of `G1`.
    pub fn neg(self: Self) -> Self {
        let mut y = self.y.neg();

        if self.infinity {
            y = BigNum::one();
        }

        Self { x: self.x, y, infinity: self.infinity }
    }
}
