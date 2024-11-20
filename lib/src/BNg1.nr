use crate::BNfp2::Fp;

/// BN254 G1 Curve in Affine Space
pub struct G1Affine {
    /// Element of `Fp`
    pub x: Fp,
    /// Element of `Fp`
    pub y: Fp,
    /// Defined as `true` if point is at infinity
    pub infinity: bool,
}

impl G1Affine {
    /// Returns the generator point in `G1`.
    pub fn generator() -> Self {
        // todo:
        Self {
            x: Fp::new(),
            y: Fp::new(),
            infinity: false,
        }
    }

    /// Performs Y-axis negation for an element of `G1`.
    pub fn neg(self: Self) -> Self {
        let mut y = self.y.neg();

        if self.infinity {
            y = Fp::one();
        }

        Self { x: self.x, y, infinity: self.infinity }
    }
}
