use crate::bn::fp2::Fp2;

/// BN254 G2 Curve in Affine Space
pub struct G2Affine {
    /// X-coordinate in `Fp2`
    pub x: Fp2,
    /// Y-coordinate in `Fp2`
    pub y: Fp2,
    /// Whether the point is the point at infinity (identity)
    pub infinity: bool,
}

/// BN254 G2 Curve in Projective Space
pub struct G2Projective {
    /// X-coordinate in `Fp2`
    pub x: Fp2,
    /// Y-coordinate in `Fp2`
    pub y: Fp2,
    /// Z-coordinate in `Fp2`
    pub z: Fp2,
}

impl G2Affine {
    /// Returns the generator point in `G2`.
    pub fn generator() -> Self {
        G2Affine {
            x: Fp2::new(
                [0x1800deef121f1e76426a00665e5c44, 0x79674322d4f75edadd46debd5cd992, 0xf6ed],
                [0x198e9393920d483a7260bfb731fb5d, 0x25f1aa493335a9e71297e485b7aef3, 0x12c2],
            ),
            y: Fp2::new(
                [0x12c85ea5db8c6deb4aab71808dcb40, 0x8fe3d1e7690c43d37b4ce6cc0166fa, 0x7daa],
                [0x90689d0585ff075ec9e99ad690c339, 0x5bc4b313370b38ef355acdadcd1229, 0x75b],
            ),
            infinity: false,
        }
    }

    /// Performs Y-axis negation for an element in `G2`.
    pub fn neg(self: Self) -> Self {
        let mut y = self.y.neg();
        if self.infinity {
            y = Fp2::one();
        }
        G2Affine { x: self.x, y, infinity: self.infinity }
    }
}

impl G2Projective {
    /// Transforms a point from Affine Space to Projective Space
    pub fn from_affine(p: G2Affine) -> Self {
        let mut z = Fp2::one();

        if (p.infinity) {
            z = Fp2::zero();
        }

        G2Projective { x: p.x, y: p.y, z }
    }

    pub fn neg(self: Self) -> Self {
        G2Projective { x: self.x, y: self.y.neg(), z: self.z }
    }
}
