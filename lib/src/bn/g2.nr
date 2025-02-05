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
    // reference: https://github.com/arkworks-rs/curves/blob/master/bn254/src/curves/g2.rs
    // xQ = 0x19B0BEA4AFE4C330DA93CC3533DA38A9F430B471C6F8A536E81962ED967909B5
    // + 0xA1CF 585585A61C6E9880B1F2A5C539F7D906FFF238FA6341E1DE1A2E45C3F72u,
    // yQ = 0x17ABD366EBBD65333E49C711A80A0CF6D24ADF1B9B3990EEDCC91731384D2627
    // + 0xEE97D6DE9902A27D00E952232A78700863BC9AA9BE960C32F5BF9FD0A32D345u.
    pub fn generator() -> Self {
        G2Affine {
            x: Fp2::new(
                [348590364376280642591967515799385837, 1157538834227791282941611274091657575, 6144],
                [884175139096380689962467002185618114, 766260723838707333385654073554707953, 6542],
            ),
            y: Fp2::new(
                [1089883602805103349575321108707179946, 491439896715252369454387559521095651, 4808],
                [390420177976373117946371462389733211, 715570412586907199292313128083363260, 2310],
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
