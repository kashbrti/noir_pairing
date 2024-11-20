use crate::fp2::Fp2;

/// BN254 G2 Curve in Affine Space
pub struct G2Affine {
    /// Element of `Fp2`
    pub x: Fp2,
    /// Element of `Fp2`
    pub y: Fp2,
    /// Defined as `true` if point is at infinity
    pub infinity: bool,
}

impl G2Affine {
    /// Returns the generator point in `G2`.
    pub fn generator() -> Self {
        // todo:
        Self {
            x: Fp2::zero(),
            y: Fp2::zero(),
            infinity: false,
        }
    }

    /// Performs Y-axis negation for an element of `G1`.
    pub fn neg(self) -> Self {
        let mut y = self.y.neg();

        if self.infinity {
            y = Fp2::one();
        }

        G2Affine { x: self.x, y, infinity: self.infinity }
    }
}

/// BN254 G2 Curve in Projective Space
pub struct G2Projective {
    /// Element of `Fp2`
    pub x: Fp2,
    /// Element of `Fp2`
    pub y: Fp2,
    /// Element of `Fp2`
    pub z: Fp2,
}

impl G2Projective {
    /// Constructs a projective point from an affine point.
    ///
    /// > TODO: isn't point at infinity for projective space (0, 1, 0)?
    pub fn from_affine(p: G2Affine) -> Self {
        let mut z = Fp2::one();

        if (p.infinity) {
            z = Fp2::zero();
        }

        G2Projective { x: p.x, y: p.y, z }
    }

    /// Performs the doubling step of the Miller Loop.
    ///
    /// Reference Document: https://eprint.iacr.org/2010/354.pdf (Algorithm 26).
    ///
    /// > TODO: rename variables to not be bad (WIP)
    pub fn doubling_step(&mut self) -> (Fp2, Fp2, Fp2) {
        let x_squared = self.x.square();
        let y_squared = self.y.square();
        let z_squared = self.z.square();

        let y_exp_four = y_squared.square();

        // idk what to name these
        let tmp3 = (y_squared + self.x).square() - x_squared - y_exp_four;
        let tmp3 = tmp3 + tmp3;
        let tmp4 = x_squared + x_squared + x_squared;
        let tmp5 = tmp4.square();

        let out_z = self.x + tmp4;

        self.x = tmp5 - tmp3 - tmp3;
        self.z = (self.z + self.y).square() - y_squared - z_squared;
        self.y = (tmp3 - self.x) * tmp4;

        let two_y_exp_four = y_exp_four + y_exp_four;
        let eight_y_exp_four = two_y_exp_four + two_y_exp_four + two_y_exp_four + two_y_exp_four;
        self.y -= eight_y_exp_four;

        let out_y = tmp4 * z_squared;
        let out_y = out_y + out_y;
        let out_y = out_y.neg();

        let out_z = out_z.square() - x_squared - tmp5;
        let four_y_squared = y_squared + y_squared + y_squared + y_squared;
        let out_z = out_z - four_y_squared;
        let out_x = self.z * z_squared;
        let out_x = out_x + out_x;

        (out_x, out_y, out_z)
    }
}
