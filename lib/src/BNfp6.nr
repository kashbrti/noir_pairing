use crate::BNfp2::Fp;
use crate::BNfp2::Fp2;
use bignum::BigNum;

// c0 + c1*v + c2*v^2
// F_{p^6} =  F_{p^2}[v] / (v^3 - (9+u))
pub struct Fp6 {
    pub(crate) c0: Fp2,
    pub(crate) c1: Fp2,
    pub(crate) c2: Fp2,
}

// TODO (if possible) make underlying (Fp2) function unconstrained and then constrain for better performance

impl Fp6 {
    pub fn new(c0: Fp2, c1: Fp2, c2: Fp2) -> Self {
        Fp6 { c0, c1, c2 }
    }

    pub fn zero() -> Self {
        Fp6 { c0: Fp2::zero(), c1: Fp2::zero(), c2: Fp2::zero() }
    }

    pub fn one() -> Self {
        Fp6 { c0: Fp2::one(), c1: Fp2::zero(), c2: Fp2::zero() }
    }

    pub fn neg(self: Self) -> Self {
        Self { c0: self.c0.neg(), c1: self.c1.neg(), c2: self.c2.neg() }
    }

    pub fn mul_by_1(self: Self, c1: Fp2) -> Fp6 {
        Fp6 { c0: (self.c2 * c1).mul_by_non_residue(), c1: self.c0 * c1, c2: self.c1 * c1 }
    }

    pub fn mul_by_01(self: Self, c0: Fp2, c1: Fp2) -> Fp6 {
        let a_a = self.c0 * c0;
        let b_b = self.c1 * c1;

        let t1 = ((self.c2 * c1).mul_by_non_residue()) + a_a;

        let t2 = (c0 + c1) * (self.c0 + self.c1) - a_a - b_b;

        let t3 = self.c2 * c0 + b_b;

        Fp6 { c0: t1, c1: t2, c2: t3 }
    }

    // v is a quadratic residue in Fp6
    // returns self * v
    pub fn mul_by_nonresidue(self: Self) -> Self {
        // Multiply c0 + c1 * v + c2 * v^2 by v
        // (c0 + c1 * v + c2 * v^2) * v
        // = c0 * v + c1 * v^2 + c2 * v^3
        // And v^3 equals (u + 1):
        // = c0 * v + c1 * v^2 + c2 * (u + 1)
        // = c0 * v + c1 * v^2 + c2 * u + c2
        // rewrite according to the coefficients
        // = (c2 * u + c2) + (c0 * v) + (c1 * v^2)
        // = (c2 * (u + 1)) + (c0 * v) + (c1 * v^2)
        Fp6 {
            c0: self.c2.mul_by_non_residue(), // c2 * (u + 1)
            c1: self.c0, // c0
            c2: self.c1, // c1
        }
    }

    // https://eprint.iacr.org/2022/367
    // for inputs `a = (a0, a1, a2)` and `b = (b0, b1, b2)`,
    // return multiplication a*b
    fn mul_interleaved(a: Self, b: Self) -> Self {
        unsafe {
            // Unconstrained calculations first
            //   result_c0_0 = a00 * b00 - a01 * b01 + a10 * b20 - a11 * b21 + a20 * b10 - a21 * b11
            //                                        - a10 * b21 - a11 * b20 - a20 * b11 - a21 * b10.
            //               = a00 * b00 - a01 * b01 + a10 * (b20 - b21) - a11 * (b20 + b21)
            //                                        + a20 * (b10 - b11) - a21 * (b10 + b11).
            let a00b00 = a.c0.c0.__mul(b.c0.c0);
            let a01b01 = a.c0.c1.__mul(b.c0.c1);
            let b20_sub_b21 = b.c2.c0.__sub(b.c2.c1);
            let a10b20_sub_b21 = a.c1.c0.__mul(b20_sub_b21);
            let b20_add_b21 = b.c2.c0.__add(b.c2.c1);
            let a11b20_add_b21 = a.c1.c1.__mul(b20_add_b21);
            let b10_sub_b11 = b.c1.c0.__sub(b.c1.c1);
            let a20b10_sub_b11 = a.c2.c0.__mul(b10_sub_b11);
            let b10_add_b11 = b.c1.c0.__add(b.c1.c1);
            let a21b10_add_b11 = a.c2.c1.__mul(b10_add_b11);

            let mut result_c0_0 = a00b00.__sub(a01b01);
            result_c0_0 = result_c0_0.__add(a10b20_sub_b21);
            result_c0_0 = result_c0_0.__sub(a11b20_add_b21);
            result_c0_0 = result_c0_0.__add(a20b10_sub_b11);
            result_c0_0 = result_c0_0.__sub(a21b10_add_b11);

            //   result_c0_1 = a00 * b01 + a01 * b00 + a10 * b21 + a11 * b20 + a20 * b11 + a21 * b10
            //                                        + a10 * b20 - a11 * b21 + a20 * b10 - a21 * b11.
            //               = a00 * b01 + a01 * b00 + a10 * (b20 + b21) + a11 * (b20 - b21)
            //                                        + a20 * (b10 + b11) + a21 * (b10 - b11).
            let a00b01 = a.c0.c0.__mul(b.c0.c1);
            let a01b00 = a.c0.c1.__mul(b.c0.c0);
            let b20_add_b21 = b.c2.c0.__add(b.c2.c1);
            let a10b20_add_b21 = a.c1.c0.__mul(b20_add_b21);
            let b20_sub_b21 = b.c2.c0.__sub(b.c2.c1);
            let a11b20_sub_b21 = a.c1.c1.__mul(b20_sub_b21);
            let b10_add_b11 = b.c1.c0.__add(b.c1.c1);
            let a20b10_add_b11 = a.c2.c0.__mul(b10_add_b11);
            let b10_sub_b11 = b.c1.c0.__sub(b.c1.c1);
            let a21b10_sub_b11 = a.c2.c1.__mul(b10_sub_b11);

            let mut result_c0_1 = a00b01.__add(a01b00);
            result_c0_1 = result_c0_1.__add(a10b20_add_b21);
            result_c0_1 = result_c0_1.__add(a11b20_sub_b21);
            result_c0_1 = result_c0_1.__add(a20b10_add_b11);
            result_c0_1 = result_c0_1.__add(a21b10_sub_b11);

            //   result_c1_0 = a00 * b10 - a01 * b11 + a10 * b00 - a11 * b01 + a20 * b20 - a21 * b21
            //                                                          - a20 * b21 - a21 * b20.
            //               = a00 * b10 - a01 * b11 + a10 * b00 - a11 * b01 + a20 * (b20 - b21)
            //                                                          - a21 * (b20 + b21).
            let a00b10 = a.c0.c0.__mul(b.c1.c0);
            let a01b11 = a.c0.c1.__mul(b.c1.c1);
            let a10b00 = a.c1.c0.__mul(b.c0.c0);
            let a11b01 = a.c1.c1.__mul(b.c0.c1);
            let b20_sub_b21 = b.c2.c0.__sub(b.c2.c1);
            let a20b20_sub_b21 = a.c2.c0.__mul(b20_sub_b21);
            let b20_add_b21 = b.c2.c0.__add(b.c2.c1);
            let a21b20_add_b21 = a.c2.c1.__mul(b20_add_b21);

            let mut result_c1_0 = a00b10.__sub(a01b11);
            result_c1_0 = result_c1_0.__add(a10b00);
            result_c1_0 = result_c1_0.__sub(a11b01);
            result_c1_0 = result_c1_0.__add(a20b20_sub_b21);
            result_c1_0 = result_c1_0.__sub(a21b20_add_b21);

            //   result_c1_1 = a00 * b11 + a01 * b10 + a10 * b01 + a11 * b00 + a20 * b21 + a21 * b20
            //                                                          + a20 * b20 - a21 * b21.
            //               = a00 * b11 + a01 * b10 + a10 * b01 + a11 * b00 + a20 * (b20 + b21)
            //                                                          + a21 * (b20 - b21).
            let a00b11 = a.c0.c0.__mul(b.c1.c1);
            let a01b10 = a.c0.c1.__mul(b.c1.c0);
            let a10b01 = a.c1.c0.__mul(b.c0.c1);
            let a11b00 = a.c1.c1.__mul(b.c0.c0);
            let b20_add_b21 = b.c2.c0.__add(b.c2.c1);
            let a20b20_add_b21 = a.c2.c0.__mul(b20_add_b21);
            let b20_sub_b21 = b.c2.c0.__sub(b.c2.c1);
            let a21b20_sub_b21 = a.c2.c1.__mul(b20_sub_b21);

            let mut result_c1_1 = a00b11.__add(a01b10);
            result_c1_1 = result_c1_1.__add(a10b01);
            result_c1_1 = result_c1_1.__add(a11b00);
            result_c1_1 = result_c1_1.__add(a20b20_add_b21);
            result_c1_1 = result_c1_1.__add(a21b20_sub_b21);

            //   result_c2_0 = a00 * b20 - a01 * b21 + a10 * b10 - a11 * b11 + a20 * b00 - a21 * b01.
            let a00b20 = a.c0.c0.__mul(b.c2.c0);
            let a01b21 = a.c0.c1.__mul(b.c2.c1);
            let a10b10 = a.c1.c0.__mul(b.c1.c0);
            let a11b11 = a.c1.c1.__mul(b.c1.c1);
            let a20b00 = a.c2.c0.__mul(b.c0.c0);
            let a21b01 = a.c2.c1.__mul(b.c0.c1);

            let mut result_c2_0 = a00b20.__sub(a01b21);
            result_c2_0 = result_c2_0.__add(a10b10);
            result_c2_0 = result_c2_0.__sub(a11b11);
            result_c2_0 = result_c2_0.__add(a20b00);
            result_c2_0 = result_c2_0.__sub(a21b01);

            //   result_c2_1 = a00 * b21 + a01 * b20 + a10 * b11 + a11 * b10 + a20 * b01 + a21 * b00.
            let a00b21 = a.c0.c0.__mul(b.c2.c1);
            let a01b20 = a.c0.c1.__mul(b.c2.c0);
            let a10b11 = a.c1.c0.__mul(b.c1.c1);
            let a11b10 = a.c1.c1.__mul(b.c1.c0);
            let a20b01 = a.c2.c0.__mul(b.c0.c1);
            let a21b00 = a.c2.c1.__mul(b.c0.c0);

            let mut result_c2_1 = a00b21.__add(a01b20);
            result_c2_1 = result_c2_1.__add(a10b11);
            result_c2_1 = result_c2_1.__add(a11b10);
            result_c2_1 = result_c2_1.__add(a20b01);
            result_c2_1 = result_c2_1.__add(a21b00);

            /////////// CONSTRAINTS ///////////
            // Constraints for result_c0_0
            // = a00 * b00 - a01 * b01 + a10 * (b20 - b21) - a11 * (b20 + b21)
            //                                        + a20 * (b10 - b11) - a21 * (b10 + b11).
            let lhs_terms_c0_0 = [[a.c0.c0], [a.c0.c1], [a.c1.c0], [a.c1.c1], [a.c2.c0], [a.c2.c1]];
            let lhs_flags_c0_0 = [[false], [true], [false], [true], [false], [true]];
            let rhs_terms_c0_0 =
                [[b.c0.c0], [b.c0.c1], [b20_sub_b21], [b20_add_b21], [b10_sub_b11], [b10_add_b11]];
            let rhs_flags_c0_0 = [[false], [false], [false], [false], [false], [false]];
            let res_c0_0 = [result_c0_0];
            let res_flag_c0_0 = [true];
            BigNum::evaluate_quadratic_expression(
                lhs_terms_c0_0,
                lhs_flags_c0_0,
                rhs_terms_c0_0,
                rhs_flags_c0_0,
                res_c0_0,
                res_flag_c0_0,
            );

            // Constraints for result_c0_1
            //   = a00 * b01 + a01 * b00 + a10 * (b20 + b21) + a11 * (b20 - b21)
            //                                        + a20 * (b10 + b11) + a21 * (b10 - b11).
            let lhs_terms_c0_1 = [[a.c0.c0], [a.c0.c1], [a.c1.c0], [a.c1.c1], [a.c2.c0], [a.c2.c1]];
            let lhs_flags_c0_1 = [[false], [false], [false], [false], [false], [false]];
            let rhs_terms_c0_1 =
                [[b.c0.c1], [b.c0.c0], [b20_add_b21], [b20_sub_b21], [b10_add_b11], [b10_sub_b11]];
            let rhs_flags_c0_1 = [[false], [false], [false], [false], [false], [false]];
            let res_c0_1 = [result_c0_1];
            let res_flag_c0_1 = [true];
            BigNum::evaluate_quadratic_expression(
                lhs_terms_c0_1,
                lhs_flags_c0_1,
                rhs_terms_c0_1,
                rhs_flags_c0_1,
                res_c0_1,
                res_flag_c0_1,
            );

            // Constraints for result_c1_0
            // = a00 * b10 - a01 * b11 + a10 * b00 - a11 * b01 + a20 * (b20 - b21)
            //                                                          - a21 * (b20 + b21).
            let lhs_terms_c1_0 = [[a.c0.c0], [a.c0.c1], [a.c1.c0], [a.c1.c1], [a.c2.c0], [a.c2.c1]];
            let lhs_flags_c1_0 = [[false], [true], [false], [true], [false], [true]];
            let rhs_terms_c1_0 =
                [[b.c1.c0], [b.c1.c1], [b.c0.c0], [b.c0.c1], [b20_sub_b21], [b20_add_b21]];
            let rhs_flags_c1_0 = [[false], [false], [false], [false], [false], [false]];
            let res_c1_0 = [result_c1_0];
            let res_flag_c1_0 = [true];
            BigNum::evaluate_quadratic_expression(
                lhs_terms_c1_0,
                lhs_flags_c1_0,
                rhs_terms_c1_0,
                rhs_flags_c1_0,
                res_c1_0,
                res_flag_c1_0,
            );

            // Constraints for result_c1_1
            // = a00 * b11 + a01 * b10 + a10 * b01 + a11 * b00 + a20 * (b20 + b21)
            //                                                          + a21 * (b20 - b21).
            let lhs_terms_c1_1 = [[a.c0.c0], [a.c0.c1], [a.c1.c0], [a.c1.c1], [a.c2.c0], [a.c2.c1]];
            let lhs_flags_c1_1 = [[false], [false], [false], [false], [false], [false]];
            let rhs_terms_c1_1 =
                [[b.c1.c1], [b.c1.c0], [b.c0.c1], [b.c0.c0], [b20_add_b21], [b20_sub_b21]];
            let rhs_flags_c1_1 = [[false], [false], [false], [false], [false], [false]];
            let res_c1_1 = [result_c1_1];
            let res_flag_c1_1 = [true];
            BigNum::evaluate_quadratic_expression(
                lhs_terms_c1_1,
                lhs_flags_c1_1,
                rhs_terms_c1_1,
                rhs_flags_c1_1,
                res_c1_1,
                res_flag_c1_1,
            );

            // Constraints for result_c2_0
            //   result_c2_0 = a00 * b20 - a01 * b21 + a10 * b10 - a11 * b11 + a20 * b00 - a21 * b01.
            let lhs_terms_c2_0 = [[a.c0.c0], [a.c0.c1], [a.c1.c0], [a.c1.c1], [a.c2.c0], [a.c2.c1]];
            let lhs_flags_c2_0 = [[false], [true], [false], [true], [false], [true]];
            let rhs_terms_c2_0 = [[b.c2.c0], [b.c2.c1], [b.c1.c0], [b.c1.c1], [b.c0.c0], [b.c0.c1]];
            let rhs_flags_c2_0 = [[false], [false], [false], [false], [false], [false]];
            let res_c2_0 = [result_c2_0];
            let res_flag_c2_0 = [true];
            BigNum::evaluate_quadratic_expression(
                lhs_terms_c2_0,
                lhs_flags_c2_0,
                rhs_terms_c2_0,
                rhs_flags_c2_0,
                res_c2_0,
                res_flag_c2_0,
            );

            // Constraints for result_c2_1
            //   result_c2_1 = a00 * b21 + a01 * b20 + a10 * b11 + a11 * b10 + a20 * b01 + a21 * b00.
            let lhs_terms_c2_1 = [[a.c0.c0], [a.c0.c1], [a.c1.c0], [a.c1.c1], [a.c2.c0], [a.c2.c1]];
            let lhs_flags_c2_1 = [[false], [false], [false], [false], [false], [false]];
            let rhs_terms_c2_1 = [[b.c2.c1], [b.c2.c0], [b.c1.c1], [b.c1.c0], [b.c0.c1], [b.c0.c0]];
            let rhs_flags_c2_1 = [[false], [false], [false], [false], [false], [false]];
            let res_c2_1 = [result_c2_1];
            let res_flag_c2_1 = [true];
            BigNum::evaluate_quadratic_expression(
                lhs_terms_c2_1,
                lhs_flags_c2_1,
                rhs_terms_c2_1,
                rhs_flags_c2_1,
                res_c2_1,
                res_flag_c2_1,
            );

            Fp6 {
                c0: Fp2 { c0: result_c0_0, c1: result_c0_1 },
                c1: Fp2 { c0: result_c1_0, c1: result_c1_1 },
                c2: Fp2 { c0: result_c2_0, c1: result_c2_1 },
            }
        }
    }

    pub fn square(self: Self) -> Self {
        // s0 = self.c0^2
        let a_0 = self.c0;
        let a_1 = self.c1;
        let a_2 = self.c2;

        let c_3 = (a_1 * a_2 + a_2 * a_1).mul_by_non_residue();
        let c_4 = (a_2 * a_2).mul_by_non_residue();

        let c_0 = a_0 * a_0 + c_3;
        let c_1 = a_0 * a_1 + a_1 * a_0 + c_4;
        let c_2 = a_0 * a_2 + a_1 * a_1 + a_2 * a_0 + c_3;

        Self { c0: c_0 , c1: c_1, c2: c_2 }
    }

    // Returns inverse and errors if there is no inverse
    // pub fn invert(self: Self) -> Self {
    //     unsafe {
    //         let mut c0 = self.c1.mul(self.c2).mul_by_nonresidue();
    //         c0 = self.c0.square() - c0;

    //         let mut c1 = self.c2.square().mul_by_nonresidue();
    //         c1 = c1 - (self.c0 * self.c1);

    //         let mut c2 = self.c1.square();
    //         c2 = c2 - (self.c0 * self.c2);

    //         let tmp = self.c1.mul(c2).add(self.c2.mul(c1)).mul_by_nonresidue().add(self.c0.mul(c0));
    //         // This will error if there is no inverse
    //         let inv_tmp = tmp.invert();

    //         let res_c0 = c0 * inv_tmp;
    //         let res_c1 = c1 * inv_tmp;
    //         let res_c2 = c2 * inv_tmp;

    //         Self { c0: res_c0, c1: res_c1, c2: res_c2 }
    //     }
    // }
    pub fn invert(self: Self) -> Self {
        let c0 = self.c2;
        let c0 = c0.mul_by_non_residue();
        let c0 = c0.mul(self.c1);
        let c0 = c0.neg();

        let c0s = self.c0;
        let c0s = c0s.square();
        let c0 = c0.add(c0s);

        let c1 = self.c2;
        let c1 = c1.square();
        let c1 = c1.mul_by_non_residue();

        let c01 = self.c0;
        let c01 = c01.mul(self.c1);
        let c1 = c1.sub(c01);

        let c2 = self.c1;
        let c2 = c2.square();
        let c02 = self.c0;
        let c02 = c02.mul(self.c2);
        let c2 = c2.sub(c02);

        let tmp1 = self.c2;
        let tmp1 = tmp1.mul(c1);
        let tmp2 = self.c1;
        let tmp2 = tmp2.mul(c2);
        let tmp1 = tmp1.add(tmp2);
        let tmp1 = tmp1.mul_by_non_residue();
        let tmp2 = self.c0;
        let tmp2 = tmp2.mul(c0);
        let tmp1 = tmp1.add(tmp2);

        let tmp = tmp1.invert();
        let c0 = tmp.mul(c0);
        let c1 = tmp.mul(c1);
        let c2 = tmp.mul(c2);

        Self {
            c0,
            c1,
            c2,
        }         
    }
}

impl std::ops::Add for Fp6 {
    fn add(self: Self, other: Self) -> Self {
        Self { c0: self.c0 + other.c0, c1: self.c1 + other.c1, c2: self.c2 + other.c2 }
    }
}

impl std::ops::Sub for Fp6 {
    fn sub(self: Self, other: Self) -> Self {
        Self { c0: self.c0 - other.c0, c1: self.c1 - other.c1, c2: self.c2 - other.c2 }
    }
}

impl std::ops::Mul for Fp6 {
    // TODO: mudify this to use the optimized mul_interleaved
    // fn mul(self: Self, other: Self) -> Self {
        // let a_0 = self.c0;
        // let a_1 = self.c1;
        // let a_2 = self.c2;

        // let b_0 = other.c0;
        // let b_1 = other.c1;
        // let b_2 = other.c2;

        // let c_3 = (a_1 * b_2 + a_2 * b_1).mul_by_nonresidue();
        // let c_4 = (a_2 * b_2).mul_by_nonresidue();

        // let c_0 = a_0 * b_0 + c_3;
        // let c_1 = a_0 * b_1 + a_1 * b_0 + c_4;
        // let c_2 = a_0 * b_2 + a_1 * b_1 + a_2 * b_0 + c_3;

        // Self { c0: c_0 , c1: c_1 , c2: c_2 }
// }
        fn mul(self: Self, other: Self) -> Self {
        let a_a = self.c0.mul(other.c0);
        let b_b = self.c1.mul(other.c1);
        let c_c = self.c2.mul(other.c2);

        // Construct t1
        let t1 = other.c1;
        let t1 = t1.add(other.c2);

        let tmp = self.c1.add(self.c2);
        let t1 = t1.mul(tmp);
        let t1 = t1.sub(b_b);
        let t1 = t1.sub(c_c);
        let t1 = t1.mul_by_non_residue();
        let t1 = t1.add(a_a);

        // Construct t3
        let t3 = other.c0;
        let t3 = t3.add(other.c2);

        let tmp = self.c0.add(self.c2);
        let t3 = t3.mul(tmp);
        let t3 = t3.sub(a_a);
        let t3 = t3.add(b_b);
        let t3 = t3.sub(c_c);

        // Construct t2
        let t2 = other.c0;
        let t2 = t2.add(other.c1);

        let tmp = self.c0.add(self.c1);
        let t2 = t2.mul(tmp);
        let t2 = t2.sub(a_a);
        let t2 = t2.sub(b_b);
        let c_c_mul_by_nonres = c_c.mul_by_non_residue();
        let t2 = t2.add(c_c_mul_by_nonres);

        Self{
            c0: t1,
            c1: t2,
            c2: t3,
        }
    }
}

impl std::cmp::Eq for Fp6 {
    fn eq(self: Self, other: Self) -> bool {
        (self.c0 == other.c0) & (self.c1 == other.c1) & (self.c2 == other.c2)
    }
}

use std::hash::keccak256;

comptime global SEED_0: [u8; 32] = [
    0x4a, 0x2b, 0x57, 0xf1, 0xd0, 0xd9, 0xe4, 0xbe, 0x39, 0x8a, 0xa9, 0x50, 0x9d, 0x25, 0x36, 0xbe,
    0x60, 0x77, 0xbf, 0x21, 0x82, 0x2e, 0x6d, 0x40, 0x5d, 0x65, 0xac, 0x2d, 0x0f, 0x0b, 0x3e, 0x75,
];

comptime global SEED_1 = keccak256(SEED_0, 32);
comptime global SEED_2 = keccak256(SEED_1, 32);
comptime global SEED_3 = keccak256(SEED_2, 32);
comptime global SEED_4 = keccak256(SEED_3, 32);
comptime global SEED_5 = keccak256(SEED_4, 32);

#[test]
fn test_add_identity() {
    let zero = Fp6::zero();

    let a = Fp6::new(
        Fp2 { c0: Fp::derive_from_seed(SEED_0), c1: Fp::derive_from_seed(SEED_1) },
        Fp2 { c0: Fp::derive_from_seed(SEED_2), c1: Fp::derive_from_seed(SEED_3) },
        Fp2 { c0: Fp::derive_from_seed(SEED_4), c1: Fp::derive_from_seed(SEED_5) },
    );

    assert(a + zero == a);
}

#[test]
fn test_add_negative() {
    let a = Fp6::new(
        Fp2 { c0: Fp::derive_from_seed(SEED_0), c1: Fp::derive_from_seed(SEED_1) },
        Fp2 { c0: Fp::derive_from_seed(SEED_2), c1: Fp::derive_from_seed(SEED_3) },
        Fp2 { c0: Fp::derive_from_seed(SEED_4), c1: Fp::derive_from_seed(SEED_5) },
    );

    let negative_a = a.neg();

    assert(negative_a + a == Fp6::zero());
}

#[test]
fn test_mul_identity() {
    let one = Fp6::one();

    let a = Fp6::new(
        Fp2 { c0: Fp::derive_from_seed(SEED_0), c1: Fp::derive_from_seed(SEED_1) },
        Fp2 { c0: Fp::derive_from_seed(SEED_2), c1: Fp::derive_from_seed(SEED_3) },
        Fp2 { c0: Fp::derive_from_seed(SEED_4), c1: Fp::derive_from_seed(SEED_5) },
    );

    assert(a * one == a);
}

#[test]
fn test_mul_zero() {
    let zero = Fp6::zero();

    let a = Fp6::new(
        Fp2 { c0: Fp::derive_from_seed(SEED_0), c1: Fp::derive_from_seed(SEED_1) },
        Fp2 { c0: Fp::derive_from_seed(SEED_2), c1: Fp::derive_from_seed(SEED_3) },
        Fp2 { c0: Fp::derive_from_seed(SEED_4), c1: Fp::derive_from_seed(SEED_5) },
    );

    assert(a * zero == zero);
}


#[test]
fn test_mul(){
    let a = Fp6 {
        c0: Fp2::new(
            [395488726691442830429832050603406204, 821075477437255834874989749241978678, 11793],
            [201601712830434592179941596462639970, 1250782793790474347909172763646462344, 9430],
        ),
        c1: Fp2::new(
            [1127469633954639076687122832280634501, 384973290833791765662151864406027283, 6577],
            [240807479084837380482481895969134031, 615221729286056607213692266989250466, 3412],
        ),
        c2: Fp2::new(
            [648583634891115237590656777438366888, 196330182240574212351377511437222610, 111],
            [415311780186006812910310129930084039, 207811584817201665093984844957365187, 6876],
        ),
    };



    let b = Fp6 {
        c0: Fp2::new(
            [695402583260866822984243440826061124, 126654564887251644705390175425443902, 8760],
            [152771845825067828801743595962419732, 3232475459128800318034314468993454, 11456],
        ),
        c1: Fp2::new(
            [266656749817180616673516288090981873, 268257111370955002784716808410768013, 4121],
            [842609025697615925249905878814309698, 1304291885525753312721670863999755284, 12039],
        ),
        c2: Fp2::new(
            [340264409667458652850640609510042144, 489147825056556503523154647718538571, 4954],
            [718683732117605645373316913303678770, 228102735559184690191635657008341956, 4458],
        ),
    };

    let c = Fp6 {
        c0: Fp2::new(
            [910388087919099045573295818704136616, 698885803714215798214598663899750780, 7850],
            [1006653378442988779829050963625961836, 836653650915028521313752353638665004, 4618],
        ),
        c1: Fp2::new(
            [197501221668325652832307237194994509, 1163481128857328662740495301069208208, 3687],
            [653110633990805892195234371661870930, 1213730868419835900120920981261017067, 7195],
        ),
        c2: Fp2::new(
            [700429773182954713755658134267344599, 791962368285800061640234898941041415, 9916],
            [1326819526840081681552200089815661525, 1205015809171828786186472600187623679, 609],
        ),
    };
    assert(a * b == c);


}

#[test]
fn test_eq() {
    let a = Fp6::new(
        Fp2 { c0: Fp::derive_from_seed(SEED_0), c1: Fp::derive_from_seed(SEED_1) },
        Fp2 { c0: Fp::derive_from_seed(SEED_2), c1: Fp::derive_from_seed(SEED_3) },
        Fp2 { c0: Fp::derive_from_seed(SEED_4), c1: Fp::derive_from_seed(SEED_5) },
    );

    let b = Fp6::new(
        Fp2 { c0: Fp::derive_from_seed(SEED_0), c1: Fp::derive_from_seed(SEED_1) },
        Fp2 { c0: Fp::derive_from_seed(SEED_2), c1: Fp::derive_from_seed(SEED_3) },
        Fp2 { c0: Fp::derive_from_seed(SEED_4), c1: Fp::derive_from_seed(SEED_5) },
    );

    assert(a == b);
}

#[test]
fn test_neg() {
    let a = Fp6::new(
        Fp2 { c0: Fp::derive_from_seed(SEED_0), c1: Fp::derive_from_seed(SEED_1) },
        Fp2 { c0: Fp::derive_from_seed(SEED_2), c1: Fp::derive_from_seed(SEED_3) },
        Fp2 { c0: Fp::derive_from_seed(SEED_4), c1: Fp::derive_from_seed(SEED_5) },
    );

    let b = Fp6::new(
        Fp2 { c0: Fp::derive_from_seed(SEED_0), c1: Fp::derive_from_seed(SEED_1) },
        Fp2 { c0: Fp::derive_from_seed(SEED_2), c1: Fp::derive_from_seed(SEED_3) },
        Fp2 { c0: Fp::derive_from_seed(SEED_4), c1: Fp::derive_from_seed(SEED_5) },
    );

    let negative_a = a.neg();

    assert(negative_a + a == Fp6::zero());
    assert(b - a == b + negative_a);
}

// #[test]
// fn test_invert() {
//     let a = Fp6 {
//         c0: Fp2::new(
//             [
//                 0x71825b46349b307298025eea77770c,
//                 0x225cbdc399bb960534036f625a3592,
//                 0xb32ff3d468eb1aba35968a2bbdb9d0,
//                 0x66e47,
//             ],
//             [
//                 0x3fc9a610ece44f3c738a358aa60150,
//                 0x602cf426b135e0e0d837a8cf5d9ec7,
//                 0x4aee9cb5468ad40ee5dccbf0264407,
//                 0xdc722,
//             ],
//         ),
//         c1: Fp2::new(
//             [
//                 0xb77accf6a4e5e337a85d2a62baf070,
//                 0x5e5a337d7db6344484087a18586c7c,
//                 0xb3f92a8e29eabd22358787933cd419,
//                 0xa4aa,
//             ],
//             [
//                 0xb33c802bede6916fa9f3d1d08a643b,
//                 0x4537aef2966b539d7364999eb35707,
//                 0xde1e413ba79349fbe6cc88298b9f18,
//                 0x1776be,
//             ],
//         ),
//         c2: Fp2::new(
//             [
//                 0x8a449468e00f21f8f78f53b57f2f81,
//                 0x618036c76149747cdb368c1d1ccf7a,
//                 0x55fe3e909dbe81264b06f976721808,
//                 0x122bf8,
//             ],
//             [
//                 0x157c33763f6eb6ec8697a38c6cb1d2,
//                 0xbae7a1aac9a8926b246110708adc2b,
//                 0x8c3ab0276bcf5791fb656da17e25c6,
//                 0xe008d,
//             ],
//         ),
//     };

//     let a_inv = a.invert();
//     let check = a_inv.mul_interleaved(a);
//     assert(check.eq(Fp6::one()));
// }

// #[test]
// fn test_arithmetic() {
//     let a = Fp6 {
//         c0: Fp2::new(
//             [
//                 0x716afa99ba71dc4ebb9630cea5ec71,
//                 0xd75cc10ca49025b37bcd6a1d508b3f,
//                 0x3194a03544948b7acc2dabfe9d6111,
//                 0x130426,
//             ],
//             [
//                 0xa927e8aecb68b2d55389a02df3fe47,
//                 0x472caff2eec27e41a561c40ddd3188,
//                 0x6367ae56d3fd535d686c06bf8dfcb7,
//                 0x8ad61,
//             ],
//         ),
//         c1: Fp2::new(
//             [
//                 0xae4a11560fd2b96a09231943c6dab,
//                 0xd18b53078ed5d3f5cadea0c65f2734,
//                 0xff3ed67578ddeae5439da2b02eab31,
//                 0x8ccec,
//             ],
//             [
//                 0x38f87f1c9644538b4bf7ec8dee1be9,
//                 0x5fb2ff8fddcf6fbb096adcc2b04bfc,
//                 0x56170fab6f325e416abe497d9e055d,
//                 0xbd355,
//             ],
//         ),
//         c2: Fp2::new(
//             [
//                 0x4579a52cfa8b4f01e3fc00c5962be5,
//                 0x45b3857ffd05157321b22d5a04c013,
//                 0xb7e3b626ce947b122db3425e68ae97,
//                 0xfc9d8,
//             ],
//             [
//                 0x8fd0966dcadba0b5f86e4c5ffa1da2,
//                 0x18864f9cef0cf179640918b4d31a73,
//                 0xf2af61549b350ae4025242be1008a6,
//                 0x325f3,
//             ],
//         ),
//     };

//     let b = Fp6 {
//         c0: Fp2::new(
//             [
//                 0x71825b46349b307298025eea77770c,
//                 0x225cbdc399bb960534036f625a3592,
//                 0xb32ff3d468eb1aba35968a2bbdb9d0,
//                 0x66e47,
//             ],
//             [
//                 0x3fc9a610ece44f3c738a358aa60150,
//                 0x602cf426b135e0e0d837a8cf5d9ec7,
//                 0x4aee9cb5468ad40ee5dccbf0264407,
//                 0xdc722,
//             ],
//         ),
//         c1: Fp2::new(
//             [
//                 0xb77accf6a4e5e337a85d2a62baf070,
//                 0x5e5a337d7db6344484087a18586c7c,
//                 0xb3f92a8e29eabd22358787933cd419,
//                 0xa4aa,
//             ],
//             [
//                 0xb33c802bede6916fa9f3d1d08a643b,
//                 0x4537aef2966b539d7364999eb35707,
//                 0xde1e413ba79349fbe6cc88298b9f18,
//                 0x1776be,
//             ],
//         ),
//         c2: Fp2::new(
//             [
//                 0x8a449468e00f21f8f78f53b57f2f81,
//                 0x618036c76149747cdb368c1d1ccf7a,
//                 0x55fe3e909dbe81264b06f976721808,
//                 0x122bf8,
//             ],
//             [
//                 0x157c33763f6eb6ec8697a38c6cb1d2,
//                 0xbae7a1aac9a8926b246110708adc2b,
//                 0x8c3ab0276bcf5791fb656da17e25c6,
//                 0xe008d,
//             ],
//         ),
//     };

//     let c = Fp6 {
//         c0: Fp2::new(
//             [
//                 0x73077f94bdbbac2eb3081372118ec2,
//                 0xebd1f3f534eff775c93619ee45d822,
//                 0x9422703b2b18bd664b85022f395f19,
//                 0xe2dbc,
//             ],
//             [
//                 0x1c91d6a1211c55b5984826c4b87dea,
//                 0xaf48b3d291c8c05191ac1195e7596a,
//                 0x2835f2efdacbba250c88683ca4ded0,
//                 0x8e2df,
//             ],
//         ),
//         c1: Fp2::new(
//             [
//                 0xa2fabdedc6584f17cbe1861e5576ef,
//                 0x9bc6801a1ff2d07dccc3361506b1a9,
//                 0xe2d66551040faf17b455385b218482,
//                 0x5181a,
//             ],
//             [
//                 0x9520077c93f6e531930b6d8bf16d4f,
//                 0xc6acd448b87eb8df27705139c33b6d,
//                 0xdc994b13d737d3921fed971d068992,
//                 0x9ac79,
//             ],
//         ),
//         c2: Fp2::new(
//             [
//                 0x157c33763f6eb6ec8697a38c6cb1d2,
//                 0xbae7a1aac9a8926b246110708adc2b,
//                 0x8c3ab0276bcf5791fb656da17e25c6,
//                 0xe008d,
//             ],
//             [
//                 0x38f87f1c9644538b4bf7ec8dee1be9,
//                 0x5fb2ff8fddcf6fbb096adcc2b04bfc,
//                 0x56170fab6f325e416abe497d9e055d,
//                 0xbd355,
//             ],
//         ),
//     };

//     let d = Fp6 {
//         c0: Fp2::new(
//             [
//                 0x1bf849df18fa4b181e5ccd88d00322,
//                 0xa873e3f09b7c2381a8fb576607f08,
//                 0xb93bad407b1a2137e9139fd1d02701,
//                 0x11fec3,
//             ],
//             [
//                 0x6fc4a038fe6225c716e715a8693747,
//                 0x42cca14b682ed6eadd856234d489fe,
//                 0xd0bb60bb40be57f29e4f5c7bae0c1b,
//                 0x1bdb31,
//             ],
//         ),
//         c1: Fp2::new(
//             [
//                 0xfd26131772a6aacc89a8241f7251f,
//                 0xc45e17aa2cd060e622199ec81811b6,
//                 0x4ea2f446bcd8021699feb458b9f171,
//                 0x94424,
//             ],
//             [
//                 0x6d68ddd1cd3385ae17b186d969614,
//                 0x62f219bc2846c05b32bfa068733b7,
//                 0x22b83f770d5e2099c98826a54e7567,
//                 0x99d2a,
//             ],
//         ),
//         c2: Fp2::new(
//             [
//                 0x949ec32b351b88e631bd6142d69083,
//                 0x1a207b24276583ac59736aed91e250,
//                 0xaa69afa41da815ae5115132007b458,
//                 0x4baae,
//             ],
//             [
//                 0xe1c4dc1304bf3b3e48eb83a8b5d3d3,
//                 0xe41bccfa304778371063dbcad7d778,
//                 0x7c6b7dca8744013ed4503ba6c430f3,
//                 0x17e9f5,
//             ],
//         ),
//     };

//     let a_square = a.square();
//     let b_square = b.square();
//     let c_square = c.square();

//     let square_mul_a = a.mul_interleaved(a);
//     let square_mul_b = b.mul_interleaved(b);
//     let square_mul_c = c.mul_interleaved(c);
//     assert(a_square.eq(square_mul_a));
//     assert(b_square.eq(square_mul_b));
//     assert(c_square.eq(square_mul_c));

//     // (a + b) * c.square(), (c * c * a) + (c * c * b)
//     let lhs = (a.add(b)).mul_interleaved(c.square());
//     let rhs =
//         (c.mul_interleaved(c).mul_interleaved(a)).add(c.mul_interleaved(c).mul_interleaved(b));

//     assert(lhs.eq(rhs));

//     // (a+b) * d = a*d + b*d
//     let lhs1 = (a.add(b)).mul(d);
//     let rhs1 = (a.mul(d)).add(b.mul(d));
//     assert(lhs1.eq(rhs1));
// }

#[test]
fn test_square() {
    let a = Fp6::new(
        Fp2::new(
            [281382237812545379193668594878174908, 90824900576507601223966148779479394, 1640],
            [1199704082799405279767970336007399382, 105838412715217748738343588203455976, 3470],
        ),
        Fp2::new(
            [699517922600938094581844456622821606, 766244793208285317105753995770523332, 5391],
            [1174771148458207518820683564604158404, 890744614640420316156493451275032886, 4653],
        ),
        Fp2::new(
            [489655109890455545367468964680482035, 398027084094548452922359112094161148, 6728],
            [496665123426422059639752844234206160, 350210767051572081629044293495878811, 1877],
        ),
    );
    println("a^2 is");
    println(a.square());
    println("a^2 should be ");
    println(a * a);
    let b = a.square(); 
    assert(b.eq(a * a));
}


#[test]
fn test_invert() {
    let a = Fp6::new(
        Fp2::new(
            [281382237812545379193668594878174908, 90824900576507601223966148779479394, 1640],
            [1199704082799405279767970336007399382, 105838412715217748738343588203455976, 3470],
        ),
        Fp2::new(
            [699517922600938094581844456622821606, 766244793208285317105753995770523332, 5391],
            [1174771148458207518820683564604158404, 890744614640420316156493451275032886, 4653],
        ),
        Fp2::new(
            [489655109890455545367468964680482035, 398027084094548452922359112094161148, 6728],
            [496665123426422059639752844234206160, 350210767051572081629044293495878811, 1877],
        ),
    );
    let b = a.invert(); 
    assert((a*b).eq(Fp6::one()));
}
