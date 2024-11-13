use crate::fp2::Fp2;
use crate::fp6::Fp6;
use crate::fp12::Fp12;
use crate::g1::{G1Affine};
use crate::g2::{G2Affine, G2Projective};

// Ref: https://github.com/zkcrypto/bls12_381/blob/main/src/pairings.rs

// https://github.com/zkcrypto/bls12_381/blob/main/src/pairings.rs#L709
pub(crate) fn doubling_step(r: G2Projective) -> (Fp2, Fp2, Fp2, G2Projective) {
    let mut tmp0 = r.x.square();
    let mut tmp1 = r.y.square();
    let mut tmp2 = tmp1.square();
    let mut tmp3 = ((tmp1 + r.x).square()).sub(tmp0).sub(tmp2);
    tmp3 = tmp3.add(tmp3);
    let tmp4 = tmp0.add(tmp0).add(tmp0);
    let mut tmp6 = r.x.add(tmp4);
    let tmp5 = tmp4.square();
    let zsquared = r.z.square();
    let new_r_x = tmp5.sub(tmp3).sub(tmp3);
    let new_r_z = ((r.z.add(r.y)).square()).sub(tmp1).sub(zsquared);
    let mut new_r_y = (tmp3.sub(new_r_x)).mul(tmp4);
    tmp2 = tmp2.add(tmp2);
    tmp2 = tmp2.add(tmp2);
    tmp2 = tmp2.add(tmp2);
    new_r_y = new_r_y.sub(tmp2);
    tmp3 = tmp4.mul(zsquared);
    tmp3 = tmp3.add(tmp3);
    tmp3 = tmp3.neg();
    tmp6 = (tmp6.square()).sub(tmp0).sub(tmp5);
    tmp1 = tmp1.add(tmp1);
    tmp1 = tmp1.add(tmp1);
    tmp6 = tmp6.sub(tmp1);
    tmp0 = new_r_z.mul(zsquared);
    tmp0 = tmp0.add(tmp0);

    (tmp0, tmp3, tmp6, G2Projective { x: new_r_x, y: new_r_y, z: new_r_z })
}

pub(crate) fn addition_step(r: G2Projective, q: G2Affine) -> (Fp2, Fp2, Fp2, G2Projective) {
    let zsquared = r.z.square();
    let ysquared = q.y.square();
    let mut t0 = zsquared.mul(q.x);
    let mut t1 = (((q.y.add(r.z)).square()).sub(ysquared).sub(zsquared)).mul(zsquared);
    let t2 = t0.sub(r.x);
    let t3 = t2.square();

    let mut t4 = t3.add(t3);
    t4 = t4.add(t4);
    let t5 = t4.mul(t2);
    let mut t6 = t1.sub(r.y).sub(r.y);
    let mut t9 = t6.mul(q.x);
    let t7 = t4.mul(r.x);

    let new_r_x = (t6.square()).sub(t5).sub(t7).sub(t7);

    let new_r_z = ((r.z.add(t2)).square()).sub(zsquared).sub(t3);

    let mut t10 = q.y.add(new_r_z);
    let t8 = (t7.sub(new_r_x)).mul(t6);
    t0 = r.y.mul(t5);
    t0 = t0.add(t0);
    let new_r_y = t8.sub(t0);

    t10 = t10.square().sub(ysquared);
    let ztsquared = new_r_z.square();
    t10 = t10.sub(ztsquared);
    t9 = t9.add(t9).sub(t10);

    let t10_double = new_r_z.add(new_r_z);
    let t6_neg = t6.neg();
    let t1_final = t6_neg.add(t6_neg);
    (t10_double, t1_final, t9, G2Projective { x: new_r_x, y: new_r_y, z: new_r_z })
}

pub(crate) fn ell(f: Fp12, coeffs: (Fp2, Fp2, Fp2), p: G1Affine) -> Fp12 {
    let c0c0 = coeffs.0.c0 * p.y;
    let c0c1 = coeffs.0.c1 * p.y;

    let c1c0 = coeffs.1.c0 * p.x;
    let c1c1 = coeffs.1.c1 * p.x;

    f.mul_by_014(
        coeffs.2,
        Fp2 { c0: c1c0, c1: c1c1 },
        Fp2 { c0: c0c0, c1: c0c1 }
    )
}

pub fn miller_loop(p: G1Affine, q: G2Affine) -> Fp12 {
    // https://github.com/zkcrypto/bls12_381/blob/main/src/pairings.rs#L672
    // Got array from zkcrypto impl
    let booleans = [
        false, true, true, false, true, false, false, true, false, false, false, false, false, false, false, false,
        true, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false,
        false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false,
        true, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false
    ];

    let mut adder_cur = G2Projective::from_affine(q);
    let mut adder_base = q;
    let mut adder_p = p;
    let mut f = Fp12::one();
    let mut found_one = false;

    for i in 0..64 {
        let i_bool = booleans[i];
        if !found_one {
            found_one = i_bool;
        } else {
            // Doubling step receives f. Does 2 things: doubling_step function & ell function
            let (coeffs_c0, coeffs_c1, coeffs_c2, cur_updated) = doubling_step(adder_cur);
            adder_cur = cur_updated;
            // ell updates f
            f = ell(f, (coeffs_c0, coeffs_c1, coeffs_c2), adder_p);

            if i_bool {
                // addition_step receives f & does 2 things: addition_step & ell function
                let (coeffs_c0, coeffs_c1, coeffs_c2, cur_updated) = addition_step(adder_cur, adder_base);
                adder_cur = cur_updated;
                // ell updates f
                f = ell(f, (coeffs_c0, coeffs_c1, coeffs_c2), adder_p);
            }
            f = f.square();
        }
    }
    // Doubling step receives f. Does 2 things: doubling_step function & ell function
    let (coeffs_c0, coeffs_c1, coeffs_c2, cur_updated) = doubling_step(adder_cur);
    adder_cur = cur_updated;
    // ell updates f
    f = ell(f, (coeffs_c0, coeffs_c1, coeffs_c2), adder_p);

    f = f.conjugate();
    f
}

// https://github.com/zkcrypto/bls12_381/blob/main/src/pairings.rs#L50
fn fp4_square(a: Fp2, b: Fp2) -> (Fp2, Fp2) {
    let t0 = a.square();
    let t1 = b.square();
    let mut t2 = t1.mul_by_nonresidue();
    let c0 = t2.add(t0);
    t2 = a.add(b);
    t2 = t2.square();
    t2 = t2.sub(t0);
    let c1 = t2.sub(t1);

    (c0, c1)
}

// https://github.com/zkcrypto/bls12_381/blob/main/src/pairings.rs#L66
fn cyclotomic_square(f: Fp12) -> Fp12 {
    let mut z0 = f.c0.c0;
    let mut z4 = f.c0.c1;
    let mut z3 = f.c0.c2;
    let mut z2 = f.c1.c0;
    let mut z1 = f.c1.c1;
    let mut z5 = f.c1.c2;

    let (t0, t1) = fp4_square(z0, z1);

    z0 = t0.sub(z0);
    z0 = z0.add(z0).add(t0);
    z1 = t1.add(z1);
    z1 = z1.add(z1).add(t1);

    let (mut t0, t1) = fp4_square(z2, z3);
    let (t2, t3) = fp4_square(z4, z5);

    z4 = t0.sub(z4);
    z4 = z4.add(z4).add(t0);

    z5 = t1.add(z5);
    z5 = z5.add(z5).add(t1);

    t0 = t3.mul_by_nonresidue();
    z2 = t0.add(z2);
    z2 = z2.add(z2).add(t0);

    z3 = t2.sub(z3);
    z3 = z3.add(z3).add(t2);

    Fp12 { c0: Fp6 { c0: z0, c1: z4, c2: z3 }, c1: Fp6 { c0: z2, c1: z1, c2: z5 } }
}

fn cycolotomic_exp(f: Fp12) -> Fp12 {
    // https://github.com/zkcrypto/bls12_381/blob/main/src/pairings.rs#L119
    // Got array from zkcrypto impl
    let booleans = [
        true, true, false, true, false, false, true, false, false, false, false, false, false, false, false, true,
        false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false,
        false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, true,
        false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false
    ];

    let mut tmp = Fp12::one();
    let mut found_one = false;
    for i in 0..64 {
        if found_one {
            tmp = cyclotomic_square(tmp);
        } else {
            found_one = booleans[i];
        }

        if booleans[i] {
            tmp = tmp.mul(f);
        }
    }

    tmp.conjugate()
}

impl Fp12 {
    pub(crate) fn final_exponentiation(self: Self) -> Self {
        let mut f = self;
        let mut t0 = f.frobenius_map().frobenius_map().frobenius_map().frobenius_map().frobenius_map().frobenius_map();
        // This will error if inverse doesn't exist
        let mut t1 = f.invert();
        let mut t2 = t0.mul(t1);
        t1 = t2;
        t2 = t2.frobenius_map().frobenius_map();
        t2 = t2.mul(t1);
        t1 = cyclotomic_square(t2).conjugate();
        let mut t3 = cycolotomic_exp(t2);
        let mut t4 = cyclotomic_square(t3);
        let mut t5 = t1.mul(t3);
        t1 = cycolotomic_exp(t5);
        t0 = cycolotomic_exp(t1);
        let mut t6 = cycolotomic_exp(t0);
        t6 = t6.mul(t4);
        t4 = cycolotomic_exp(t6);
        t5 = t5.conjugate();
        t4 = t4.mul(t5).mul(t2);
        t5 = t2.conjugate();
        t1 = t1.mul(t2);
        t1 = t1.frobenius_map().frobenius_map().frobenius_map();
        t6 = t6.mul(t5);
        t6 = t6.frobenius_map();
        t3 = t3.mul(t0);
        t3 = t3.frobenius_map().frobenius_map();
        t3 = t3.mul(t1);
        t3 = t3.mul(t6);
        f = t3.mul(t4);
        f
    }
}

pub fn pairing(p: G1Affine, q: G2Affine) -> Fp12 {
    // If any of the inputs is infinity, the result is 1
    if p.infinity | q.infinity {
        Fp12::one()
    } else {
        let tmp = miller_loop(p, q);
        tmp.final_exponentiation()
    }
}

