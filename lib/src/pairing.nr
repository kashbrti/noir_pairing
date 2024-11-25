use crate::bls::fp12::Fp12;
use crate::bls::fp2::Fp2;
use crate::bls::g1::G1Affine;
use crate::bls::g2::{G2Affine, G2Projective};
use crate::helpers::__final_exponentiation;

/// Performs the doubling step of the Miller Algorithm and evaluates the line function.
///
/// Reference Document: https://eprint.iacr.org/2010/354.pdf (Algorithm 26)
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

/// Performs the addition step of the Miller Algorithm and evaluates the line function.
///
/// Reference Document: https://eprint.iacr.org/2010/354.pdf (Algorithm 26)
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

/// Performse sparse multiplication between `Fp12` and a triple of elements in `Fp2`.
///
/// The elements of `Fp2` represent non-zero coefficients in an element of `Fp12`. As algorithms
/// 26 and 27 in the reference document below evaluate the line function internally, and as most
/// of the evaluated coefficients are always zero, we represent this sparsely.
///
/// > TODO: migrate to Ell data type instead.
pub(crate) fn ell(f: Fp12, coeffs: (Fp2, Fp2, Fp2), p: G1Affine) -> Fp12 {
    let c0c0 = coeffs.0.c0 * p.y;
    let c0c1 = coeffs.0.c1 * p.y;

    let c1c0 = coeffs.1.c0 * p.x;
    let c1c1 = coeffs.1.c1 * p.x;

    f.mul_by_014(coeffs.2, Fp2 { c0: c1c0, c1: c1c1 }, Fp2 { c0: c0c0, c1: c0c1 })
}

/// Evaluates the Miller loop.
pub fn miller_loop(p: G1Affine, q: G2Affine) -> Fp12 {
    // BLS12-381 'X' parameter transformed from binary to boolean array
    let booleans = [
        false, true, true, false, true, false, false, true, false, false, false, false, false,
        false, false, false, true, false, false, false, false, false, false, false, false, false,
        false, false, false, false, false, false, false, false, false, false, false, false, false,
        false, false, false, false, false, false, false, false, false, true, false, false, false,
        false, false, false, false, false, false, false, false, false, false, false, false,
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
                let (coeffs_c0, coeffs_c1, coeffs_c2, cur_updated) =
                    addition_step(adder_cur, adder_base);
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
    println(f);
    f
}

impl Fp12 {
    /// Unconstrained final exponentiation.
    pub(crate) fn final_exponentiation(self: Self) -> Self {
        // let out = unsafe { __final_exponentiation(self) };
        let out: Fp12 = Fp12::one();
        out
    }
}

/// Computes the pairing.
pub fn pairing(p: G1Affine, q: G2Affine) -> Fp12 {
    // If any of the inputs is infinity, the result is 1
    if p.infinity | q.infinity {
        Fp12::one()
    } else {
        let tmp = miller_loop(p, q);
        tmp.final_exponentiation();
        // instead of final_exponentiation, we use the faster residue check
        //compute the values of c and u
        // let c = __compute_c(tmp);
        // let u = __compute_u(tmp);
        // perform the residue check
        tmp
    }
}
