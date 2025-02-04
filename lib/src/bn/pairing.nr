use crate::bn::fp2::Fp2;
use crate::bn::g1::G1Affine;
use crate::bn::g2::G2Affine;
use crate::bn::g2::G2Projective;
use crate::bn::fp12::Fp12;
use crate::bn::constants::NAF_DIGIT; 
use crate::helpers::final_exponentiation;
use crate::helpers::__get_pairing_witnesses;

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

pub(crate) impl Ell {
    /// Constructs an Ell instance with coefficients in Fp2.
    pub(crate) fn new(c0_c0: Fp2, c1_c0: Fp2, c1_c1: Fp2) -> Self {
        Self { c0_c0, c1_c0, c1_c1 }
    }

    /// Scales the coefficients of an Ell instance by an affine point.
    pub(crate) fn scale_by_affine(self, p: G1Affine) -> Self {
        Self { c0_c0: self.c0_c0, c1_c0: self.c1_c0.scale(p.x), c1_c1: self.c1_c1.scale(p.y) }
    }
}


pub(crate) fn addition_step(r: G2Projective, q: G2Affine) -> (Ell, G2Projective) {
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
    (Ell::new(t9, t1_final, t10_double), G2Projective { x: new_r_x, y: new_r_y, z: new_r_z })
}

// these functions are the same for BLS and BN, we should move this outside the specific type and use generics with the extension field trait
pub(crate) fn doubling_step(r: G2Projective) -> (Ell, G2Projective) {
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

    (Ell::new(tmp6, tmp3, tmp0), G2Projective { x: new_r_x, y: new_r_y, z: new_r_z })
}



/// Evaluates the Miller loop.
pub fn miller_loop(p: G1Affine, q: G2Affine) -> Fp12 {
    // BLS12-381 'X' parameter transformed from binary to boolean array
    // let booleans = [
    //     false, true, true, false, true, false, false, true, false, false, false, false, false,
    //     false, false, false, true, false, false, false, false, false, false, false, false, false,
    //     false, false, false, false, false, false, false, false, false, false, false, false, false,
    //     false, false, false, false, false, false, false, false, false, true, false, false, false,
    //     false, false, false, false, false, false, false, false, false, false, false, false,
    // ];

    let mut adder_cur = G2Projective::from_affine(q);
    let mut adder_base = q;
    let mut adder_p = p;
    let mut f = Fp12::one();
    let mut found_one = false;

    for i in 0..66 {
        let i_NAF = NAF_DIGIT[i];
        if !found_one {
            found_one = (i_NAF != 0);
        } else {
            // Doubling step receives f. Does 2 things: doubling_step function & ell function
            let (ell_coeffs, cur_updated) = doubling_step(adder_cur);

            adder_cur = cur_updated;

            // updates f
            f = f.sparse_mul(ell_coeffs.scale_by_affine(adder_p));

            if i_NAF == 1 {
                // addition_step receives f & does 2 things: addition_step & ell function
                // negate the G2 point before adding
                let (ell_coeffs, cur_updated) = addition_step(adder_cur, adder_base);

                adder_cur = cur_updated;
                // ell updates f
                f = f.sparse_mul(ell_coeffs.scale_by_affine(adder_p));
            }
            else if i_NAF == 2 {

                let (ell_coeffs, cur_updated) = addition_step(adder_cur, adder_base.neg());

                adder_cur = cur_updated;
                // ell updates f
                f = f.sparse_mul(ell_coeffs.scale_by_affine(adder_p));
            }

            f = f.square();
        }
    }

    // Doubling step receives f. Does 2 things: doubling_step function & ell function
    let (ell_coeffs, cur_updated) = doubling_step(adder_cur);
    adder_cur = cur_updated;
    // ell updates f
    f = f.sparse_mul(ell_coeffs.scale_by_affine(adder_p));

    f = f.conjugate();
    println(f);
    f
}



// now we get to the main function 
pub fn pairing(p: G1Affine, q: G2Affine) -> Fp12 {
    let f = miller_loop(p, q); 
    final_exponentiation(f)
}

// the pairing equality check 
pub fn pairing_equality(p: G1Affine, q: G2Affine, r: G1Affine, s: G2Affine) -> bool {
    // we have to check that there exists c and u such that  f_(r,Q0)(P0)* f_(r,Q1)(-P1) = c^\lambda * u and u is a d^i-th root of unity.
    let f = miller_loop(p, q); 
    let g = miller_loop(r.neg(), s); 
    let miller_output = f*g; 
    let (c,u) = unsafe{__get_pairing_witnesses(miller_output)}; 
    // now we do the verification check on the witnesses 
    let witness_res = c.exp_lambda()*u; 
    witness_res == miller_output
}


#[test]
fn test_miller_loop_ghe() {
    let p = G1Affine::generator(); 
    let q = G2Affine::generator(); 
    let f = miller_loop(p, q); 
    println(f); 
}

