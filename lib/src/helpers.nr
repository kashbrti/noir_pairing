use crate::bls::g1::G1Affine;
use crate::bls::g2::{G2Affine, G2Projective};
use crate::bn::exp_mod::ExpModParams;
use crate::bn::fp12::Fp12;
use crate::bn::fp2::Fp2;
use crate::bn::fp6::Fp6;
use crate::bn::mod_p12m1::ModP12M1Params;
use bignum::BigNum;
use bignum::fields::bn254Fq::BN254_Fq_Params;
use bignum::fns::unconstrained_ops::__pow;

// pub type ExponentsType = BigNum<4, 381>;
// TODO: we need bignum params for the exponents mod p^12-1
// we should add an exponentiation function that takes an Fp12 element and a bignum and does the exp using fancy frobenius
// we need 3rd root and r-th root algorithms
// we need to compute the value of lambda = \sum_i a_i p^i such that a_i is smooth
pub type Fp = BigNum<3, 254, BN254_Fq_Params>;
// the type of integers modulo p^12-1 /r
pub type Mod_h = BigNum<24, 2790, ExpModParams>;
// the type of integers modulo p^12 -1
pub type Mod_p12m1 = BigNum<26, 3044, ModP12M1Params>;

/// Unconstrained final exponentiation

pub(crate) fn final_exponentiation(input: Fp12) -> Fp12 {
    let mut f = input;
    let mut t0 = f
        .frobenius_map()
        .frobenius_map()
        .frobenius_map()
        .frobenius_map()
        .frobenius_map()
        .frobenius_map();
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


pub unconstrained fn __get_pairing_witnesses(a: Fp12)-> (Fp12, Fp12){
    let (c_coeffs,u_coeffs) = __get_pairing_witnesses_oracle(a);
    let c = Fp12::from_coeff(c_coeffs);
    let u = Fp12::from_coeff(u_coeffs);
    (c,u)
}

#[oracle(get_pairing_witnesses)]
pub unconstrained fn __get_pairing_witnesses_oracle(a: Fp12)-> ([[Field; 3]; 12], [[Field; 3]; 12]){
}

// /// Unconstrained `Fp12::exponentiation`
// pub (crate) unconstrained fn __fp12_exponentiation<let N: u32, let params:BigNumParams>(input: Fp12, exponent: [Field; N]) -> Fp12 {
//     let mut res = Fp12::one();
//     for i in 0..N {
//         if exponent[i] == 1 {
//             res = __fp12_square(res);
//         }
//         res = __fp12_mul(res,input);
//     }
//     res
// }

// pub(crate) fn res_reciprocity(c: Field, u: Fp12, miller_loop_output: Fp12) {
//     // we have to check that miller_loop_output = c^lambda u
//     assert(true);
// }

/// Computes the square of an element in `Fp4` composed of elements in `Fp2`.
///
/// Reference Document: https://eprint.iacr.org/2022/367.pdf (Algorithm 9)
fn fp4_square(a: Fp2, b: Fp2) -> (Fp2, Fp2) {
    let t0 = a.square();
    let t1 = b.square();
    let mut t2 = t1.mul_by_non_residue();
    let c0 = t2.add(t0);
    t2 = a.add(b);
    t2 = t2.square();
    t2 = t2.sub(t0);
    let c1 = t2.sub(t1);

    (c0, c1)
}



/// Efficient `Fp12` square in the cyclotomic subgroup $C_{\phi^6}$
///
/// Reference Document: https://www.math.u-bordeaux.fr/~damienrobert/csi/book/book.pdf
/// (Algorithm 5.5.4)
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

    t0 = t3.mul_by_non_residue();
    z2 = t0.add(z2);
    z2 = z2.add(z2).add(t0);

    z3 = t2.sub(z3);
    z3 = z3.add(z3).add(t2);

    Fp12 { c0: Fp6 { c0: z0, c1: z4, c2: z3 }, c1: Fp6 { c0: z2, c1: z1, c2: z5 } }
}

/// Returns the cyclotomic exponentiation of an element in `Fp12`.
///
/// We use a square-and-multiply algorithm for exponentiation.
fn cycolotomic_exp(f: Fp12) -> Fp12 {
    // BLS12-381 'X' parameter transformed from binary to boolean array
    let booleans = [
        true, true, false, true, false, false, true, false, false, false, false, false, false,
        false, false, true, false, false, false, false, false, false, false, false, false, false,
        false, false, false, false, false, false, false, false, false, false, false, false, false,
        false, false, false, false, false, false, false, false, true, false, false, false, false,
        false, false, false, false, false, false, false, false, false, false, false, false,
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

// #[test]
// fn test_fp12_exponentiation(){
//     let input = Fp12::one();
//     let exponent = [1,0,1,0,1,0];
//     let res = __fp12_exponentiation::<6>(input,exponent);
//     assert(res == input);
// }


#[test]
fn test_witness_generator(){
    let a = Fp12::one();
    let (c,u) = unsafe{__get_pairing_witnesses(a)};
    println("c: ");
    println(c);
    println("u: ");
    println(u);
}

