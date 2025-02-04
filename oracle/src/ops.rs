use ark_ff::{BigInteger768, BitIteratorBE, Field, PrimeField}; // Add BitIteratorBE to the imports
                                                               // We'll use a field associated with the BLS12-381 pairing-friendly
                                                               // group for this example.
use ark_bn254::Fr;
use ark_bn254::{Fq, Fq12};
use ark_test_curves::bls12_381::Fq2 as F;
// use ark_grumpkin::Fq as Fr;
// use ark_test_curves::bn254::Fq12 as F2;
// `ark-std` is a utility crate that enables `arkworks` libraries
// to easily support `std` and `no_std` workloads, and also re-exports
// useful crates that should be common across the entire ecosystem, such as `rand`.
use ark_std::{One, UniformRand};
use num_bigint::BigUint;

use ark_ff::BigInt;
use num_bigint::BigInt as BI;


pub struct WitnessGenerator {
}

pub trait WitnessGeneratorTrait {
    fn witness_generator(input: Fq12) -> (Fq12, Fq12);
    fn is_third_root(a: &Fq12) -> bool;
    fn mp_th_root_of_c(c: Fq12) -> Fq12;
    fn r_th_root_of_f(f: Fq12) -> Fq12;
    fn pow_p12_minus_one_div_3(a: &Fq12) -> Fq12;
    fn pow_p12_minus_one_div_27(a: &Fq12) -> Fq12;
    fn find_27th_root() -> Fq12;
    fn find_third_non_residue() -> (Fq12, Fq12);
    fn get_order(a: Fq12) -> u32;
    fn tonelli_shanks_third_root(a: Fq12) -> Fq12;
    fn invert(a: &BigUint, modulus: &BigUint) -> BI;
    fn extended_gcd(_a: &BI, _b: &BI) -> (BI, BI, BI);
    fn rand_third_root() -> Fq12;
}

impl WitnessGeneratorTrait for WitnessGenerator {
    fn witness_generator(input: Fq12) -> (Fq12, Fq12) {
    let w = Self::find_27th_root();   
    let mut f = input;
    let s: u32 = if Self::pow_p12_minus_one_div_3(&f) == Fq12::one() {
        0
    } else if Self::pow_p12_minus_one_div_3(&(f * w)) == Fq12::one() {
        f = f * w;
        1
    } else {
        f = f * w.pow([2]);
        2
    };
    // set c to be the r-th root of f
    let mut c = Self::r_th_root_of_f(f);
    // take the m'-th root of c
    c = Self::mp_th_root_of_c(c);
    // take the third root of the result
    c = Self::tonelli_shanks_third_root(c);
    (c, w.pow([3^s as u64]))
}

 fn mp_th_root_of_c(c: Fq12) -> Fq12 {
    // m = lambda / r
    // lambda = 2 + q - q^2 + q^3 + 6*x
    // x is 65 bits long
    // x-1 /2 = 2482830683596424440
    let q: BigUint = Fq::MODULUS.into();
    let xm1div2 = BigUint::from(2482830683596424440u64);
    let x = xm1div2 * BigUint::from(2u64) + BigUint::one();
    let mut lambda = BigUint::from(6u64) * x + BigUint::from(2u64);
    lambda = lambda + &q;
    lambda = lambda + &q.pow(3);
    lambda = lambda - &q.pow(2);
    // mp = m/3
    let r: BigUint = Fr::MODULUS.into();
    let m = lambda / &r;
    // // mp = m/3
    let mp = m / BigUint::from(3u64);
    // // now we invert mp mod h, where h = q^12 - 1 / r
    let h = (q.pow(12) - BigUint::one()) / &r;
    let mpp: BigUint = Self::invert(&mp, &h).to_biguint().unwrap();
    // cast mpp to a bigint
    let mpp_bigint: BigInt<50> = BigInt::<50>::try_from(mpp).unwrap();
    let res = c.pow(mpp_bigint);
    res
}

 fn r_th_root_of_f(f: Fq12) -> Fq12 {
    // the value of r is
    //21888242871839275222246405745257275088548364400416034343698204186575808495617
    let r: BigUint = Fr::MODULUS.into();
    // the value of h is
    // p^12 - 1 / r
    let p = Fq::MODULUS;
    let p_biguint: BigUint = p.into();
    let p12 = p_biguint.pow(12);
    let h = (p12 - BigUint::one()) / &r;
    // now invert r mod h
    let inv_r: BI = Self::invert(&r, &h);
    let inv_r = inv_r.to_biguint().unwrap();
    let inv_r_bigint: BigInt<50> = BigInt::<50>::try_from(inv_r).unwrap();
    // compute the inverse of r mod h
    let res = f.pow(inv_r_bigint);
    res
}



 fn extended_gcd(_a: &BI, _b: &BI) -> (BI, BI, BI) {
    let (mut x, mut y) = (BI::from(0), BI::from(1));
    let (mut u, mut v) = (BI::from(1), BI::from(0));
    let mut a = _a.clone();
    let mut b = _b.clone();
    while a != BI::from(0) {
        let q = &b / &a;
        let r = &b % &a;
        let m = &x - &u * &q;
        let n = &y - &v * &q;
        b = a.clone();
        a = r;
        x = u;
        y = v;
        u = m;
        v = n;
    }
    (b.clone(), x.clone(), y.clone())
}

 fn invert(a: &BigUint, modulus: &BigUint) -> BI {
    // perform a -> BI conversion
    let a_bigint = BI::from(a.clone());
    let modulus_bigint = BI::from(modulus.clone());
    let (gcd, r, _) = Self::extended_gcd(&a_bigint, &modulus_bigint);
    assert_eq!(gcd, BI::from(1), "input and modulus are not coprime");
    let mut res = r.clone();
    while res < BI::from(0) {
        res = res + &modulus_bigint;
    }
    res
}

 fn tonelli_shanks_third_root(a: Fq12) -> Fq12 {
    let w = Self::find_27th_root();
    // s = p^12 -1 / 27
    let modulus: BigUint = Fq::MODULUS.into();
    let s = (modulus.pow(12u32) - BigUint::one()) / BigUint::from(27u64);
    let exp_biguint: BigUint = (s + BigUint::one()) / BigUint::from(3u64);
    let exp: BigInt<50> = BigInt::<50>::try_from(exp_biguint).unwrap();
    let mut x = a.pow(exp);
    while x.pow([3]) * (a.inverse().unwrap()) != Fq12::one() {
        x = x * w;
    }
    x
}

 fn get_order(a: Fq12) -> u32 {
    let mut a = a;
    let mut t = 0;
    while a != Fq12::one() {
        t += 1;
        a = a.pow(&[3, 0, 0, 0]);
    }

    return t;
}

 fn find_27th_root() -> Fq12 {
    let one = Fq12::one();
    let mut rng = ark_std::test_rng();
    loop {
        let x = Fq12::rand(&mut rng);
        let w = Self::pow_p12_minus_one_div_27(&x);
        if w != one && w.pow([9]) != one {
            return w;
        }
    }
}

 fn find_third_non_residue() -> (Fq12, Fq12) {
    loop {
        let mut rng = ark_std::test_rng();
        let a = Fq12::rand(&mut rng);
        let b = Self::pow_p12_minus_one_div_27(&a);
        if b.pow([3]) != Fq12::one() {
            return (a, b);
        }
    }
}

 fn pow_p12_minus_one_div_27(a: &Fq12) -> Fq12 {
    let p = Fq::MODULUS;
    let p_biguint: BigUint = p.into();
    let p12 = p_biguint.pow(12);
    let p12m1div27 = (p12 - BigUint::one()) / BigUint::from(27u64);
    // move this value back to ff_bigint
    let p12m1div27_bigint: BigInt<48> = BigInt::<48>::try_from(p12m1div27).unwrap();
    let res = a.pow(p12m1div27_bigint);
    res
}

 fn pow_p12_minus_one_div_3(a: &Fq12) -> Fq12 {
    let p = Fq::MODULUS;
    let p_biguint: BigUint = p.into();
    let p12 = p_biguint.pow(12);
    let p12m1div3 = (p12 - BigUint::one()) / BigUint::from(3u64);
    // move this value back to ff_bigint
    let p12m1div3_bigint: BigInt<48> = BigInt::<48>::try_from(p12m1div3).unwrap();
    let res = a.pow(p12m1div3_bigint);
    res
}

fn is_third_root(a: &Fq12) -> bool {
    let res = Self::pow_p12_minus_one_div_3(a);
    res == Fq12::one()
}

fn rand_third_root() -> Fq12 {
    let mut rng = ark_std::test_rng();
    let mut a = Fq12::rand(&mut rng);
    while !Self::is_third_root(&a) {
        a = Fq12::rand(&mut rng);
    }
    a
}

}



#[test]
fn test_mp_th_root_of_c() {
    let mut rng = ark_std::test_rng();
    let a = Fq12::rand(&mut rng);
    let res = WitnessGenerator::mp_th_root_of_c(a);
    println!("res: {:?}", res);
}

#[test]
fn test_r_th_root_of_f() {
    let mut rng = ark_std::test_rng();
    let a = Fq12::rand(&mut rng);
    let res = WitnessGenerator::r_th_root_of_f(a);
    println!("res: {:?}", res);
}



#[test]
fn test_pow_p12_minus_one_div_3() {
    let mut rng = ark_std::test_rng();
    let a = Fq12::rand(&mut rng);
    let res = WitnessGenerator::pow_p12_minus_one_div_3(&a);
    assert_eq!(res.pow([3]), Fq12::one());
}

#[test]
fn test_pow_p12_minus_one_div_27() {
    let mut rng = ark_std::test_rng();
    let a = Fq12::rand(&mut rng);
    let res = WitnessGenerator::pow_p12_minus_one_div_27(&a);
    // now we exponentiate by 27
    let twenty_seven = BigUint::from(27u64);
    let twenty_seven_bigint: BigInt<3> = BigInt::<3>::try_from(twenty_seven).unwrap();
    let res27 = res.pow(twenty_seven_bigint);
    assert_eq!(res27, Fq12::one());
}

#[test]
fn test_get_order() {
    let mut rng = ark_std::test_rng();
    let a = Fq12::rand(&mut rng);
    let b = WitnessGenerator::pow_p12_minus_one_div_27(&a);
    let res = WitnessGenerator::get_order(b);
}

#[test]
fn test_pow_fq12() {
    let mut rng = ark_std::test_rng();
    let a = Fq12::rand(&mut rng);
    let p = Fq::MODULUS;
    let p_biguint: BigUint = p.into();
    let p2 = p_biguint.pow(12);
    let p12m1 = p2 - BigUint::one();
    let p12m1_bigint = BigInt::<50>::try_from(p12m1).unwrap();
    let res = a.pow(p12m1_bigint);
    assert_eq!(res, Fq12::one());
}


#[test]
fn test_find_third_non_residue() {
    let res = WitnessGenerator::find_third_non_residue();
    let a = res.0;
    let b = res.1;
    let c = b.pow([3]);
    assert_ne!(c, Fq12::one());
    assert_eq!(c.pow([3]), Fq12::one());
}

#[test]
fn test_find_27th_root() {
    let res = WitnessGenerator::find_27th_root();
    assert_ne!(res.pow([9]), Fq12::one());
    assert_eq!(res.pow([27]), Fq12::one());
}

#[test]
fn test_tonelli_shanks_third_root() {
    let mut rng = ark_std::test_rng();
    let mut a = Fq12::rand(&mut rng);
    // check that a is a third root of unity
    // we do this by powering it up to (p^12-1)/3 and checking that it is 0
    while WitnessGenerator::pow_p12_minus_one_div_3(&a) != Fq12::one() {
        a = Fq12::rand(&mut rng);
    }
    let res = WitnessGenerator::tonelli_shanks_third_root(a);
    assert_eq!(res.pow([3]), a);
}

#[test]
fn test_invert() {
    // get a random field element
    let mut rng = ark_std::test_rng();
    let a = Fq::rand(&mut rng);
    // cast the field element to a biguint
    let a_biguint: BigUint = a.into();
    let modulus: BigUint = Fq::MODULUS.into();
    let res: BigUint = WitnessGenerator::invert(&a_biguint, &modulus).to_biguint().unwrap();
    let res_fq = Fq::from(res);
    assert_eq!(res_fq * a, Fq::one());
    assert_eq!(res_fq, a.inverse().unwrap());
}

#[test]
fn test_witness_generator() {
    let mut rng = ark_std::test_rng();
    let a = Fq12::rand(&mut rng);
    let res = WitnessGenerator::witness_generator(a);
    println!("res: {:?}", res);
}

#[test]
fn test_r_th_root_of_f_2() {
    let mut rng = ark_std::test_rng();
    let a = Fq12::rand(&mut rng);
    let res = WitnessGenerator::r_th_root_of_f(a);
    let r : BigUint = Fr::MODULUS.into();
    println!("res: {:?}", res);
    let res2 = res.pow(BigInt::<50>::try_from(r).unwrap());
    println!("res2: {:?}", res2);
}

// export the witness generator as a module
pub mod witness_generator {
    pub use super::WitnessGenerator;
    pub use super::WitnessGeneratorTrait;
}