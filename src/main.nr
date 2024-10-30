mod fp2;

use fp2::Fp2;
use bignum::BigNum;
use noir_bigcurve::curves::bls12_381::BLS12_381Fq;

fn main(x: Fp2, y: pub Fp2) {
    // assert(x != y);
}

#[test]
fn test_main() {
    let x: BLS12_381Fq = BigNum::new();
    let y: BLS12_381Fq = BigNum::new();
    let fp2_1 = Fp2 { c0: x, c1: y };
    let fp2_2 = Fp2 { c0: x, c1: y };
    main(fp2_1, fp2_2);
}
