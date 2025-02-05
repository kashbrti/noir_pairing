use dep::pairing_lib::bls::g1::G1Affine;
use dep::pairing_lib::bls::g2::G2Affine;
use dep::pairing_lib::bls::pairing::pairing;

fn main(p: G1Affine, q: G2Affine) {
    let res = pairing(p, q);
}
