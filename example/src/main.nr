use dep::pairing::pairing::{G1Affine, G2Affine, pairing};

fn main(p: G1Affine, q: G2Affine) {
    let res = pairing(p, q);
}
