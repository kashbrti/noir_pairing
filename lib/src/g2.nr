use crate::fp2::{Fp2};

// Ref: https://github.com/zkcrypto/bls12_381/blob/main/src/pairings.rs

pub struct G2Affine {
  pub x: Fp2,
  pub y: Fp2,
  pub infinity: bool,
}

pub struct G2Projective {
  pub x: Fp2,
  pub y: Fp2,
  pub z: Fp2,
}

impl G2Affine {
    pub fn neg(self: Self) -> Self {
        let mut y = self.y.neg();
        if self.infinity {
            y = Fp2::one();
        }
        G2Affine { x: self.x, y, infinity: self.infinity }
    }

    pub fn generator() -> Self {
        G2Affine {
            x: Fp2::new(
                [
                0xac0326a805bbefd48056c8c121bdb8,
                0x7ad4fa403b02b4510b647ae3d1770b,
                0xb2f08f0a91260805272dc51051c6e4,
                0x24aa2
            ],
                [
                0x4cf11213945d57e5ac7d055d042b7e,
                0xd0d09920b61ab5da61bbdc7f504933,
                0x6052719f607dacd3a088274f65596b,
                0x13e02b
            ]
            ),
            y: Fp2::new(
                [
                0x3ac9cc3baca289e193548608b82801,
                0x9baa8cbdd3a76d429a695160d12c92,
                0x27727d6e118cc9cdc6da2e351aadfd,
                0xce5d5
            ],
                [
                0x370d275cec1da1aaa9075ff05f79be,
                0x287e85a763af267492ab572e99ab3f,
                0xa02ea734cc32acd2b02bc28b99cb3e,
                0x606c4
            ]
            ),
            infinity: false
        }
    }
}

impl G2Projective {
    pub fn from_affine(p: G2Affine) -> Self {
        let mut z = Fp2::one();
        if (p.infinity) {
            z = Fp2::zero();
        }
        G2Projective { x: p.x, y: p.y, z }
    }
}
