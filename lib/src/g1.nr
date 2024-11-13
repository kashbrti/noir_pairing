use crate::fp2::{Fp};
use bignum::BigNum;

// Ref: https://github.com/zkcrypto/bls12_381/blob/main/src/pairings.rs

pub struct G1Affine {
  pub x: Fp,
  pub y: Fp,
  pub infinity: bool,
}

impl G1Affine {
    pub fn generator() -> Self {
        Self {
            x: BigNum {
                limbs: [
                    0x55e83ff97a1aeffb3af00adb22c6bb,
                    0x8c4f9774b905a14e3a3f171bac586c,
                    0xa73197d7942695638c4fa9ac0fc368,
                    0x17f1d3
                ]
            },
            y: BigNum {
                limbs: [
                    0x3cc744a2888ae40caa232946c5e7e1,
                    0xe095d5d00af600db18cb2c04b3edd0,
                    0x81e3aaa0f1a09e30ed741d8ae4fcf5,
                    0x8b3f4
                ]
            },
            infinity: false
        }
    }

    pub fn neg(self: Self) -> Self {
        let mut y = self.y.neg();
        if self.infinity {
            y = BigNum::one();
        }
        Self { x: self.x, y, infinity: self.infinity }
    }
}