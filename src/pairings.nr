use noir_bigcurve::curves::bls12_381::{BLS12_381Fq, BLS12_381_Params};
use crate::fp2::Fp2;
use crate::fp6::Fp6;
use crate::fp12::Fp12;
use bignum::BigNum;

// Ref: https://github.com/zkcrypto/bls12_381/blob/main/src/pairings.rs

pub struct G2Projective {
  x: Fp2,
  y: Fp2,
  z: Fp2
}

impl G2Projective {
    fn from_affine(p: G2Affine) -> Self {
        let mut z = Fp2::one();
        if (p.x.eq(Fp2::one()) & p.y.eq(Fp2::zero())) {
            z = Fp2::zero();
        }
        G2Projective { x: p.x, y: p.x, z }
    }
}

pub struct G2Affine {
  x: Fp2,
  y: Fp2,
  infinity: bool
}

pub type Fp = BLS12_381Fq;

pub struct G1Affine {
  x: Fp,
  y: Fp,
  infinity: bool
}

// https://github.com/zkcrypto/bls12_381/blob/main/src/pairings.rs#L709
fn doubling_step(r: G2Projective) -> (Fp2, Fp2, Fp2, G2Projective) {
    let mut tmp0 = r.x.square();
    let mut tmp1 = r.y.square();
    let mut tmp2 = tmp1.square();
    let mut tmp3 = ((tmp1.add(r.x)).square()).sub(tmp0).sub(tmp2);
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

fn addition_step(r: G2Projective, q: G2Affine) -> (Fp2, Fp2, Fp2, G2Projective) {
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

fn ell(f: Fp12, coeffs: (Fp2, Fp2, Fp2), p: G1Affine) -> Fp12 {
    let c0c0 = coeffs.0.c0.mul(p.y);
    let c0c1 = coeffs.0.c1.mul(p.y);

    let c1c0 = coeffs.1.c0.mul(p.x);
    let c1c1 = coeffs.1.c1.mul(p.x);

    f.mul_by_014(
        coeffs.2,
        Fp2 { c0: c0c0, c1: c0c1 },
        Fp2 { c0: c1c0, c1: c1c1 }
    )
}

fn miller_loop(p: G1Affine, q: G2Affine) -> Fp12 {
    let BLS_X = [
        true, true, false, true, false, false, true, false, false, false, false, false, false, false, false, true, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, true, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false
    ];

    let mut adder_cur = G2Projective::from_affine(q);
    let mut adder_base = q;
    let mut adder_p = p;
    let mut f = Fp12::one();
    let mut found_one = false;

    for i in 0..64 {
        let i_bool = BLS_X[i];
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

#[test]
fn test_doubling_step() {
    let generator = G2Projective {
        x: Fp2::new(
            BigNum {
            limbs: [
                0xac0326a805bbefd48056c8c121bdb8, 0x7ad4fa403b02b4510b647ae3d1770b, 0xb2f08f0a91260805272dc51051c6e4, 0x24aa2
            ]
        },
            BigNum {
            limbs: [
                0x4cf11213945d57e5ac7d055d042b7e, 0xd0d09920b61ab5da61bbdc7f504933, 0x6052719f607dacd3a088274f65596b, 0x13e02b
            ]
        }
        ),
        y: Fp2::new(
            BigNum {
            limbs: [
                0x3ac9cc3baca289e193548608b82801, 0x9baa8cbdd3a76d429a695160d12c92, 0x27727d6e118cc9cdc6da2e351aadfd, 0xce5d5
            ]
        },
            BigNum {
            limbs: [
                0x370d275cec1da1aaa9075ff05f79be, 0x287e85a763af267492ab572e99ab3f, 0xa02ea734cc32acd2b02bc28b99cb3e, 0x606c4
            ]
        }
        ),
        z: Fp2::one()
    };

    let a = Fp2::new(
        BigNum {
        limbs: [
            0x3f27323d5e8a27cc4e521822e0f559,
            0x23253f723bde4dd997044ed24e8e2a,
            0xb39075d1abe80b8f65256d2793537f,
            0x199642
        ]
    },
        BigNum {
        limbs: [
            0xdc349d73b07686aaa41d7fc17de6f8,
            0xa1fa169d8ebc99d24aad5cba66acfc,
            0x80ba9cd330cab34ac0af0a2e672cf8,
            0x181b12
        ]
    }
    );
    let b = Fp2::new(
        BigNum {
        limbs: [
            0xb8f2795d21163ce39975990c57da6d,
            0x1f7a219c6ac2d020fff505604caf82,
            0x2c22ea6825bb8178176c3a6ed7df18,
            0x291a
        ]
    },
        BigNum {
        limbs: [
            0xe781898edb7994dd769c7c2c6a8447,
            0x690e2cdf9728c6c7a4c90f9ad207a7,
            0x4f2c59babfdd1438e64394468e1eaf,
            0x8a764
        ]
    }
    );
    let c = Fp2::new(
        BigNum {
        limbs: [
            0x3cd6e1c4403fdca07f3855ab48b7d,
            0xe23363105a3956fc734cd655c20465,
            0x4989c84ac56fcf81ef6ef0c51408b3,
            0x116916
        ]
    },
        BigNum {
        limbs: [
            0x2136e9d5f30d406e4145ecc62c1f99,
            0xed1365d17357536a2d488e2a528552,
            0xac0820ccbd5704388bf7fc5530c4d9,
            0x1b696
        ]
    }
    );
    let doubled = G2Projective {
        x: Fp2::new(
            BigNum {
            limbs: [
                0xd403df5b1ac137922b915f1f38ec37, 0xf87ae23fe110c5aa48bc208fbd1a0e, 0xdd65245f8908a03288adcd396f489e, 0x177d52
            ]
        },
            BigNum {
            limbs: [
                0x8e26f1e5febc57b750d1c3dc4261d8, 0x5f169f0d4ddf46e23f15bb4f24304a, 0x8b9e689553d58194f79863fe02902c, 0xcf815
            ]
        }
        ),
        y: Fp2::new(
            BigNum {
            limbs: [
                0x730c92c30a69d0602a8804c2a7c65f, 0x9cdb85d9080b2c9b4766434c33d103, 0x215806e8a55fd2d9ec4af9d2d44859, 0x65ae9
            ]
        },
            BigNum {
            limbs: [
                0x361c973962955ba6405f0e51beb98b, 0xe41a04569ae33184419e66bbc0271c, 0x2d8a6d4b3a1cbd02c7bdc0e0aa304d, 0xe9c34
            ]
        }
        ),
        z: Fp2::new(
            BigNum {
            limbs: [
                0x75939877594513c326a90c11705002, 0x3755197ba74eda8534d2a2c1a25924, 0x4ee4fadc2319939b8db45c6a355bfb, 0x19cbaa
            ]
        },
            BigNum {
            limbs: [
                0x6e1a4eb9d83b4355520ebfe0bef37c, 0x50fd0b4ec75e4ce92556ae5d33567e, 0x405d4e69986559a56057851733967c, 0xc0d89
            ]
        }
        )
    };

    // Result are 3 Fp2's and the adjusted projective point 
    // (in Rust this is changed in place, but this is not supported in Noir so it is an additional return value)
    let res = doubling_step(generator);

    // Check changed point
    assert(res.3.x.eq(doubled.x));
    assert(res.3.y.eq(doubled.y));
    assert(res.3.z.eq(doubled.z));

    // Check Fp2's 
    assert(res.0.eq(a));
    assert(res.1.eq(b));
    assert(res.2.eq(c));
}

#[test]
fn test_addition_step() {
    /*
  Values from zkcrypto impl
  #[test]
  fn test_addition_step() {
    let a = Scalar::from_raw([1, 2, 3, 4]).invert().unwrap().square();
    let test_value = G2Affine::generator()*a;
    let mut test_val_proj = G2Projective::from(test_value);
    println!("{}\n", test_val_proj);

    let t = test_value.into();
    println!("{}\n", t);
    let dbl = addition_step(&mut test_val_proj, &t);
    println!("{:?}\n", dbl);
    println!("{}\n", test_val_proj);
  }
  */
    let point1 = G2Projective {
        x: Fp2::new(
            BigNum {
            limbs: [
                0x335fee64bd06db6573fdcb9ee98319, 0x593fac9c1b9b2480ae1c7f4e45f9bb, 0x6da041aadfe26730abaac871996971, 0x7c7d4
            ]
        },
            BigNum {
            limbs: [
                0xf3f398d2eec5613d7e8fb8eea72f72, 0xdbaafceeb0accc8abf24d95553a45c, 0x15344fea28e0fd9ab57be98b054672, 0x2531e
            ]
        }
        ),
        y: Fp2::new(
            BigNum {
            limbs: [
                0x69c8bd9c9d71502f1318fe1be21cc0, 0xc079d0edfc279c0052d38038f28fc5, 0xadb0d4365ff48c2cd32891ac2615a8, 0x75e8e
            ]
        },
            BigNum {
            limbs: [
                0xc3230bf99a869a5f4a6627fe3851b9, 0xa833d8fb5f0e5559cd8e03356c9e4d, 0x889e30f914fe674802ed8803f0f0d7, 0x590bf
            ]
        }
        ),
        z: Fp2::new(
            BigNum {
            limbs: [
                0x91320d77ab2c01f0cf4b79f8644bc8, 0x375ffeade6cac1c62cf9ab9052de7, 0x706f3ced1794bc6fd493745e2340a5, 0xda79
            ]
        },
            BigNum {
            limbs: [
                0xc39955cb1638df156cbab123f66dbc, 0xac805f7513b6a5a8e1e1258d57f404, 0x8940c37a963a03b135428745f2340c, 0x18c088
            ]
        }
        )
    };

    let affine = G2Affine {
        x: Fp2::new(
            BigNum {
            limbs: [
                0xb140828914d4659e103008b0f99ce2, 0xab1484bd130dc52e9a993b56bfa96e, 0xad2ce8990405d4fafad4edbce84d24, 0x71593
            ]
        },
            BigNum {
            limbs: [
                0x7866ec1320541b8e914140cada9eac, 0xecc253bae630f295cb97d79685e9ee, 0x2a956ba1f0c00f7ac188eeda032acb, 0x6322e
            ]
        }
        ),
        y: Fp2::new(
            BigNum {
            limbs: [
                0x688f932215126dacf1612804d3dac3, 0xd4fb713edca1dbb4cbfa81de69b297, 0x423f261bb8f9441937551b8c3cc029, 0xcacc
            ]
        },
            BigNum {
            limbs: [
                0x3840064ad1355cd762cc26266d678f, 0x5ff57854ba0a2e170720dc73b35053, 0xbff73fc2bd237da821284216f8dc9f, 0xaab8c
            ]
        }
        ),
        infinity: false
    };

    let a = Fp2::new(
        BigNum {
        limbs: [
            0x1bf08d8f6c261980aa6ec9377eeb5a,
            0x560c420d56d1257cad1c8daa6ffd6a,
            0x15ac001d9750dc97595e05a3ab81a1,
            0x19e61c
        ]
    },
        BigNum {
        limbs: [
            0x24330a557692fa7b0c56585c2c0696,
            0x5eb50c857006e9d130313da2ccbbc6,
            0xfaac085ce0a7842b2aa3044ca02fe7,
            0x10048e
        ]
    }
    );
    let b = Fp2::new(
        BigNum {
        limbs: [
            0x589dad22b7c810c06691568df9e5c5,
            0xb21119bda48ce10170e15e1d6b435c,
            0xff859681f3e36d0ca9803e70a898f7,
            0x8e196
        ]
    },
        BigNum {
        limbs: [
            0xe9c4e582ffd6d9076c4caeefea9722,
            0x2895927e77470516516482f3284d77,
            0x85492b6009b98ac5a4990e3225a290,
            0x80188
        ]
    }
    );
    let c = Fp2::new(
        BigNum {
        limbs: [
            0x6be25edc8c984a0727703e81c7369e,
            0x556f550564310a5e0521371de6e410,
            0xeed02b82a7cfb8bbe3547aed326bb7,
            0x8e06c
        ]
    },
        BigNum {
        limbs: [
            0xe8e4c2f187c9ba393c2f72314a8bba,
            0x573eb0a205bc6a8eb23577b3f87754,
            0x14d1f9e61982c4efa5fdc5d9ae9afe,
            0x153622
        ]
    }
    );
    let addition_res_expected = G2Projective {
        x: Fp2::new(
            BigNum {
            limbs: [
                0xe60e8ea5133cd50e9bbd821aebf8b0, 0x264be53d3aba569e1414ee1cbdeed6, 0xf27b9d49940bc78483c46b66300132, 0x18fe28
            ]
        },
            BigNum {
            limbs: [
                0xd85e3b56371e75a3e1ab6f4de4684c, 0x645f74f93b077a6d50bb261c4fb3f8, 0x7a0d54da0e220533084ab75e195f31, 0x78ff4
            ]
        }
        ),
        y: Fp2::new(
            BigNum {
            limbs: [
                0xdc24946be6b95b618be1503266315e, 0x477f1fdb9da707a951559beccf63d3, 0xdfb35908d8a5d1f5de994c0ab7dd28, 0x123c72
            ]
        },
            BigNum {
            limbs: [
                0x61b1cac14443e5fc724d1f17117692, 0x3df3f01f65da37a4816ee5269ad530, 0x8744f8a6346d7242623e875b47e064, 0x17ab85
            ]
        }
        ),
        z: Fp2::new(
            BigNum {
            limbs: [
                0xdf846c7b6130cc05537649bbf75ad, 0xab062106ab6892be568e46d537feb5, 0xad6000ecba86e4bacaf02d1d5c0d0, 0xcf30e
            ]
        },
            BigNum {
            limbs: [
                0x1219852abb497d3d862b2c2e16034b, 0xaf5a8642b80374e898189ed1665de3, 0x7d56042e7053c215955182265017f3, 0x80247
            ]
        }
        )
    };

    // Result are 3 Fp2's and the adjusted projective point 
    // (in Rust this is changed in place, but this is not supported in Noir so it is an additional return value)
    let res: (Fp2, Fp2, Fp2, G2Projective) = addition_step(point1, affine);

    // Check changed point
    assert(res.3.x.eq(addition_res_expected.x));
    assert(res.3.y.eq(addition_res_expected.y));
    assert(res.3.z.eq(addition_res_expected.z));

    // Check Fp2's 
    assert(res.0.eq(a));
    assert(res.1.eq(b));
    assert(res.2.eq(c));
}
