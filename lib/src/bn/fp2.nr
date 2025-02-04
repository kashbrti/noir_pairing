use bignum::BigNum;
use bignum::fields::bn254Fq::BN254_Fq_Params;

/// Field Element in BN254
///
/// Uses `BigNum` with 3 limbs, 254 bit modulus, and BN254 parameters defined in `noir-bignum`.
pub type Fp = BigNum<3, 254, BN254_Fq_Params>;

// c0 + c1*u
// F_p / u^2 + 1
/// Polynomial Extension Field of Degree `2`
///
/// Representable as $F_p[u] / u^2 + 1$, we implement `c0 + c1 * u` such that `c0` and `c1` are
/// elements of `Fp`.
pub struct Fp2 {
    /// Constant Coefficient
    pub(crate) c0: Fp,
    /// Exponent `1` Coefficient
    pub(crate) c1: Fp,
}

pub trait ExtensionFieldTrait {
    fn zero() -> Self;
    fn one() -> Self;
    fn neg(self: Self) -> Self;
    fn double(self: Self) -> Self;
    fn square(self: Self) -> Self;
    fn invert(self: Self) -> Self;
    fn frobenius_map(self: Self) -> Self;
}


impl ExtensionFieldTrait for Fp2 {
    /// Constructs the additive unit (zero) element.
    fn zero() -> Fp2 {
        Fp2 { c0: BigNum::new(), c1: BigNum::new() }
    }

    /// Constructs the multiplicative unit (one) element.
    fn one() -> Fp2 {
        Fp2 { c0: BigNum::one(), c1: BigNum::new() }
    }

    /// Returns the element's additive inverse (negation).
    fn neg(self: Self) -> Self {
        Self { c0: self.c0.neg(), c1: self.c1.neg() }
    }

    

    // Returns the double of the element
    fn double(self: Self) -> Self {
        let c0 = self.c0; 
        let c1 = self.c1;
        let res_c0 = c0 + c0;
        let res_c1 = c1 + c1;
        Self { c0: res_c0, c1: res_c1 }
    }

    /// Returns the square of the element.
    fn square(self: Self) -> Self {
        unsafe {
            // a = c0 + c1
            // b = c0 - c1
            // c = c0 + c0
            let a = self.c0.__add(self.c1);
            let b = self.c0.__sub(self.c1);
            let c = self.c0.__add(self.c0);

            // res_c0 = (c0 + c1) * (c0 - c1)
            let res_c0 = a.__mul(b);
            // res_c1 = (c0 + c0) * c1
            let res_c1 = c.__mul(self.c1);

            // Constrain res_c0
            let lhs_terms = [[self.c0, self.c1]];
            let lhs_flags = [[false, false]];
            let rhs_terms = [[self.c0, self.c1]];
            let rhs_flags = [[false, true]];
            let res = [res_c0];
            let res_flag = [true];
            BigNum::evaluate_quadratic_expression(
                lhs_terms,
                lhs_flags,
                rhs_terms,
                rhs_flags,
                res,
                res_flag,
            );

            // Constrain res_c1
            let lhs_terms1 = [[self.c0, self.c0]];
            let lhs_flags1 = [[false, false]];
            let rhs_terms1 = [[self.c1, BigNum::new()]];
            let rhs_flags1 = [[false, false]];
            let res1 = [res_c1];
            let res_flag1 = [true];
            BigNum::evaluate_quadratic_expression(
                lhs_terms1,
                lhs_flags1,
                rhs_terms1,
                rhs_flags1,
                res1,
                res_flag1,
            );

            Self { c0: res_c0, c1: res_c1 }
        }
    }

    /// Returns the elements multiplicative inverse.
    ///
    /// ## Panics
    ///
    /// - if there is no multiplicative inverse for the element.
    fn invert(self: Self) -> Self {
        unsafe {
            // (self.c0 * self.c0 + self.c1 * self.c1).invert()
            let c0_square = self.c0.__mul(self.c0);
            let c1_square = self.c1.__mul(self.c1);
            let temp = (c0_square.__add(c1_square));

            // Constrain temp = self.c0 * self.c0 + self.c1 * self.c1
            let lhs_terms = [[self.c0], [self.c1]];
            let lhs_flags = [[false], [false]];
            let rhs_terms = [[self.c0], [self.c1]];
            let rhs_flags = [[false], [false]];
            let res = [temp];
            let res_flag = [true];
            BigNum::evaluate_quadratic_expression(
                lhs_terms,
                lhs_flags,
                rhs_terms,
                rhs_flags,
                res,
                res_flag,
            );

            // __invmod asserts that multiplicative inverse exists according to parameters
            // in this case that should pass
            let inv = temp.__invmod();
            // Constrain inverse
            // TODO can this be done more efficiently?
            let check_identity = inv * temp;
            assert(check_identity == BigNum::one());

            // Computing c0, c1 constrained.
            // TODO how can this be improved?
            // c0 = self.c0 * inv
            // c1 = self.c1 * (-inv)
            let res_c0 = self.c0 * inv;
            let res_c1 = self.c1 * inv.neg();

            Self { c0: res_c0, c1: res_c1 }
        }
    }



    /// Returns the result of the frobenius endomorphism over `Fp2`
    ///
    /// Frobenius exponentiates the element to the prime field modulus `p` such that raising the
    /// constant term, `c0` to the `p` returns `c0`.
    ///
    /// `(c0 ** p) % p == c0`
    ///
    /// However, raising the exponent `1` term `c1 * u` to the `p` returns `-c1 * u`.
    ///
    /// `(c1 ** p) % p == c1`
    /// `(u ** p) % p == -u`
    ///
    /// Thus we return the result of `Fp2::conjugate`.
    fn frobenius_map(self: Self) -> Self {
        // (c0 + c1*u)^p = c0^p + (c1*u)^p
        // and c0^p=c0, c1^p=c1, and u^p = -u, so the result is:
        // c0 - c1
        self.conjugate()
    }
}

impl Fp2 {
    /// Constructs an element of `Fp2` from two field arrays.
    ///
    /// > TODO: Should this be parameterized by `Fp` instances instead?
    pub fn new(c0: [Field; 3], c1: [Field; 3]) -> Self {
        Fp2 { c0: BigNum { limbs: c0 }, c1: BigNum { limbs: c1 } }
    }

    /// Returns true if both coefficients are zero.
    pub fn is_zero(self: Self) -> bool {
        (self.c0 == BigNum::new()) & (self.c1 == BigNum::new())
    }


    /// Multiplies the element by the quadratic non-residue `u` in `Fp2`.
    ///
    /// This is used for optimized multiplication over `Fp2`.
    ///
    /// `c0 + c1*u` simplifies to `(c0 - c1) + (c0 + c1)*u`
    pub fn mul_by_non_residue(self: Self) -> Self {
        // Multiply c0 + c1*u by u + 9:
        let mult_by_nine = self.double().double().double().add(self); 
        let mult_by_u = Self { c0: self.c1.neg(), c1: self.c0 };
        mult_by_nine + mult_by_u
    }

    /// Multiplies the element by a non-residue of power 1
    ///
    /// Reference Document: https://hackmd.io/@Wimet/ry7z1Xj-2#Frobenius-Operator
    fn mul_by_non_residue_power_1(self: Self) -> Self {
        let y = Fp2::new(
            [426377461259128476197171212816344176, 950761433168210257943908058976182364, 4740],
            [631238624071726213905819414132122284, 783787493034767597529245196718149236, 9321],
        );

        self.mul(y)
    }

    /// Multiplies the element by a non-residue of power 2
    ///
    /// Reference Document: https://hackmd.io/@Wimet/ry7z1Xj-2#Frobenius-Operator
    fn mul_by_non_residue_power_2(self: Self) -> Self {
        let y = Fp2::new(
            [730877957994261685262562863694697789, 371742299709972197061924180894761911, 12211],
            [904313815760113414035932210111957922, 1190666331444014455945460071356929074, 5833],
        );

        self.mul(y)
    }

    /// Multiplies the element by a non-residue of power 3
    ///
    /// Reference Document: https://hackmd.io/@Wimet/ry7z1Xj-2#Frobenius-Operator
    fn mul_by_non_residue_power_3(self: Self) -> Self {
        let y = Fp2::new(
            [887252555323038411621405899505013594, 1261835301057603723149303145102424539, 1596],
            [667165268406215768209883118386262243, 315364551883617340788418358721588001, 1984],
        );
        self.mul(y)
    }

    /// Multiplies the element by a non-residue of power 4
    ///
    /// Reference Document: https://hackmd.io/@Wimet/ry7z1Xj-2#Frobenius-Operator
    fn mul_by_non_residue_power_4(self: Self) -> Self {
        let y = Fp2::new(
            [266875291493575529387049817426798434, 412105995003785042396868112665183443, 1461],
            [694853336217345651450901527565926694, 492536132578680416344610139344245302, 11284],
        );
        self.mul(y)
    }

    /// Multiplies the element by a non-residue of power 5
    ///
    /// Reference Document: https://hackmd.io/@Wimet/ry7z1Xj-2#Frobenius-Operator
    fn mul_by_non_residue_power_5(self: Self) -> Self {
        let y = Fp2::new(
            [479103104124920441723708103336618287, 1006804826950535178307367240931996738, 387],
            [668743480999689269792587498111380555, 1260642313736059160892374157579222448, 4780],
        );
        self.mul(y)
    }

    /// Returns the conjugate of the element.
    ///
    /// Functionally, this inverts the exponent `1` coefficient, not changing the constant.
    ///
    /// See the `Fp2::frobenius_map` for more information.
    pub fn conjugate(self: Self) -> Self {
        Fp2 { c0: self.c0, c1: self.c1.neg() }
    }


    /// Scales the element by an element of Fp.
    pub fn scale(self, scalar: Fp) -> Self {
        Self { c0: self.c0 * scalar, c1: self.c1 * scalar }
    }
}

impl std::ops::Add for Fp2 {
    fn add(self: Self, other: Self) -> Self {
        Self { c0: self.c0 + other.c0, c1: self.c1 + other.c1 }
    }
}

impl std::ops::Sub for Fp2 {
    fn sub(self: Self, other: Self) -> Self {
        Self { c0: self.c0 - other.c0, c1: self.c1 - other.c1 }
    }
}

impl std::ops::Mul for Fp2 {
    fn mul(self: Self, other: Self) -> Self {
        // TODO: don't wrap the whole thing in unsafe.
        unsafe {
            // Calculate unconstrained
            // self.c0 * other.c0
            let temp = self.c0.__mul(other.c0);
            // self.c1 * other.c1
            let temp1 = self.c1.__mul(other.c1);
            // c0 = self.c0 * other.c0 - self.c1 * other.c1
            let c0 = temp.__sub(temp1);

            // (self.c0 + self.c1)
            let temp2 = self.c0.__add(self.c1);
            // (other.c0 + other.c1)
            let temp3 = other.c0.__add(other.c1);
            let temp4 = temp2.__mul(temp3);
            // self.c0 * other.c0
            let temp5 = self.c0.__mul(other.c0);
            // self.c1 * other.c1
            let temp6 = self.c1.__mul(other.c1);
            // c1 = (self.c0 + self.c1) * (other.c0 + other.c1) - self.c0 * other.c0 - self.c1 * other.c1
            let c1 = temp4.__sub(temp5).__sub(temp6);

            // Constrain c0
            // c0 = self.c0 * other.c0 - self.c1 * other.c1
            let lhs_terms = [[self.c0], [self.c1]];
            let lhs_flags = [[false], [true]];
            let rhs_terms = [[other.c0], [other.c1]];
            let rhs_flags = [[false], [false]];
            let res = [c0];
            let res_flag = [true];
            BigNum::evaluate_quadratic_expression(
                lhs_terms,
                lhs_flags,
                rhs_terms,
                rhs_flags,
                res,
                res_flag,
            );

            // Constrain c1
            // c1 = (self.c0 + self.c1) * (other.c0 + other.c1) - self.c0 * other.c0 - self.c1 * other.c1
            // (adding zeroes for compilation)
            let lhs_terms1 =
                [[self.c0, self.c1], [self.c0, BigNum::new()], [self.c1, BigNum::new()]];
            let lhs_flags1 = [[false, false], [true, false], [true, false]];
            let rhs_terms1 =
                [[other.c0, other.c1], [other.c0, BigNum::new()], [other.c1, BigNum::new()]];
            let rhs_flags1 = [[false, false], [false, false], [false, false]];
            let res1 = [c1];
            let res_flag1 = [true];
            BigNum::evaluate_quadratic_expression(
                lhs_terms1,
                lhs_flags1,
                rhs_terms1,
                rhs_flags1,
                res1,
                res_flag1,
            );

            // c0 = self.c0 * other.c0 - self.c1 * other.c1
            // c1 = (self.c0 + self.c1) * (other.c0 + other.c1) - self.c0 * other.c0 - self.c1 * other.c1
            Self { c0, c1 }
        }
    }
}

impl std::cmp::Eq for Fp2 {
    fn eq(self: Self, other: Self) -> bool {
        (self.c0 == other.c0) & (self.c1 == other.c1)
    }
}

// #[test]
// fn test_mul_1() {
//     let x = [
//         0x55e83ff97a1aeffb3af00adb22c6bb,
//         0x8c4f9774b905a14e3a3f171bac586c,
//         0xa73197d7942695638c4fa9ac0fc368,
//         0x17f1d3
//     ];
//     let y = [0, 0, 0, 0];

//     let zero = Fp2::zero();
//     let test_val = Fp2::new(x, y);

//     // test multiplication by 1
//     let mult_by_1 = test_val.mul(Fp2::one());

//     assert(mult_by_1.c0 == test_val.c0);
//     assert(mult_by_1.c1 == test_val.c1);

//     // test multiplication by 0
//     let mult_by_0 = test_val.mul(zero);
//     assert(mult_by_0.c0 == BigNum::new());
//     assert(mult_by_0.c1 == BigNum::new());
// }

#[test]
fn test_one_mul() {
    let one = Fp2::one();
    let test_val = Fp2::new([1, 0, 0], [0, 0, 0]);
    let res = one.mul(test_val);
    assert(res == test_val);
}

#[test]
fn test_zero_mul() {
    let zero = Fp2::zero();
    let test_val = Fp2::new([1, 0, 0], [0, 0, 0]);
    let res = zero.mul(test_val);
    assert(res == zero);
}

#[test]
fn test_add_0() {
    let a = Fp2::zero();
    let b = Fp2::new([1, 0, 0], [0, 0, 0]);
    let res = a.add(b);
    assert(res == b);
}

#[test]
fn test_add_1() {
    let a = Fp2::new(
        [1313602686734452469215463949774943363, 952032783256258840911419402804250364, 2982],
        [1106169269385670298369265154383545602, 492176543493153977447485047093792002, 8430],
    );
    let b: _ = Fp2::new(
        [1058505771725674244616784819913523546, 1157098086699933031536880181494638523, 1958],
        [273210913073766385990275248560883300, 864307921850106992811755483719817132, 10077],
    );
    let c = Fp2::new(
        [1042880462675210840928441709408122333, 779902874171275999544492524018544312, 4941],
        [50152186674520811455733342664084326, 27256469558345097355433470533264559, 18508],
    );
    let res = a.add(b);
    assert(res.c0 == c.c0);
    assert(res.c1 == c.c1);
}

#[test]
fn test_eq() {
    let a = Fp2::one();
    let b = Fp2::new([1, 0, 0], [0, 0, 0]);
    assert(a == b);
}

#[test]
fn test_sub_0() {
    let a = Fp2::zero();
    let b = Fp2::new([1, 0, 0], [0, 0, 0]);
    let res = b - a;
    assert(res == b);
}

#[test]
fn test_neg() {
    let a = Fp2::new(
        [1313602686734452469215463949774943363, 952032783256258840911419402804250364, 2982],
        [1106169269385670298369265154383545602, 492176543493153977447485047093792002, 8430],
    );
    let b: _ = Fp2::new(
        [1058505771725674244616784819913523546, 1157098086699933031536880181494638523, 1958],
        [273210913073766385990275248560883300, 864307921850106992811755483719817132, 10077],
    );
    let bneg = b.neg();
    assert(a + bneg == a - b);
    assert(b + bneg == Fp2::zero());
}

#[test]
fn test_mul_by_nonresidue() {
    let a = Fp2::one();
    let b = a.mul_by_non_residue();
    let c: Fp2 = Fp2::new([1, 0, 0], [1, 0, 0]);
    let d = c.mul_by_non_residue();
    assert(b == Fp2::new([9, 0, 0], [1, 0, 0]));
    println(d);
    assert(d == Fp2::new([8, 0, 0], [10, 0, 0]));
}

// #[test]
// fn test_mul_by_non_residue_power() {
//     let a = Fp2::one();
//     let b = a.mul_by_non_residue();
//     assert(b == Fp2::new([9, 0, 0], [1, 0, 0]));
// }

// TODO: write tests for everything elses

// #[test]
// fn test_subtraction() {
//     let a = Fp2::new(
//         [
//         0xae4a11560fd2b96a09231943c6dab,
//         0xd18b53078ed5d3f5cadea0c65f2734,
//         0xff3ed67578ddeae5439da2b02eab31,
//         0x8ccec
//     ],
//         [
//         0x38f87f1c9644538b4bf7ec8dee1be9,
//         0x5fb2ff8fddcf6fbb096adcc2b04bfc,
//         0x56170fab6f325e416abe497d9e055d,
//         0xbd355
//     ]
//     );
//     let b = Fp2::new(
//         [
//         0x4579a52cfa8b4f01e3fc00c5962be5,
//         0x45b3857ffd05157321b22d5a04c013,
//         0xb7e3b626ce947b122db3425e68ae97,
//         0xfc9d8
//     ],
//         [
//         0x8fd0966dcadba0b5f86e4c5ffa1da2,
//         0x18864f9cef0cf179640918b4d31a73,
//         0xf2af61549b350ae4025242be1008a6,
//         0x325f3
//     ]
//     );
//     let c = Fp2::new(
//         [
//         0x716afa99ba71dc4ebb9630cea5ec71,
//         0xd75cc10ca49025b37bcd6a1d508b3f,
//         0x3194a03544948b7acc2dabfe9d6111,
//         0x130426
//     ],
//         [
//         0xa927e8aecb68b2d55389a02df3fe47,
//         0x472caff2eec27e41a561c40ddd3188,
//         0x6367ae56d3fd535d686c06bf8dfcb7,
//         0x8ad61
//     ]
//     );
//     let res = a.sub(b);
//     assert(res.c0 == c.c0);
//     assert(res.c1 == c.c1);
// }

// #[test]
// fn test_multiplication() {
//     // zkcrypto test
//     let a = Fp2::new(
//         [
//         0xae4a11560fd2b96a09231943c6dab,
//         0xd18b53078ed5d3f5cadea0c65f2734,
//         0xff3ed67578ddeae5439da2b02eab31,
//         0x8ccec
//     ],
//         [
//         0x38f87f1c9644538b4bf7ec8dee1be9,
//         0x5fb2ff8fddcf6fbb096adcc2b04bfc,
//         0x56170fab6f325e416abe497d9e055d,
//         0xbd355
//     ]
//     );
//     let b = Fp2::new(
//         [
//         0x4579a52cfa8b4f01e3fc00c5962be5,
//         0x45b3857ffd05157321b22d5a04c013,
//         0xb7e3b626ce947b122db3425e68ae97,
//         0xfc9d8
//     ],
//         [
//         0x8fd0966dcadba0b5f86e4c5ffa1da2,
//         0x18864f9cef0cf179640918b4d31a73,
//         0xf2af61549b350ae4025242be1008a6,
//         0x325f3
//     ]
//     );
//     let c = Fp2::new(
//         [
//         0x77b92dc10d7dc4df06bbf3f0fe3a47,
//         0x68be63dff75f16e1ac619556ce0f96,
//         0xeade1ad5acb2d25bcc90c425869481,
//         0x80690
//     ],
//         [
//         0xca7635cc60dd923e8c7db21115127e,
//         0x66780a57f3ae3477f44f30c19df036,
//         0xb168c449a974a86ef04b08f869126,
//         0x6e753
//     ]
//     );

//     let res = a.mul(b);
//     assert(res.c0 == c.c0);
//     assert(res.c1 == c.c1);
// }

// #[test]
// fn test_square() {
//     let a = Fp2::new(
//         [
//         0xae4a11560fd2b96a09231943c6dab,
//         0xd18b53078ed5d3f5cadea0c65f2734,
//         0xff3ed67578ddeae5439da2b02eab31,
//         0x8ccec
//     ],
//         [
//         0x38f87f1c9644538b4bf7ec8dee1be9,
//         0x5fb2ff8fddcf6fbb096adcc2b04bfc,
//         0x56170fab6f325e416abe497d9e055d,
//         0xbd355
//     ]
//     );
//     let c = Fp2::new(
//         [
//         0x4579a52cfa8b4f01e3fc00c5962be5,
//         0x45b3857ffd05157321b22d5a04c013,
//         0xb7e3b626ce947b122db3425e68ae97,
//         0xfc9d8
//     ],
//         [
//         0x8fd0966dcadba0b5f86e4c5ffa1da2,
//         0x18864f9cef0cf179640918b4d31a73,
//         0xf2af61549b350ae4025242be1008a6,
//         0x325f3
//     ]
//     );
//     let res = a.square();
//     assert(res.c0 == c.c0);
//     assert(res.c1 == c.c1);
// }

// #[test]
// fn test_inversion() {
//     let a = Fp2::new(
//         [
//         0xb77accf6a4e5e337a85d2a62baf070,
//         0x5e5a337d7db6344484087a18586c7c,
//         0xb3f92a8e29eabd22358787933cd419,
//         0xa4aa
//     ],
//         [
//         0xb33c802bede6916fa9f3d1d08a643b,
//         0x4537aef2966b539d7364999eb35707,
//         0xde1e413ba79349fbe6cc88298b9f18,
//         0x1776be
//     ]
//     );
//     let c = Fp2::new(
//         [
//         0x71825b46349b307298025eea77770c,
//         0x225cbdc399bb960534036f625a3592,
//         0xb32ff3d468eb1aba35968a2bbdb9d0,
//         0x66e47
//     ],
//         [
//         0x3fc9a610ece44f3c738a358aa60150,
//         0x602cf426b135e0e0d837a8cf5d9ec7,
//         0x4aee9cb5468ad40ee5dccbf0264407,
//         0xdc722
//     ]
//     );

//     let res = a.invert();
//     assert(res.c0 == c.c0);
//     assert(res.c1 == c.c1);

//     let identity_check = res.mul(a);
//     assert(identity_check.eq(Fp2::one()));
// }

