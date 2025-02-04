use crate::bn::fp2::ExtensionFieldTrait;
use crate::bn::fp2::Fp2;
use crate::bn::fp2::Fp;
use crate::bn::fp6::Fp6;
use crate::bn::pairing::Ell;
use crate::bn::mod_p12m1::ModP12M1Params;
use bignum::BigNum;
use bignum::BigNumTrait;
use bignum::utils::u60_representation::U60Repr;




/// Polynomial Extension Field Tower of Degree `12`.
///
/// Representable as $F_{p^{12}}[w] = F_{p^6}[w] / (w^2 - v)$, we implement `c0 + c1 * w` such that
/// `c0` and `c1` are elements of `Fp6`.
pub struct Fp12 {
    /// Constant Coefficient
    pub(crate) c0: Fp6,
    /// Exponent `1` Coefficient
    pub(crate) c1: Fp6,
}



impl ExtensionFieldTrait for Fp12 {
    /// Constructs the additive unit (zero) element.
    fn zero() -> Self {
        Fp12 { c0: Fp6::zero(), c1: Fp6::zero() }
    }

    /// Constructs the multiplicative unit (one) element.
    fn one() -> Self {
        Fp12 { c0: Fp6::one(), c1: Fp6::zero() }
    }

    /// Returns the element's additive inverse (negation).
    fn neg(self: Self) -> Self {
        Self { c0: self.c0.neg(), c1: self.c1.neg() }
    }

    /// Returns the square of the element.
    fn square(self: Fp12) -> Fp12 {
        // ab = self.c0 * self.c1
        let ab = self.c0 * self.c1;

        // c0_plus_c1 = self.c0 + self.c1
        let c0_plus_c1 = self.c0 + self.c1;

        // temp = self.c1 * nonresidue
        let temp = self.c1.mul_by_nonresidue();

        // temp_with_c0 = temp + self.c0
        // temp_with_c0 = (self.c1 * nonresidue) + self.c0
        let temp_with_c0 = temp + self.c0;

        // step 1: c0 = temp_with_c0 * c0_plus_c1 - ab
        // c0 = (self.c1 * nonresidue + self.c0) * (self.c0 + self.c1) - ab
        let mut c0 = (temp_with_c0 * c0_plus_c1) - ab;
        // step 2: c0 = ((self.c1 * nonresidue + self.c0) * (self.c0 + self.c1) - ab) - (ab * nonresidue)
        c0 = c0 - ab.mul_by_nonresidue();

        // c1 = 2 * ab
        let c1 = ab + ab;

        // Return the final Fp12 structure
        Fp12 { c0, c1 }
    }

    /// Returns the result of the frobenius endomorphism over `Fp12`.
    ///
    /// For optimization purposes, `Fp12` being composed of `Fp6` elements compmosed of `Fp2`
    /// elements enables a faster frobenius endomorphism over `Fp2` (see `Fp2::frobenius_map`) and
    /// `Fp6`, finally we perform multiplication over `c1` by the constant defined below.
    fn frobenius_map(self: Self) -> Self {
        // Fp12 = Fp2[x]/X^6 - u
        // decompose the point into its Fp2 coefficients
        // so if we have f = g + h w where h , g are in Fp6
        // and h = h_0 + h_1 * v + h_2 * v^2
        // and g = g_0 + g_1 * v + g_2 * v^2
        // then f = g_0 + h_0 w + g_1 w^2 + h_1 w^3 + g_2 w^4 + h_2 w^5
        // so to present Fp12 elements as t_0 + t_1 w + t_2 w^2 + t_3 w^3 + t_4 w^4 + t_5 w^5
        // we can use the following map
        // t_i s are elements in Fp2
        // we have to apply the frobenius on each coefficient in FP2 which is a simple conjugation map
        let mut t0 = self.c0.c0.conjugate();
        let mut t1 = self.c1.c0.conjugate();
        let mut t2 = self.c0.c1.conjugate();
        let mut t3 = self.c1.c1.conjugate();
        let mut t4 = self.c0.c2.conjugate();
        let mut t5 = self.c1.c2.conjugate();
        // now we multiply by the precomputed powers of the nonresidue
        t1 = t1.mul_by_non_residue_power_1();
        t2 = t2.mul_by_non_residue_power_2();
        t3 = t3.mul_by_non_residue_power_3();
        t4 = t4.mul_by_non_residue_power_4();
        t5 = t5.mul_by_non_residue_power_5();

        Fp12 { c0: Fp6::new(t0, t2, t4), c1: Fp6::new(t1, t3, t5) }
    }

    fn double(self) -> Self {
        self + self
    }

    /// Returns the elements multiplicative inverse.
    ///
    /// ## Panics
    ///
    /// - if there is no multiplicative inverse for the element.
    pub fn invert(self: Self) -> Self {
        let c02 = self.c0.mul(self.c0);
        let c12 = self.c1.mul(self.c1);
        let c12 = c12.mul_by_nonresidue();
        let c02 = c02.sub(c12);
        let c0 = c02.invert();
        let c1 = c02.invert();
        let c0 = c0.mul(self.c0);
        let c1 = c1.mul(self.c1);
        let c1 = c1.neg();

        Self { c0, c1 }
    }



}

impl Fp12 {
    /// Returns the conjugate of the element.
    ///
    /// Functionally, this inverts the exponent `1` coefficient, not changing the constant.
    pub fn from_coeff(limbs: [[Field; 3]; 12]) -> Self {
        let c00 = Fp2 {c0: Fp {limbs: limbs[0]}, c1: Fp {limbs: limbs[1]}}; 
        let c01 = Fp2 {c0: Fp {limbs: limbs[2]}, c1: Fp {limbs: limbs[3]}}; 
        let c02 = Fp2 {c0: Fp {limbs: limbs[4]}, c1: Fp {limbs: limbs[5]}}; 
        let c10 = Fp2 {c0: Fp {limbs: limbs[6]}, c1: Fp {limbs: limbs[7]}}; 
        let c11 = Fp2 {c0: Fp {limbs: limbs[8]}, c1: Fp {limbs: limbs[9]}}; 
        let c12 = Fp2 {c0: Fp {limbs: limbs[10]}, c1: Fp {limbs: limbs[11]}}; 
        let c0 = Fp6 {c0: c00, c1: c01, c2: c02}; 
        let c1 = Fp6 {c0: c10, c1: c11, c2: c12}; 
        Fp12 {c0, c1}
    }

    pub fn conjugate(self: Self) -> Self {
        Self { c0: self.c0, c1: self.c1.neg() }
    }

    /// Returns the result of the frobenius endomorphism two times over `Fp12`
    ///
    /// > TODO: switch the second degree frobenius map to use the precomputation table
    pub fn second_degree_frobenius_map(self: Self) -> Self {
        self.frobenius_map().frobenius_map()
    }

    /// Returns the result of the frobenius endomorphism three times over `Fp12`
    ///
    /// > TODO: switch the third degree frobenius map to use the precomputation table
    pub fn third_degree_frobenius_map(self: Self) -> Self {
        self.frobenius_map().frobenius_map().frobenius_map()
    }

    /// Exponentiate an element of `Fp12` to a constant `lambda`.
    ///
    /// $q = 36x^4 + 36x^3 + 24x^2 + 6x + 1$
    ///
    /// $\lambda = 6x + 2 + q - q^2 + q^3$
    ///
    /// > Note: Using a loop over the bits of $\lambda$ does not automatically unroll itself, so we
    /// > manually unroll the loop, performing `exp*=a; a*=a` on `1` bits and `a*=a` on `0` bits.
    pub fn exp_lambda(self) -> Self {
        let mut a = self;
        let mut exp = Fp12::one();

        exp *= a;
        a *= a;
        a *= a;
        a *= a;
        a *= a;
        exp *= a;
        a *= a;
        exp *= a;
        a *= a;
        exp *= a;
        a *= a;
        exp *= a;
        a *= a;
        exp *= a;
        a *= a;
        a *= a;
        a *= a;
        exp *= a;
        a *= a;
        a *= a;
        a *= a;
        a *= a;
        a *= a;
        exp *= a;
        a *= a;
        a *= a;
        a *= a;
        exp *= a;
        a *= a;
        a *= a;
        exp *= a;
        a *= a;
        exp *= a;
        a *= a;
        a *= a;
        a *= a;
        exp *= a;
        a *= a;
        a *= a;
        exp *= a;
        a *= a;
        a *= a;
        a *= a;
        exp *= a;
        a *= a;
        a *= a;
        a *= a;
        a *= a;
        exp *= a;
        a *= a;
        a *= a;
        exp *= a;
        a *= a;
        exp *= a;
        a *= a;
        a *= a;
        exp *= a;
        a *= a;
        a *= a;
        exp *= a;
        a *= a;
        a *= a;
        a *= a;
        exp *= a;
        a *= a;
        a *= a;
        a *= a;
        exp *= a;
        a *= a;
        exp *= a;
        a *= a;
        a *= a;
        a *= a;
        exp *= a;
        a *= a;
        a *= a;
        exp *= a;
        a *= a;
        exp *= a;
        a *= a;
        exp *= a;
        a *= a;
        a *= a;
        a *= a;
        exp *= a;
        a *= a;
        a *= a;
        a *= a;
        a *= a;
        exp *= a;
        a *= a;

        let first_frobenius = self.frobenius_map();
        let second_frobenius = first_frobenius.frobenius_map();
        let third_frobenius = second_frobenius.frobenius_map();

        exp * first_frobenius * second_frobenius.invert() * third_frobenius
    }


    /// Performs a sparse multiplication between `Fp12` and three elements of `Fp2`.
    ///
    /// > TODO: migrate to Ell data type instead.
    pub(crate) fn mul_by_014(self: Self, c0: Fp2, c1: Fp2, c4: Fp2) -> Fp12 {
        let aa = self.c0.mul_by_01(c0, c1);
        let bb = self.c1.mul_by_1(c4);
        let o = c1 + c4;
        let mut c1 = self.c1 + self.c0;
        c1 = c1.mul_by_01(c0, o);
        c1 = c1.sub(aa).sub(bb);
        let mut c0 = bb;
        c0 = c0.mul_by_nonresidue();
        c0 = c0 + aa;

        Fp12 { c0, c1 }
    }


    pub fn sparse_mul(self, ell: Ell) -> Fp12 {
        let aa = self.c0.mul_by_01(ell.c0_c0, ell.c1_c0);
        let bb = self.c1.mul_by_1(ell.c1_c1);
        let o = ell.c1_c0.add(ell.c1_c1);

        let mut c1 = self.c1.add(self.c0);
        c1 = c1.mul_by_01(ell.c0_c0, o);
        c1 = c1.sub(aa).sub(bb);

        let mut c0 = bb;
        c0 = c0.mul_by_nonresidue();
        c0 = c0.add(aa);

        Fp12 { c0, c1 }
    }


    // an unconstrained function that returns the third root of the element 
    pub unconstrained fn __third_root(self: Self) -> Self {
        let limbs: [[Field; 3]; 12] = __third_root_oracle(self); 
        Self::from_coeff(limbs)
    }

    pub unconstrained fn __is_third_root(self: Self) -> bool {
        let out = unsafe{__is_third_root_oracle(self)};
        out
    }

    pub unconstrained fn __random_third_root() -> Self {
        let limbs = unsafe{__random_third_root_oracle()};
        Self::from_coeff(limbs)
    }
}

impl std::ops::Add for Fp12 {
    fn add(self, other: Self) -> Self {
        Self { c0: self.c0 + other.c0, c1: self.c1 + other.c1 }
    }
}

impl std::ops::Sub for Fp12 {
    fn sub(self: Self, other: Self) -> Self {
        Self { c0: self.c0.sub(other.c0), c1: self.c1.sub(other.c1) }
    }
}

impl std::ops::Mul for Fp12 {
    fn mul(self, other: Self) -> Self {
        let a0 = self.c0 * other.c0;
        let a1 = self.c0 * other.c1 + self.c1 * other.c0;
        let a2 = self.c1 * other.c1;

        let c0 = a0 + a2.mul_by_nonresidue();
        let c1 = a1;
        // Return the final Fp12 structure
        Self { c0, c1 }
    }
}

impl std::cmp::Eq for Fp12 {
    fn eq(self, other: Self) -> bool {
        (self.c0 == other.c0) & (self.c1 == other.c1)
    }
}

#[oracle(third_root)]
pub unconstrained fn __third_root_oracle(a: Fp12) -> [[Field; 3];12] {
}


#[oracle(random_third_root)]
pub unconstrained fn __random_third_root_oracle() -> [[Field; 3];12] {
}

#[oracle(is_third_root)]
pub unconstrained fn __is_third_root_oracle(a: Fp12) -> bool {
}



// the oracle tests ==================================================
#[test]
fn test_random_third_root() {
    let out = unsafe{Fp12::__random_third_root()}; 
    assert(unsafe{out.__is_third_root() == true});
}

#[test]
fn test_third_root() {
    let a = unsafe{Fp12::__random_third_root()};
    let out = unsafe{a.__third_root()};
    assert((out* out * out) == a);
}



#[test]
fn test_is_third_root() {
    let a = Fp12::one();
    let out = unsafe{a.__is_third_root()};
    assert(out);
}

//==================================================


// #[test]
// fn test_arithmetic() {
//     let a = Fp12 {
//         c0: Fp6 {
//             c0: Fp2::new(
//                 [
//                     0xb3af6c44e48eb9b7d1b4b1c41e97b8,
//                     0x2b689fdd5da6619b2e16bf2482973d,
//                     0x4eb1b4adb3916c37153fa843bc1724,
//                     0x1ddb59,
//                 ],
//                 [
//                     0x687898c89dac3938647a0337390e91,
//                     0xfb34d9652fc9a5b7632f105f949e70,
//                     0xd4f2477096ad91b69a7a646bec7844,
//                     0xed52c,
//                 ],
//             ),
//             c1: Fp2::new(
//                 [
//                     0x8d42a6eee24c21ad6d25f3b895a378,
//                     0xaf34b84ed0165f4bc59f666809a305,
//                     0xa9a43882ae8a5cba9a9f45ace2559c,
//                     0x1a2632,
//                 ],
//                 [
//                     0x836efb730e4bf48a0ec0a552727f8f,
//                     0x6de2efd5fb2c80577e407d56e3584d,
//                     0x10275ac5fff269bdd4506932e43405,
//                     0x72a61,
//                 ],
//             ),
//             c2: Fp2::new(
//                 [
//                     0x85d961dd2944c554f41dbdfa55dc3e,
//                     0x94dd5f8f62f58e54520ec7ab9f7bdc,
//                     0xd76c86b02647531ee625f52f154f17,
//                     0x12bebd,
//                 ],
//                 [
//                     0x2f83f04cbabb9a42140b1cb7871cff,
//                     0x2d6ae61a88834a22708ea42961c300,
//                     0xcdd27a7dae8e0f9a9612d0a86d41d6,
//                     0x915aa,
//                 ],
//             ),
//         },
//         c1: Fp6 {
//             c0: Fp2::new(
//                 [
//                     0x175068bbfe0b1a2bbff8b7a5351255,
//                     0x1fa519cb5837aa2e1c4671ef9d5374,
//                     0xa3023daea697df7b912350025858bd,
//                     0x1676c5,
//                 ],
//                 [
//                     0x89a7261e177c96b86bb084e8844c12,
//                     0x23900e93f79297dc34656842864f9f,
//                     0x95b4e58f05d199c2c7760527c09bbe,
//                     0x16b0ae,
//                 ],
//             ),
//             c1: Fp2::new(
//                 [
//                     0x8350cdccab87c3593f38f0df87efc1,
//                     0x8925cda7f91b443d571d7e2bbf8b10,
//                     0x6f39338c7c39fb7434af47f709ae2b,
//                     0xe95c5,
//                 ],
//                 [
//                     0xa5ff1ac8f2aa1d45a9bd44f53561ce,
//                     0xe52df8692146afee27747a9ec2d33e,
//                     0xe3165b35d21f915f3332d58ce673c7,
//                     0x12cfd8,
//                 ],
//             ),
//             c2: Fp2::new(
//                 [
//                     0x38fb1613304b1b18a981ee7ea62d12,
//                     0xce7d54f7127829d3b51d6a0d90c481,
//                     0xd46daeed9ca25aa97c964fe9d14b69,
//                     0x163a1a,
//                 ],
//                 [
//                     0x823f42cc9582b3199546dbec99ad49,
//                     0xff1210846645f54b39d0cfc0742bd9,
//                     0x5b0e86ab5bdd5ada02591ced1e0280,
//                     0x19ece0,
//                 ],
//             ),
//         },
//     };

//     let b = Fp12 {
//         c0: Fp6 {
//             c0: Fp2::new(
//                 [
//                     0x6418e17bb08a26d334754a88af9ed8,
//                     0x5c4d4bf9a95873660b982f3b389486,
//                     0xe6b9301b716296eceb90902265e9c2,
//                     0x344eb,
//                 ],
//                 [
//                     0x2c4ff8f95c0d905b6460185b3881b2,
//                     0x48281080ec95757ca0e3cd5ce42ddf,
//                     0xe53b49ce6208454794fc46ac3383f3,
//                     0x4c61f,
//                 ],
//             ),
//             c1: Fp2::new(
//                 [
//                     0x1a0f38f91126c4385f50f0d87a3bf9,
//                     0xdfb2b1a80a69779abdccea85ff97c9,
//                     0xe9fc25a78a348ae518826015446871,
//                     0x1c49d1,
//                 ],
//                 [
//                     0x640bec9596ee7f23bf29bda216eb07,
//                     0x673020119d960490f956dd027725dd,
//                     0xb27c4fe60439d127a1c7a41318dcad,
//                     0x1bd151,
//                 ],
//             ),
//             c2: Fp2::new(
//                 [
//                     0xba255b575ff9286c9388bca4019c61,
//                     0x4c51a6b0b47e7d7598a652f0a73bd8,
//                     0xfdc9e2d0587d113f13f4826beb11f9,
//                     0x1761ea,
//                 ],
//                 [
//                     0xea9ac885f1c1478407e854749e0022,
//                     0xac5c8ba611a851b9d14acd557e2174,
//                     0x37319d961520599c35245529ac93b5,
//                     0x11f751,
//                 ],
//             ),
//         },
//         c1: Fp6 {
//             c0: Fp2::new(
//                 [
//                     0xd12f1144e528d5cce5f7eb98ea5d43,
//                     0xad4ce3734b537475320ace5bf8e001,
//                     0xc0be42a2a9066aa38e777bbd4497a8,
//                     0x10951b,
//                 ],
//                 [
//                     0x270a62aab9b6a4fa161969069340a,
//                     0x71b7c596b5fe7ea5a04c14f539c618,
//                     0xc0e27d563e716251e497b1851e5443,
//                     0x1c3fae,
//                 ],
//             ),
//             c1: Fp2::new(
//                 [
//                     0xb1b0ff6adab8c07a0309ac928eab17,
//                     0x202bb2a36371c30dee11befaad53ec,
//                     0x5c6128c31815fa8f11adff21bc0e93,
//                     0x5931,
//                 ],
//                 [
//                     0x38340422d00d5a2fd80c18732c1166,
//                     0xe56145364add39e1433e7666bec7ff,
//                     0x7bd111cb2e1ffd57621bfa33f6baec,
//                     0x3e8ca,
//                 ],
//             ),
//             c2: Fp2::new(
//                 [
//                     0x9951a3096171aa6b2b8f17ab3f4e97,
//                     0x4c6a214832ad14421d6205a401a29,
//                     0x746a9493f6e699f74f16cb44af029,
//                     0x1e06a7,
//                 ],
//                 [
//                     0x41bca809ec9e7f1fc87ac6e3a3cce9,
//                     0xb8a9fd6b10a1fa0a41984d71297189,
//                     0x2606b9b5501e92e16fba99ef7b27aa,
//                     0x143564,
//                 ],
//             ),
//         },
//     };

//     let c = Fp12 {
//         c0: Fp6 {
//             c0: Fp2::new(
//                 [
//                     0x720c1a648a6dc4377626d77b03a764,
//                     0x6e3db0d99300ef2265f0995ef97dc,
//                     0x435d78cbd6c7416bb0370f21ff21cf,
//                     0x13a44e,
//                 ],
//                 [
//                     0x84bfec07a5800c9c7c55d23ec6671b,
//                     0x2fdc513cdeffe29adcc3cb947507b5,
//                     0x76eb72e11ce785c99eeb4fb5d8040b,
//                     0x121e34,
//                 ],
//             ),
//             c1: Fp2::new(
//                 [
//                     0xbe7fbf6e2353e21c543e4bc02e0a9c,
//                     0xf46431d46c4a56f11e303698d2eebf,
//                     0xd496eb7c1c4f6e1d05f6939c63a421,
//                     0xd0b56,
//                 ],
//                 [
//                     0xcd31244bc979628056b1d106ecc8d8,
//                     0xfdaf068872d24f31805084a57f6c62,
//                     0x62991a2ebf541d9d6ef7c91dc50f97,
//                     0x14ab7a,
//                 ],
//             ),
//             c2: Fp2::new(
//                 [
//                     0x5f048867ad692fa8210d94fabfa2fd,
//                     0x2e62e766035568f1a9721a764287ab,
//                     0xbab9bfcac9c257038f6a9f9396dbb2,
//                     0x1da83b,
//                 ],
//                 [
//                     0x53ac5d922508c72ca4362d04e32dfa,
//                     0x166c20653130bad0fb9cad59fc12fc,
//                     0x655683e232b0eb233b495982cc79f0,
//                     0x42738,
//                 ],
//             ),
//         },
//         c1: Fp6 {
//             c0: Fp2::new(
//                 [
//                     0x571bacfab9112ef16ee344f0c16268,
//                     0x7a3f5e2e9db974f0d87ea80fa1d428,
//                     0xb0104e0b6a7ef0c318eacf111802b0,
//                     0xa9736,
//                 ],
//                 [
//                     0xe41318acb7a8c0e6a3458d84fd1c44,
//                     0xaa72d9c749ee1ace125dbefcd38c5b,
//                     0xc3217370529d5d92bcb222724d76fd,
//                     0x1812f3,
//                 ],
//             ),
//             c1: Fp2::new(
//                 [
//                     0x7cd66ddcbf855f85e9692aeabff5d8,
//                     0xef4f784c404b42a9f697767fd1aaa7,
//                     0x6fdc056bd66cf95fedc12f069ba3e1,
//                     0x86933,
//                 ],
//                 [
//                     0x35205828d34eae50579abdc513cbee,
//                     0xa36ac47cdfab5279bd233f26346f95,
//                     0x1a3ab7932780db066622c871503529,
//                     0xe5913,
//                 ],
//             ),
//             c2: Fp2::new(
//                 [
//                     0xb1f4a2c5ab1fde30be6b78791d983c,
//                     0xf9cb088c4af65ede4cf0b54f5f5430,
//                     0xd49efbdab23637124bda021af5f7b8,
//                     0x1b50ba,
//                 ],
//                 [
//                     0x866177018d092807ccdc7d8e9a1ec7,
//                     0x69fcac9ef500eb8a513c902d90ad92,
//                     0x1cec57a784c5e8abfb1a4223227d2f,
//                     0x74c86,
//                 ],
//             ),
//         },
//     };

//     let a_square = a.square();
//     let a_times_a = a.mul(a);
//     assert(a_square.eq(a_times_a));

//     let b_square = b.square();
//     let b_times_b = b.mul(b);
//     assert(b_square.eq(b_times_b));

//     // (a+b)*c = a*c + b*c
//     let lhs = (a.add(b)).mul(c);
//     let rhs = (a.mul(c)).add(b.mul(c));
//     assert(lhs.eq(rhs));
// }

// // #[test]
// fn test_invert() {
//     let a = Fp12 {
//         c0: Fp6 {
//             c0: Fp2::new(
//                 [
//                     0xcd5d36546bea81490d5309111ddceb,
//                     0x60fb0821e8477be42000281597c30d,
//                     0xc5bceea125200ead53b9846d34e496,
//                     0xe6106,
//                 ],
//                 [
//                     0x8ebdbdbc10648d74bccd3a998f7d14,
//                     0xbc1f6ba6171b9872df351363aff1fe,
//                     0xae277a7948eb47594c31af60354995,
//                     0x138b40,
//                 ],
//             ),
//             c1: Fp2::new(
//                 [
//                     0x93de59c75d21d5a02ed285a867278d,
//                     0x68743f57ffe21ee7ce991a169d868c,
//                     0xc5eaf10f2f2c5b130bafcb25ceee18,
//                     0x14cf6e,
//                 ],
//                 [
//                     0x27d65314e3e04fe2e42598308f46f6,
//                     0xcc9f7ad74ebbce368ea95ae986a3dc,
//                     0x891b4de53b5c393493ba1f272679ad,
//                     0x143ce9,
//                 ],
//             ),
//             c2: Fp2::new(
//                 [
//                     0x5589d21d49d8839b50d88f7990521a,
//                     0xfb74891851fe5934d176112797d7af,
//                     0x7b93b06afebc4aa2b67b0c035d2567,
//                     0x91a36,
//                 ],
//                 [
//                     0xfe74e25003c99c8d3ed0d02c49d9bc,
//                     0xdb28b3c2cf0f4cfbad2e238c6c6c0d,
//                     0xafac7905021bc99cfb6dca725c7127,
//                     0x7ef0f,
//                 ],
//             ),
//         },
//         c1: Fp6 {
//             c0: Fp2::new(
//                 [
//                     0x364af51cbd5ff7e6819c1d34723deb,
//                     0xbf5560c1cad36ff82b7ae6b4e9b669,
//                     0xb8b43b7067d5aebbfaf30af941fb6b,
//                     0x5d606,
//                 ],
//                 [
//                     0xd6283df9238fd8892faacf218aab13,
//                     0x778fff2b53389c2467832ee0c98fae,
//                     0xbe763b594c7de966e4aee7235e03f1,
//                     0xded4f,
//                 ],
//             ),
//             c1: Fp2::new(
//                 [
//                     0xb1b1dc8139b0ea6c55d4787a32c7fd,
//                     0xb0f3f1c9b39a58c54a1d41048049b6,
//                     0xa721d13a399308f0e7f083ded43cb8,
//                     0x1c010c,
//                 ],
//                 [
//                     0x973a2a6b834b9ec30d3fa2125d635b,
//                     0xcb60de306415b58bb2b74169b3b22d,
//                     0x2d5df2d5a06c53ff0f6738cd85426e,
//                     0x15dbd0,
//                 ],
//             ),
//             c2: Fp2::new(
//                 [
//                     0x3a5579ef31cee3bf72ede9a4ab4751,
//                     0xd6e693b6dc6b8074b6832bc70ba85,
//                     0x6d2789051356b780ff5ff6d005f5dc,
//                     0x1460d,
//                 ],
//                 [
//                     0xd2c6fb516a060d3f63a44ecc0b8532,
//                     0x4d8aea5422fec4e9c93bd8b65553e9,
//                     0xd7333204133e13cdb7faf3bcab8632,
//                     0x187373,
//                 ],
//             ),
//         },
//     };

//     let a_inv = a.invert();
//     let check = a.mul(a_inv);
//     assert(check.eq(Fp12::one()));
// }

#[test]
fn test_one() {
    let a = Fp12 {
        c0: Fp6 {
            c0: Fp2::new(
                [281382237812545379193668594878174908, 90824900576507601223966148779479394, 1640],
                [1199704082799405279767970336007399382, 105838412715217748738343588203455976, 3470],
            ),
            c1: Fp2::new(
                [699517922600938094581844456622821606, 766244793208285317105753995770523332, 5391],
                [1174771148458207518820683564604158404, 890744614640420316156493451275032886, 4653],
            ),
            c2: Fp2::new(
                [489655109890455545367468964680482035, 398027084094548452922359112094161148, 6728],
                [496665123426422059639752844234206160, 350210767051572081629044293495878811, 1877],
            ),
        },
        c1: Fp6 {
            c0: Fp2::new(
                [331816679168596495947305917667961895, 1126765732926166911585746373099035494, 3978],
                [18371723599779598995506163868555812, 32240751050272653558476274846155484, 6277],
            ),
            c1: Fp2::new(
                [1278358882567528155089415933911458643, 32628173174439344556782608871133808, 5767],
                [
                    920401518443459062035096401025158681,
                    1092237696507406181510465998476272124,
                    11505,
                ],
            ),
            c2: Fp2::new(
                [638485774760063632257790124842757197, 235056720649876689441574113803220985, 6647],
                [740256602156068976136086134195702814, 1157233779917845386103985573670726856, 9415],
            ),
        },
    };
    assert((a * Fp12::one()).eq(a));
}

#[test]
fn test_invert() {
    let a = Fp12 {
        c0: Fp6 {
            c0: Fp2::new(
                [281382237812545379193668594878174908, 90824900576507601223966148779479394, 1640],
                [1199704082799405279767970336007399382, 105838412715217748738343588203455976, 3470],
            ),
            c1: Fp2::new(
                [699517922600938094581844456622821606, 766244793208285317105753995770523332, 5391],
                [1174771148458207518820683564604158404, 890744614640420316156493451275032886, 4653],
            ),
            c2: Fp2::new(
                [489655109890455545367468964680482035, 398027084094548452922359112094161148, 6728],
                [496665123426422059639752844234206160, 350210767051572081629044293495878811, 1877],
            ),
        },
        c1: Fp6 {
            c0: Fp2::new(
                [331816679168596495947305917667961895, 1126765732926166911585746373099035494, 3978],
                [18371723599779598995506163868555812, 32240751050272653558476274846155484, 6277],
            ),
            c1: Fp2::new(
                [1278358882567528155089415933911458643, 32628173174439344556782608871133808, 5767],
                [
                    920401518443459062035096401025158681,
                    1092237696507406181510465998476272124,
                    11505,
                ],
            ),
            c2: Fp2::new(
                [638485774760063632257790124842757197, 235056720649876689441574113803220985, 6647],
                [740256602156068976136086134195702814, 1157233779917845386103985573670726856, 9415],
            ),
        },
    };
    let b: Fp12 = a.invert();
    let check = a.mul(b);
    assert(check.eq(Fp12::one()));
}

#[test]
fn test_square() {
    let a = Fp12 {
        c0: Fp6 {
            c0: Fp2::new(
                [281382237812545379193668594878174908, 90824900576507601223966148779479394, 1640],
                [1199704082799405279767970336007399382, 105838412715217748738343588203455976, 3470],
            ),
            c1: Fp2::new(
                [699517922600938094581844456622821606, 766244793208285317105753995770523332, 5391],
                [1174771148458207518820683564604158404, 890744614640420316156493451275032886, 4653],
            ),
            c2: Fp2::new(
                [489655109890455545367468964680482035, 398027084094548452922359112094161148, 6728],
                [496665123426422059639752844234206160, 350210767051572081629044293495878811, 1877],
            ),
        },
        c1: Fp6 {
            c0: Fp2::new(
                [331816679168596495947305917667961895, 1126765732926166911585746373099035494, 3978],
                [18371723599779598995506163868555812, 32240751050272653558476274846155484, 6277],
            ),
            c1: Fp2::new(
                [1278358882567528155089415933911458643, 32628173174439344556782608871133808, 5767],
                [
                    920401518443459062035096401025158681,
                    1092237696507406181510465998476272124,
                    11505,
                ],
            ),
            c2: Fp2::new(
                [638485774760063632257790124842757197, 235056720649876689441574113803220985, 6647],
                [740256602156068976136086134195702814, 1157233779917845386103985573670726856, 9415],
            ),
        },
    };
    let res = a.square();
    assert(res.eq(a * a));
}

#[test]
fn test_frobenius_map() {
    let a = Fp12 {
        c0: Fp6 {
            c0: Fp2::new(
                [281382237812545379193668594878174908, 90824900576507601223966148779479394, 1640],
                [1199704082799405279767970336007399382, 105838412715217748738343588203455976, 3470],
            ),
            c1: Fp2::new(
                [699517922600938094581844456622821606, 766244793208285317105753995770523332, 5391],
                [1174771148458207518820683564604158404, 890744614640420316156493451275032886, 4653],
            ),
            c2: Fp2::new(
                [489655109890455545367468964680482035, 398027084094548452922359112094161148, 6728],
                [496665123426422059639752844234206160, 350210767051572081629044293495878811, 1877],
            ),
        },
        c1: Fp6 {
            c0: Fp2::new(
                [331816679168596495947305917667961895, 1126765732926166911585746373099035494, 3978],
                [18371723599779598995506163868555812, 32240751050272653558476274846155484, 6277],
            ),
            c1: Fp2::new(
                [1278358882567528155089415933911458643, 32628173174439344556782608871133808, 5767],
                [
                    920401518443459062035096401025158681,
                    1092237696507406181510465998476272124,
                    11505,
                ],
            ),
            c2: Fp2::new(
                [638485774760063632257790124842757197, 235056720649876689441574113803220985, 6647],
                [740256602156068976136086134195702814, 1157233779917845386103985573670726856, 9415],
            ),
        },
    };
    let res = a
        .frobenius_map()
        .frobenius_map()
        .frobenius_map()
        .frobenius_map()
        .frobenius_map()
        .frobenius_map()
        .frobenius_map()
        .frobenius_map()
        .frobenius_map()
        .frobenius_map()
        .frobenius_map()
        .frobenius_map();

    // println(res);
    // println(a);
    assert(res.eq(a));
}

#[test]
fn test_exp_lambda() {
    // let a = Fp12::one();
    // let res = a.exp_lambda();
    // assert(res.eq(a));
    let a = Fp12 {
        c0: Fp6 {
            c0: Fp2::new(
                [281382237812545379193668594878174908, 90824900576507601223966148779479394, 1640],
                [1199704082799405279767970336007399382, 105838412715217748738343588203455976, 3470],
            ),
            c1: Fp2::new(
                [699517922600938094581844456622821606, 766244793208285317105753995770523332, 5391],
                [1174771148458207518820683564604158404, 890744614640420316156493451275032886, 4653],
            ),
            c2: Fp2::new(
                [489655109890455545367468964680482035, 398027084094548452922359112094161148, 6728],
                [496665123426422059639752844234206160, 350210767051572081629044293495878811, 1877],
            ),
        },
        c1: Fp6 {
            c0: Fp2::new(
                [331816679168596495947305917667961895, 1126765732926166911585746373099035494, 3978],
                [18371723599779598995506163868555812, 32240751050272653558476274846155484, 6277],
            ),
            c1: Fp2::new(
                [1278358882567528155089415933911458643, 32628173174439344556782608871133808, 5767],
                [
                    920401518443459062035096401025158681,
                    1092237696507406181510465998476272124,
                    11505,
                ],
            ),
            c2: Fp2::new(
                [638485774760063632257790124842757197, 235056720649876689441574113803220985, 6647],
                [740256602156068976136086134195702814, 1157233779917845386103985573670726856, 9415],
            ),
        },
    };
    let res = a.exp_lambda();
    let b: Fp12 = a.invert();
    let res2: Fp12 = b.exp_lambda();
    assert((res * res2).eq(Fp12::one()));
}

