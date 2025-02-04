use ark_bn254::Fr;
use ark_bn254::{Fq, Fq12, Fq2, Fq6};
use ark_ff::Zero;
use jsonrpsee::core::params;
use num_bigint::BigInt;
use num_bigint::BigUint;
use num_traits::Num;
use serde_json::{json, Value};
// use std::str::FromStr;

use crate::foreign_call::ForeignCallParam;
use crate::ops::witness_generator::WitnessGenerator;
use crate::ops::witness_generator::WitnessGeneratorTrait;

// a struct that emulates the bignum Params params from the params.nr file
struct Params {
    has_multiplicative_inverse: bool,
    modulus: BigUint,
    double_modulus: Vec<Fr>,
    redc_param: Vec<Fr>,
}

/**** THERE'S A LOT OF BOILERPLATE INSIDE THESE "HANDLERS", THAT WE CAN PROBABLY PUT INTO COMMON HELPER FUNCTIONS ****/

pub(crate) fn handle_witness_gen(inputs: &Vec<ForeignCallParam<String>>) -> Value {
    Value::String("Hello, world!".to_string())
}

pub(crate) fn handle_third_root(inputs: &Vec<ForeignCallParam<String>>) -> Value {
    // the input has 12 elements, each a bignum representing an FP element
    // now we cast the bignums to bigUints
    let fp12 = get_fq12_from_callparam(inputs);
    let result = WitnessGenerator::tonelli_shanks_third_root(fp12);
    let results_formatted = cast_fp12_to_noir_fp12(result);
    let return_vec: Vec<Vec<String>> = vec![results_formatted];
    let json_response = json!({"values" : return_vec});
    json_response
}

pub(crate) fn handle_is_third_root(inputs: &Vec<ForeignCallParam<String>>) -> Value {
    let fp12 = get_fq12_from_callparam(inputs);
    let result: bool = WitnessGenerator::is_third_root(&fp12);

    let as_big_uint: BigUint = result.into();
    let as_hex_str = as_big_uint.to_str_radix(16);
    let mut results: Vec<String> = vec![];
    results.push(as_hex_str);
    let oracle_return_data_the_noir_program_expects = results;

    /**** FORMAT RESULT FOR NOIR CONSUMPTION, AND CONVERT RESULT TO JSON `Value` TYPE ****/
    let return_vec = oracle_return_data_the_noir_program_expects; // Notice! This is a different type from the singular handle_get_sqrt function! Hence why the `Value` is being computed inside this function, instead in the calling function.

    let json_response = json!({"values" : return_vec});
    json_response
}

pub(crate) fn handle_random_third_root(inputs: &Vec<ForeignCallParam<String>>) -> Value {
    let result = WitnessGenerator::rand_third_root();
    let results_formatted = cast_fp12_to_noir_fp12(result);
    let return_vec: Vec<Vec<String>> = vec![results_formatted];
    let json_response = json!({"values" : return_vec});
    json_response
}

pub(crate) fn handle_get_pairing_witnesses(inputs: &Vec<ForeignCallParam<String>>) -> Value {
    let fp12 = get_fq12_from_callparam(inputs);
    let (c, u) = WitnessGenerator::witness_generator(fp12);
    let c_formatted = cast_fp12_to_noir_fp12(c);
    let u_formatted = cast_fp12_to_noir_fp12(u);
    let return_vec: Vec<Vec<String>> = vec![c_formatted, u_formatted];
    let json_response = json!({"values" : return_vec});
    json_response
}

pub fn get_fq12_from_callparam(inputs: &Vec<ForeignCallParam<String>>) -> Fq12 {
    let mut biguints: Vec<Fq> = vec![];
    for input in inputs {
        let biguint = cast_to_biguint(callparam_to_string(input));
        biguints.push(Fq::from(biguint));
    }
    // now cast these into an ark_ff::Fq12 element
    let fp12 = Fq12 {
        c0: Fq6 {
            c0: Fq2 {
                c0: biguints[0],
                c1: biguints[1],
            },
            c1: Fq2 {
                c0: biguints[2],
                c1: biguints[3],
            },
            c2: Fq2 {
                c0: biguints[4],
                c1: biguints[5],
            },
        },
        c1: Fq6 {
            c0: Fq2 {
                c0: biguints[6],
                c1: biguints[7],
            },
            c1: Fq2 {
                c0: biguints[8],
                c1: biguints[9],
            },
            c2: Fq2 {
                c0: biguints[10],
                c1: biguints[11],
            },
        },
    };
    fp12
}

pub(crate) fn cast_fp12_to_noir_fp12(input: Fq12) -> Vec<String> {
    // we have a 2 over 3 over 2 field element, so we need to extract the limbs
    let limbs: Vec<BigUint> = vec![
        input.c0.c0.c0.into(),
        input.c0.c0.c1.into(),
        input.c0.c1.c0.into(),
        input.c0.c1.c1.into(),
        input.c0.c2.c0.into(),
        input.c0.c2.c1.into(),
        input.c1.c0.c0.into(),
        input.c1.c0.c1.into(),
        input.c1.c1.c0.into(),
        input.c1.c1.c1.into(),
        input.c1.c2.c0.into(),
        input.c1.c2.c1.into(),
    ];

    let mut results_formatted: Vec<String> = vec![];
    for result in limbs {
        let limbs = cast_bigint_to_bignum_limbs(&result.into(), 3);
        for limb in limbs {
            results_formatted.push(limb);
        }
    }
    results_formatted
}

pub(crate) fn cast_to_biguint(input_strings: Vec<&str>) -> BigUint {
    // split the limbs
    let mut limbs: Vec<BigUint> = vec![];
    for input_string in input_strings {
        // handle the case of a zero input
        if input_string == "" {
            let x_big_uint = BigUint::from_str_radix("0", 16).unwrap();
            limbs.push(x_big_uint);
        } else {
            let x_big_uint: BigUint = BigUint::from_str_radix(input_string, 16).unwrap();
            limbs.push(x_big_uint);
        }
    }
    // a constant 2^120 as biguint
    let base: BigUint = BigUint::from(2u32);
    let exp = 120u32;
    let shift_constant = base.pow(exp);
    let mut res = BigUint::ZERO;
    for i in 0..limbs.len() {
        res = res + &limbs[i] * &shift_constant.pow(i as u32);
    }
    res
}

// helper function to get limbs of a big num and pack them into a vector of Fr elements
pub(crate) fn gets_limbs(input_strings: Vec<&str>) -> Vec<Fr> {
    let mut limbs: Vec<Fr> = vec![];
    for input_string in input_strings {
        // handle the case of a zero input
        if input_string == "" {
            let x_big_uint = BigUint::from_str_radix("0", 16).unwrap();
            let limb: Fr = x_big_uint.into();
            limbs.push(limb);
        } else {
            let x_big_uint: BigUint = BigUint::from_str_radix(input_string, 16).unwrap();
            let limb: Fr = x_big_uint.into();
            limbs.push(limb);
        }
    }
    limbs
}

pub(crate) fn callparam_to_string(input: &ForeignCallParam<String>) -> Vec<&str> {
    match input {
        ForeignCallParam::Single(value) => vec![value.trim_start_matches('0')],
        ForeignCallParam::Array(values) => values
            .into_iter()
            .map(|v| v.trim_start_matches('0'))
            .collect(),
    }
}

pub(crate) fn get_u32_from_callparam(input: &ForeignCallParam<String>) -> u32 {
    let input_string = callparam_to_string(input)[0];
    u32::from_str_radix(input_string, 16).unwrap()
}

pub(crate) fn get_bool_from_callparam(input: &ForeignCallParam<String>) -> bool {
    let mut input_string = callparam_to_string(input)[0];
    if input_string == "" {
        input_string = "0";
    }
    let res = u32::from_str_radix(input_string, 16).unwrap();
    res == 1
}

pub(crate) fn cast_biguint_to_bignum_limbs(input: &BigUint, num_limbs: u32) -> Vec<String> {
    // a constant 2^120 as biguint
    let base: BigUint = BigUint::from(2u32);
    let exp = 120u32;
    let shift_constant = base.pow(exp);
    let mut input_copy = input.clone();
    // an empty array of size num_limbs of type hex limbs
    let mut limbs_hex: Vec<String> = vec![];
    for i in 0..num_limbs {
        let remainder = &input_copy % &shift_constant;
        limbs_hex.push(remainder.to_str_radix(16));
        let quetient: BigUint = input_copy / &shift_constant;
        input_copy = quetient.clone();
    }
    limbs_hex
}

pub(crate) fn cast_bigint_to_bignum_limbs(input: &BigInt, num_limbs: u32) -> Vec<String> {
    // a constant 2^120 as biguint
    let base: BigInt = BigInt::from(2u32);
    let exp = 120u32;
    let shift_constant = base.pow(exp);
    let mut input_copy = input.clone();
    // an empty array of size num_limbs of type hex limbs
    let mut limbs_hex: Vec<String> = vec![];
    for i in 0..num_limbs {
        let remainder = &input_copy % &shift_constant;
        limbs_hex.push(remainder.to_str_radix(16));
        let quetient: BigInt = input_copy / &shift_constant;
        input_copy = quetient.clone();
    }
    limbs_hex
}

impl Params {
    // this function takes the foreign call params and returns a Params struct
    pub fn from_foreign_call_params(inputs: &Vec<ForeignCallParam<String>>) -> Params {
        let has_multiplicative_inverse_fc = &inputs[0];
        let has_multiplicative_inverse = get_bool_from_callparam(&has_multiplicative_inverse_fc);
        let modulus_fc = &inputs[1];
        // let modulus = gets_limbs(callparam_to_string(modulus_fc));
        let modulus_str = callparam_to_string(modulus_fc);
        let modulus = cast_to_biguint(modulus_str);
        let double_modulus_fc = &inputs[4];
        let double_modulus = gets_limbs(callparam_to_string(double_modulus_fc));
        let redc_param_fc = &inputs[5];
        let redc_param = gets_limbs(callparam_to_string(redc_param_fc));

        Params {
            has_multiplicative_inverse: has_multiplicative_inverse,
            modulus: modulus,
            double_modulus: double_modulus,
            redc_param: redc_param,
        }
    }
}

impl std::fmt::Debug for Params {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "Params {{ has_multiplicative_inverse: {:?}, modulus: {:?}, double_modulus: {:?}, redc_param: {:?}", self.has_multiplicative_inverse, self.modulus, self.double_modulus, self.redc_param)
    }
}
