pub mod fp2;
pub mod fp6;
pub mod fp12;
pub mod g1;
pub mod g2;
pub mod pairing;
mod test;

// export at the top level:
pub use pairing::pairing;
pub use g1::G1Affine;
pub use g2::G2Affine;
pub use fp12::Fp12;
pub use fp6::Fp6;
pub use fp2::Fp2;
