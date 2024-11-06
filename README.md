# Pairing over BLS12-381

Implementation of pairing over BLS12-381 in Noir. This uses the new [BigNum](https://github.com/noir-lang/noir-bignum) library. 

## Implementation

This codebase uses the zkcrypto repository as a referece: https://github.com/zkcrypto/bls12_381. This has also been used for tests and test values. 
The Noir implementation passes pairing & bilinearity tests obtained from the zkcrypto repo. 

## Test

Run all tests
```
nargo test
```

Note that a good amount of the tests are commented out because they take a fair amount of time (20-30 min) to run. For example `test_pairings_1` and `test_bilinearity` in `pairings.nr`. 

## Related Noir work
- BLS12_381 Elliptic Curve Pairing and Signature Verification Library by @onurinanc: [repo](https://github.com/onurinanc/noir-bls-signature).
- Noir BigCurve library: [repo](https://github.com/noir-lang/noir_bigcurve).




