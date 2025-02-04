// This is basically a copy of noir/acvm-repo/brillig/src/foreign_call.rs

use serde::{de::value, Deserialize, Serialize};

/// Single output of a [foreign call][crate::Opcode::ForeignCall].
#[derive(Debug, PartialEq, Eq, Serialize, Deserialize, Clone)]
#[serde(untagged)]
pub enum ForeignCallParam<F> {
    Single(F),
    Array(Vec<F>),
}

impl<F> From<F> for ForeignCallParam<F> {
    fn from(value: F) -> Self {
        ForeignCallParam::Single(value)
    }
}

impl<F> From<Vec<F>> for ForeignCallParam<F> {
    fn from(values: Vec<F>) -> Self {
        ForeignCallParam::Array(values)
    }
}

impl<F> ForeignCallParam<F>
where
    F: Clone,
{
    pub fn len(&self) -> usize {
        match self {
            ForeignCallParam::Single(_) => 1,
            ForeignCallParam::Array(values) => values.len(),
        }
    }

    pub fn get_values(&self) -> Vec<F> {
        match self {
            ForeignCallParam::Single(value) => vec![value.clone()],
            ForeignCallParam::Array(values) => values.clone(),
        }
    }
}
