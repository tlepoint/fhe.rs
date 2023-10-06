//! Traits and markers distinguishing different rounds of a protocol.

/// Indicates that a type marks a particular round.
pub trait Round: sealed::Sealed {}

/// Marks the shares produced in round 1
#[derive(Debug, PartialEq, Eq, Clone)]
pub struct R1;
/// Marks the aggregated shares from round 1
#[derive(Debug, PartialEq, Eq, Clone)]
pub struct R1Aggregated;
/// Marks the shares produced in round 2
#[derive(Debug, PartialEq, Eq, Clone)]
pub struct R2;

impl Round for R1 {}
impl Round for R1Aggregated {}
impl Round for R2 {}

mod sealed {
    pub trait Sealed {}
    impl Sealed for super::R1 {}
    impl Sealed for super::R1Aggregated {}
    impl Sealed for super::R2 {}
}
