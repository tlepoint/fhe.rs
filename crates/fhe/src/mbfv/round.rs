//! Traits and markers distinguishing different rounds of a protocol.

/// Indicates that a type marks a particular round.
pub trait Round: sealed::Sealed + std::fmt::Debug + Clone + Eq {
    /// Round-specific dependency stored in a relinearization share. The second
    /// round always retains its first-round aggregation.
    type RelinDependency: std::fmt::Debug + Clone + Eq;
}

/// Marks the shares produced in round 1
#[derive(Debug, PartialEq, Eq, Clone)]
pub struct R1;
/// Marks the aggregated shares from round 1
#[derive(Debug, PartialEq, Eq, Clone)]
pub struct R1Aggregated;
/// Marks the shares produced in round 2
#[derive(Debug, PartialEq, Eq, Clone)]
pub struct R2;

impl Round for R1 {
    type RelinDependency = ();
}
impl Round for R1Aggregated {
    type RelinDependency = ();
}
impl Round for R2 {
    type RelinDependency = std::sync::Arc<super::RelinKeyShare<R1Aggregated>>;
}

mod sealed {
    pub trait Sealed {}
    impl Sealed for super::R1 {}
    impl Sealed for super::R1Aggregated {}
    impl Sealed for super::R2 {}
}
