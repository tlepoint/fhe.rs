//! Explicit timing and diagnostic permissions.

/// Evidence that the caller has classified data as public.
///
/// Constructing this value is an explicit assertion: using it for secret data
/// may expose information through timing, but does not violate Rust's memory
/// safety guarantees.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct PublicData(());

impl PublicData {
    /// Assert that the data involved in an operation is public.
    #[must_use]
    pub const fn assert_public() -> Self {
        Self(())
    }
}

/// Permission to use variable-time algorithms on data classified as public.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct VariableTime(PublicData);

impl VariableTime {
    /// Create variable-time permission from an explicit public-data assertion.
    #[must_use]
    pub const fn new(public_data: PublicData) -> Self {
        Self(public_data)
    }
}

/// Acknowledgment that a diagnostic may reveal secret-dependent information.
///
/// This is separate from public-data classification: both the returned result
/// and the running time of a diagnostic can disclose information about secrets.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct SecretDependentDiagnostics(());

impl SecretDependentDiagnostics {
    /// Accept secret-dependent output and timing leakage in a diagnostic
    /// setting. Do not expose such diagnostics to untrusted callers.
    #[must_use]
    pub const fn acknowledge_leakage() -> Self {
        Self(())
    }
}
