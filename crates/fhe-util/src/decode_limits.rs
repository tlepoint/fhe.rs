//! Resource bounds shared by contextual protobuf import APIs.

/// Limits applied before decoding wire buffers or constructing arithmetic data.
/// These bound input and expanded polynomial storage, not allocator overhead or
/// the process's total RSS. Parameter cache size is also bounded by degree and
/// modulus count. Raise individual limits explicitly for larger trusted
/// workloads.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct DecodeLimits {
    /// Maximum bytes in a complete input message (default: 64 MiB).
    pub max_input_bytes: usize,
    /// Maximum polynomial degree, including parameter imports (default: 32768).
    pub max_degree: usize,
    /// Maximum RNS modulus count (default: 32).
    pub max_moduli: usize,
    /// Maximum encoded big plaintext modulus length (default: 1024 bytes).
    pub max_plaintext_bytes: usize,
    /// Maximum polynomial slots, including seeded and nested material (4096).
    pub max_polynomials: usize,
    /// Maximum charged expanded residue storage (default: 256 MiB).
    /// Each polynomial slot is charged two u64 arrays at the full supplied
    /// modulus chain, conservatively covering NTT/Shoup storage at lower
    /// levels.
    pub max_residue_bytes: usize,
}

impl Default for DecodeLimits {
    fn default() -> Self {
        Self {
            max_input_bytes: 64 << 20,
            max_degree: 32768,
            max_moduli: 32,
            max_plaintext_bytes: 1024,
            max_polynomials: 4096,
            max_residue_bytes: 256 << 20,
        }
    }
}

/// A wire resource exceeded its configured bound.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[non_exhaustive]
pub struct DecodeLimitError {
    /// Name of the exceeded resource.
    pub resource: &'static str,
    /// Requested count; overflow is reported as `usize::MAX`.
    pub actual: usize,
    /// Configured upper bound.
    pub maximum: usize,
}
impl std::fmt::Display for DecodeLimitError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "decode {} {} exceeds limit {}",
            self.resource, self.actual, self.maximum
        )
    }
}
impl std::error::Error for DecodeLimitError {}

impl DecodeLimits {
    /// Check a named resource before allocation.
    pub fn check(
        &self,
        resource: &'static str,
        actual: usize,
        maximum: usize,
    ) -> Result<(), DecodeLimitError> {
        if actual > maximum {
            Err(DecodeLimitError {
                resource,
                actual,
                maximum,
            })
        } else {
            Ok(())
        }
    }

    /// Check input size and a supplied or decoded arithmetic context.
    pub fn check_context(
        &self,
        input_bytes: usize,
        degree: usize,
        moduli: usize,
    ) -> Result<(), DecodeLimitError> {
        self.check("input bytes", input_bytes, self.max_input_bytes)?;
        self.check("degree", degree, self.max_degree)?;
        self.check("moduli", moduli, self.max_moduli)
    }

    /// Check polynomial count and conservatively charged residue storage.
    pub fn check_polynomials(
        &self,
        count: usize,
        degree: usize,
        moduli: usize,
    ) -> Result<(), DecodeLimitError> {
        self.check("polynomials", count, self.max_polynomials)?;
        let Some(bytes) = count
            .checked_mul(degree)
            .and_then(|n| n.checked_mul(moduli))
            .and_then(|n| n.checked_mul(16))
        else {
            return Err(DecodeLimitError {
                resource: "residue size overflow",
                actual: usize::MAX,
                maximum: self.max_residue_bytes,
            });
        };
        self.check("residue bytes", bytes, self.max_residue_bytes)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn limits_include_boundaries_and_reject_overflow() {
        let limits = DecodeLimits {
            max_polynomials: 2,
            max_residue_bytes: 1024,
            ..DecodeLimits::default()
        };
        assert!(limits.check_polynomials(2, 16, 2).is_ok());
        assert_eq!(
            limits.check_polynomials(3, 16, 2).unwrap_err().resource,
            "polynomials"
        );
        assert_eq!(
            limits.check_polynomials(2, 17, 2).unwrap_err().resource,
            "residue bytes"
        );
        let huge = DecodeLimits {
            max_polynomials: usize::MAX,
            max_residue_bytes: usize::MAX,
            ..limits
        };
        assert!(huge.check_polynomials(usize::MAX, usize::MAX, 1).is_err());
        assert!(
            limits
                .check_context(limits.max_input_bytes + 1, 16, 2)
                .is_err()
        );
        assert!(limits.check_context(0, limits.max_degree + 1, 2).is_err());
        assert!(limits.check_context(0, 16, limits.max_moduli + 1).is_err());
    }
}
