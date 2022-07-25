#![crate_name = "crrust_util"]
#![crate_type = "lib"]
#![warn(missing_docs, unused_imports)]

//! Utilities for the crrust library.

use num_bigint::{prime::probably_prime, BigUint};
use std::panic::UnwindSafe;

/// Define catch_unwind to silence the panic in unit tests.
pub fn catch_unwind<F, R>(f: F) -> std::thread::Result<R>
where
	F: FnOnce() -> R + UnwindSafe,
{
	let prev_hook = std::panic::take_hook();
	std::panic::set_hook(Box::new(|_| {}));
	let r = std::panic::catch_unwind(f);
	std::panic::set_hook(prev_hook);
	r
}

/// Returns whether the modulus p is prime; this function is 100% accurate.
pub fn is_prime(p: u64) -> bool {
	probably_prime(&BigUint::from(p), 0)
}
