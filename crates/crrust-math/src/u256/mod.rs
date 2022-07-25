//! Unsigned 256-bit integer.

use std::ops::{Not, Shr, ShrAssign};

#[repr(C)]
#[derive(Eq, PartialEq, Debug, Copy, Clone)]
pub struct U256(u64, u64, u64, u64);

impl U256 {
	pub fn zero() -> U256 {
		U256 {
			0: 0,
			1: 0,
			2: 0,
			3: 0,
		}
	}

	pub fn overflowing_add(&mut self, other: U256) {
		let mut _carry = false;
		(self.0, _carry) = self.0.carrying_add(other.0, _carry);
		(self.1, _carry) = self.1.carrying_add(other.1, _carry);
		(self.2, _carry) = self.2.carrying_add(other.2, _carry);
		(self.3, _carry) = self.3.carrying_add(other.3, _carry);
	}

	pub fn overflowing_sub(&mut self, other: U256) {
		let mut _borrow = false;
		(self.0, _borrow) = self.0.borrowing_sub(other.0, _borrow);
		(self.1, _borrow) = self.1.borrowing_sub(other.1, _borrow);
		(self.2, _borrow) = self.2.borrowing_sub(other.2, _borrow);
		(self.3, _borrow) = self.3.borrowing_sub(other.3, _borrow);
	}

	/// Converts to u64.
	///
	/// Aborts if larger than 2^64 in debug mode.
	pub fn as_u64(self) -> u64 {
		debug_assert!(self.1 == 0 && self.2 == 0 && self.3 == 0);
		self.0
	}

	/// Converts to u128.
	///
	/// Aborts if larger than 2^128 in debug mode.
	pub fn as_u128(self) -> u128 {
		debug_assert!(self.2 == 0 && self.3 == 0);
		(self.0 as u128) + ((self.1 as u128) << 64)
	}

	/// Returns the most significant bit of the unsigned integer.
	pub fn msb(self) -> u64 {
		self.3 >> 63
	}
}

impl From<[u64; 4]> for U256 {
	fn from(a: [u64; 4]) -> Self {
		Self {
			0: a[0],
			1: a[1],
			2: a[2],
			3: a[3],
		}
	}
}

impl From<U256> for [u64; 4] {
	fn from(a: U256) -> [u64; 4] {
		[a.0, a.1, a.2, a.3]
	}
}

impl Not for U256 {
	type Output = Self;

	fn not(self) -> Self {
		Self {
			0: !self.0,
			1: !self.1,
			2: !self.2,
			3: !self.3,
		}
	}
}

impl Shr<usize> for U256 {
	type Output = Self;

	fn shr(self, rhs: usize) -> Self {
		let mut r = self.clone();
		r >>= rhs;
		r
	}
}

impl ShrAssign<usize> for U256 {
	fn shr_assign(&mut self, rhs: usize) {
		debug_assert!(rhs < 256);

		if rhs >= 192 {
			self.0 = self.3 >> (rhs - 192);
			self.1 = 0;
			self.2 = 0;
			self.3 = 0;
		} else if rhs > 128 {
			self.0 = (self.2 >> (rhs - 128)) | (self.3 << (192 - rhs));
			self.1 = self.3 >> (rhs - 128);
			self.2 = 0;
			self.3 = 0;
		} else if rhs == 128 {
			self.0 = self.2;
			self.1 = self.3;
			self.2 = 0;
			self.3 = 0;
		} else if rhs > 64 {
			self.0 = (self.1 >> (rhs - 64)) | (self.2 << (128 - rhs));
			self.1 = (self.2 >> (rhs - 64)) | (self.3 << (128 - rhs));
			self.2 = self.3 >> (rhs - 64);
			self.3 = 0;
		} else if rhs == 64 {
			self.0 = self.1;
			self.1 = self.2;
			self.2 = self.3;
			self.3 = 0;
		} else if rhs > 0 {
			self.0 = (self.0 >> rhs) | (self.1 << (64 - rhs));
			self.1 = (self.1 >> rhs) | (self.2 << (64 - rhs));
			self.2 = (self.2 >> rhs) | (self.3 << (64 - rhs));
			self.3 >>= rhs;
		}
	}
}

#[cfg(test)]
mod tests {

	use super::U256;

	#[test]
	fn test_zero() {
		assert_eq!(U256::zero().as_u64(), 0u64);
	}

	proptest! {
		#[test]
		fn test_u64(a: u64) {
			prop_assert_eq!(a, U256::from([a, 0, 0, 0]).as_u64());
		}

		#[test]
		fn test_u128(a: u128) {
			prop_assert_eq!(a, U256::from([a as u64, (a >> 64) as u64, 0, 0]).as_u128());
		}

		#[test]
		fn test_from_into(a: u64, b: u64, c: u64, d:u64) {
			prop_assert_eq!(<[u64; 4]>::from(U256::from([a, b, c, d])), [a, b, c, d]);
		}

		#[test]
		fn test_shift(a: u64, b: u64, c: u64, d:u64, shift in 0..256usize) {
			prop_assert_eq!(<[u64; 4]>::from(U256::from([a, b, c, d]) >> 0), [a, b, c, d]);
			prop_assert_eq!(<[u64; 4]>::from(U256::from([a, b, c, d]) >> 64), [b, c, d, 0]);
			prop_assert_eq!(<[u64; 4]>::from(U256::from([a, b, c, d]) >> 128), [c, d, 0, 0]);
			prop_assert_eq!(<[u64; 4]>::from(U256::from([a, b, c, d]) >> 192), [d, 0, 0, 0]);

			prop_assume!(shift % 64 != 0);
			if shift < 64 {
				prop_assert_eq!(<[u64; 4]>::from(U256::from([a, b, c, d]) >> shift), [(a >> shift) | (b << (64 - shift)), (b >> shift) | (c << (64 - shift)), (c >> shift) | (d << (64 - shift)), d >> shift]);
			} else if shift < 128 {
				prop_assert_eq!(<[u64; 4]>::from(U256::from([a, b, c, d]) >> shift), [(b >> (shift - 64)) | (c << (128 - shift)), (c >> (shift - 64)) | (d << (128 - shift)), (d >> (shift - 64)), 0]);
			} else if shift < 192 {
				prop_assert_eq!(<[u64; 4]>::from(U256::from([a, b, c, d]) >> shift), [(c >> (shift - 128)) | (d << (192 - shift)), (d >> (shift - 128)), 0, 0]);
			} else {
				prop_assert_eq!(<[u64; 4]>::from(U256::from([a, b, c, d]) >> shift), [(d >> (shift - 192)), 0, 0, 0]);
			}
		}

		#[test]
		fn test_shift_assign(a: u64, b: u64, c: u64, d:u64, shift in 0..256usize) {
			let mut u = U256::from([a, b, c, d]);
			
			u >>= 0;
			prop_assert_eq!(<[u64; 4]>::from(u), [a, b, c, d]);
			
			u >>= 64;
			prop_assert_eq!(<[u64; 4]>::from(u), [b, c, d, 0]);

			u = U256::from([a, b, c, d]);
			u >>= 128;
			prop_assert_eq!(<[u64; 4]>::from(u), [c, d, 0, 0]);

			u = U256::from([a, b, c, d]);
			u >>= 192;
			prop_assert_eq!(<[u64; 4]>::from(u), [d, 0, 0, 0]);
			
			prop_assume!(shift % 64 != 0);
			u = U256::from([a, b, c, d]);
			u >>= shift;
			if shift < 64 {
				prop_assert_eq!(<[u64; 4]>::from(u), [(a >> shift) | (b << (64 - shift)), (b >> shift) | (c << (64 - shift)), (c >> shift) | (d << (64 - shift)), d >> shift]);
			} else if shift < 128 {
				prop_assert_eq!(<[u64; 4]>::from(u), [(b >> (shift - 64)) | (c << (128 - shift)), (c >> (shift - 64)) | (d << (128 - shift)), (d >> (shift - 64)), 0]);
			} else if shift < 192 {
				prop_assert_eq!(<[u64; 4]>::from(u), [(c >> (shift - 128)) | (d << (192 - shift)), (d >> (shift - 128)), 0, 0]);
			} else {
				prop_assert_eq!(<[u64; 4]>::from(u), [(d >> (shift - 192)), 0, 0, 0]);
			}
		}
	}
}
