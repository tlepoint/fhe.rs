//! Unsigned 256-bit integer.

use std::ops::{Shr, ShrAssign};

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
	pub fn new(a: [u64; 4]) -> Self {
		U256 {
			0: a[0],
			1: a[1],
			2: a[2],
			3: a[3],
		}
	}
	pub fn overflowing_add(&mut self, other: U256) {
		let mut carry = false;
		(self.0, carry) = self.0.carrying_add(other.0, carry);
		(self.1, carry) = self.1.carrying_add(other.1, carry);
		(self.2, carry) = self.2.carrying_add(other.2, carry);
		(self.3, _) = self.3.carrying_add(other.3, carry);
	}
	pub fn overflowing_sub(&mut self, other: U256) {
		let mut borrow = false;
		(self.0, borrow) = self.0.borrowing_sub(other.0, borrow);
		(self.1, borrow) = self.1.borrowing_sub(other.1, borrow);
		(self.2, borrow) = self.2.borrowing_sub(other.2, borrow);
		(self.3, _) = self.3.borrowing_sub(other.3, borrow);
	}
}

impl From<U256> for u64 {
	fn from(val: U256) -> u64 {
		val.0
	}
}

impl From<u64> for U256 {
	fn from(val: u64) -> U256 {
		U256 {
			0: val,
			1: 0,
			2: 0,
			3: 0,
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
