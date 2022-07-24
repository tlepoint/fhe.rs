//! Unsigned 256-bit integer.

#![no_std]

use std::ops::{Shr, ShrAssign};

#[repr(C)]
#[derive(Eq, PartialEq, Debug, Copy, Clone)]
pub struct U256([u64; 4]);

impl U256 {
	pub fn zero() -> U256 {
		U256([0; 4])
	}
	pub fn new(a: [u64; 4]) -> Self {
		U256(a)
	}
	pub fn overflowing_add(self, other: U256) -> (U256, bool) {
		let U256(ref me) = self;
		let U256(ref you) = other;

		let mut ret = [0u64; 4];
		let mut carry = false;
		for i in 0..4 {
			(ret[i], carry) = me[i].carrying_add(you[i], carry);
		}

		(U256(ret), carry)
	}
	pub fn overflowing_sub(self, other: U256) -> (U256, bool) {
		let U256(ref me) = self;
		let U256(ref you) = other;

		let mut ret = [0u64; 4];
		let mut borrow = false;
		for i in 0..4 {
			(ret[i], borrow) = me[i].borrowing_sub(you[i], borrow);
		}

		(U256(ret), borrow)
	}
}

impl From<U256> for u64 {
	fn from(val: U256) -> u64 {
		val.0[0]
	}
}

impl From<u64> for U256 {
	fn from(val: u64) -> U256 {
		U256([0, 0, 0, val])
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
			self.0[0] = self.0[3] >> (rhs - 192);
			self.0[1] = 0;
			self.0[2] = 0;
			self.0[3] = 0;
		} else if rhs > 128 {
			self.0[0] = (self.0[2] >> (rhs - 128)) | (self.0[3] << (192 - rhs));
			self.0[1] = self.0[3] >> (rhs - 128);
			self.0[2] = 0;
			self.0[3] = 0;
		} else if rhs == 128 {
			self.0[0] = self.0[2];
			self.0[1] = self.0[3];
			self.0[2] = 0;
			self.0[3] = 0;
		} else if rhs > 64 {
			self.0[0] = (self.0[1] >> (rhs - 64)) | (self.0[2] << (128 - rhs));
			self.0[1] = (self.0[2] >> (rhs - 64)) | (self.0[3] << (128 - rhs));
			self.0[2] = self.0[3] >> (rhs - 64);
			self.0[3] = 0;
		} else if rhs == 64 {
			self.0[0] = self.0[1];
			self.0[1] = self.0[2];
			self.0[2] = self.0[3];
			self.0[3] = 0;
		} else if rhs > 0 {
			self.0[0] = (self.0[0] >> rhs) | (self.0[1] << (64 - rhs));
			self.0[1] = (self.0[1] >> rhs) | (self.0[2] << (64 - rhs));
			self.0[2] = (self.0[2] >> rhs) | (self.0[3] << (64 - rhs));
			self.0[3] >>= rhs;
		}
	}
}
