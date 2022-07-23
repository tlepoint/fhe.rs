#![warn(missing_docs, unused_imports)]

//! Internal 193-bit signed integer.

/// Structure holding a 193-bit signed integer.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Int193 {
	lo: u64,
	mi: u64,
	hi: u64,
	sign: bool,
}

impl Int193 {
	/// Create an Int193 holding the value 0.
	pub fn zero() -> Self {
		Self {
			lo: 0,
			mi: 0,
			hi: 0,
			sign: false,
		}
	}

	/// Add an Int193 represented as three u64 and a sign.
	///
	/// The Int193 is defined as `lo + 2^64 * mi + 2^128 * hi`, and `sign` indicates
	/// whether it is negative.
	pub fn add(&mut self, lo: u64, mi: u64, hi: u64, sign: bool) {
		if sign {
			let mut borrow;
			(self.lo, borrow) = self.lo.borrowing_sub(lo, false);
			(self.mi, borrow) = self.mi.borrowing_sub(mi, borrow);
			(self.hi, borrow) = self.hi.borrowing_sub(hi, borrow);
			if borrow {
				self.sign = true
			}
		} else {
			let mut carry;
			(self.lo, carry) = self.lo.carrying_add(lo, false);
			(self.mi, carry) = self.mi.carrying_add(mi, carry);
			(self.hi, carry) = self.hi.carrying_add(hi, carry);
			if carry {
				self.sign = false
			}
		}
	}

	/// Returns true if self is negative and false if the number is zero or positive.
	pub fn is_negative(&self) -> bool {
		self.sign
	}

	/// Returns the least significant bit.
	pub fn lsb(&self) -> u64 {
		self.lo & 1
	}
}

impl std::ops::Shr<usize> for &Int193 {
	type Output = Int193;

	fn shr(self, rhs: usize) -> Self::Output {
		let mut clone = self.clone();
		clone >>= rhs;
		clone
	}
}

impl std::ops::ShrAssign<usize> for Int193 {
	fn shr_assign(&mut self, rhs: usize) {
		debug_assert!(rhs <= 191);

		if rhs > 128 {
			let shift = rhs - 128;
			self.lo = if self.sign {
				(self.hi >> shift) | (u64::MAX << (64 - shift))
			} else {
				self.hi >> shift
			};
			self.mi = if self.sign { u64::MAX } else { 0 };
			self.hi = if self.sign { u64::MAX } else { 0 };
		} else if rhs == 128 {
			self.lo = self.hi;
			self.mi = if self.sign { u64::MAX } else { 0 };
			self.hi = if self.sign { u64::MAX } else { 0 };
		} else if rhs > 64 {
			let shift = rhs - 64;
			self.lo = (self.mi >> shift) | (self.hi << (64 - shift));
			self.mi = if self.sign {
				(self.hi >> shift) | (u64::MAX << (64 - shift))
			} else {
				self.hi >> shift
			};
			self.hi = if self.sign { u64::MAX } else { 0 };
		} else if rhs == 64 {
			self.lo = self.mi;
			self.mi = self.hi;
			self.hi = if self.sign { u64::MAX } else { 0 };
		} else if rhs > 0 {
			self.lo = (self.lo >> rhs) | (self.mi << (64 - rhs));
			self.mi = (self.mi >> rhs) | (self.hi << (64 - rhs));
			self.hi = if self.sign {
				(self.hi >> rhs) | (u64::MAX << (64 - rhs))
			} else {
				self.hi >> rhs
			};
		}
	}
}

impl From<&Int193> for i128 {
	fn from(a: &Int193) -> Self {
		if a.sign {
			debug_assert_eq!(a.hi, u64::MAX);
			debug_assert_eq!(a.mi >> 63, 1);
			-(1i128 + (!a.lo as i128) + ((!a.mi as i128) << 64))
		} else {
			debug_assert_eq!(a.hi, 0);
			debug_assert_eq!(a.mi >> 63, 0);
			(a.lo as i128) + (((a.mi << 1) as i128) << 63)
		}
	}
}

impl From<i128> for Int193 {
	fn from(a: i128) -> Self {
		let mut abs = a.unsigned_abs();
		if a.is_negative() {
			abs -= 1;
			Self {
				lo: !abs as u64,
				mi: !((abs >> 64) as u64),
				hi: u64::MAX,
				sign: a.is_negative(),
			}
		} else {
			Self {
				lo: abs as u64,
				mi: (abs >> 64) as u64,
				hi: 0,
				sign: false,
			}
		}
	}
}

#[cfg(test)]
mod tests {

	use super::Int193;

	#[test]
	fn test_from_into() {
		for a in [0, 1, -1, i128::MIN + 1, i128::MAX] {
			let b = Int193::from(a);
			assert_eq!(b.is_negative(), a.is_negative());

			let c: i128 = (&b).into();
			assert_eq!(c, a);
		}

		let b = Int193::zero();
		assert!(!b.is_negative());

		let c: i128 = (&b).into();
		assert_eq!(c, 0i128);
	}

	#[test]
	fn test_add_positive() {
		for lo in [0u64, 1u64, u64::MAX] {
			for mi in [0u64, 1u64, u64::MAX >> 1] {
				let mut a = Int193::zero();
				a.add(lo, mi, 0, false);
				let c: i128 = (&a).into();

				assert!(!a.is_negative());
				assert_eq!(c, (lo as i128) + ((mi as i128) << 64));
			}
		}
	}

	#[test]
	fn test_add_negative() {
		let mut a: Int193;
		let mut c: i128;

		a = Int193::zero();
		a.add(0, 0, 0, true);
		c = (&a).into();
		assert!(!a.is_negative());
		assert_eq!(c, 0i128);

		a = Int193::zero();
		a.add(1, 0, 0, true);
		c = (&a).into();
		assert!(a.is_negative());
		assert_eq!(c, -1);

		a = Int193::zero();
		a.add(1, 1, 0, true);
		assert!(a.is_negative());
		c = (&a).into();
		assert_eq!(c, -((1i128 << 64) + 1));

		for lo in [0u64, 1u64, u64::MAX] {
			for mi in [0u64, 1u64, u64::MAX >> 1] {
				let mut a = Int193::zero();
				a.add(lo, mi, 0, true);
				let c: i128 = (&a).into();

				assert_eq!(a.is_negative(), lo != 0 || mi != 0);
				assert_eq!(-c, (lo as i128) + ((mi as i128) << 64));
			}
		}
	}

	#[test]
	fn test_shift() {
		let mut a: Int193;

		a = Int193::zero();
		a.add(2, 0, u64::MAX, false);
		a >>= 1;
		assert!((a.lo == 1) && (a.mi == 1u64 << 63) && (a.hi == u64::MAX >> 1));
		a >>= 64;
		assert!((a.lo == 1u64 << 63) && (a.mi == u64::MAX >> 1) && (a.hi == 0));

		a = Int193::zero();
		a.add(2, 0, u64::MAX, false);
		a >>= 65;
		assert!((a.lo == 1u64 << 63) && (a.mi == u64::MAX >> 1) && (a.hi == 0));

		a = Int193::zero();
		a.add(2, u64::MAX - 2, u64::MAX - 5, false);
		a >>= 128;
		assert!((a.lo == u64::MAX - 5) && (a.mi == 0) && (a.hi == 0));
	}
}
