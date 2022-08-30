use fhers_traits::{FheParametrized, FhePlaintext};
use zeroize::{Zeroize, ZeroizeOnDrop};

use crate::bfv::{BfvParameters, Encoding, Plaintext};

/// A wrapper around a vector of plaintext which implements the [`FhePlaintext`]
/// trait, and therefore can be encoded to / decoded from.
pub struct PlaintextVec(pub Vec<Plaintext>);

impl FhePlaintext for PlaintextVec {
	type Encoding = Encoding;
}

impl FheParametrized for PlaintextVec {
	type Parameters = BfvParameters;
}

impl Zeroize for PlaintextVec {
	fn zeroize(&mut self) {
		self.0.zeroize()
	}
}

impl ZeroizeOnDrop for PlaintextVec {}

#[cfg(test)]
mod tests {
	use std::sync::Arc;

	use fhers_traits::{FheDecoder, FheEncoder};

	use crate::bfv::{BfvParameters, Encoding, PlaintextVec};

	#[test]
	fn encode_decode() {
		(0..40).for_each(|_| {
			for i in 1..5 {
				let params = Arc::new(BfvParameters::default(1, 8));
				let a = params.plaintext.random_vec(params.degree() * i);

				let plaintexts =
					PlaintextVec::try_encode(&a as &[u64], Encoding::poly_at_level(0), &params);
				assert!(plaintexts.is_ok());
				let plaintexts = plaintexts.unwrap();
				assert_eq!(plaintexts.0.len(), i);
				for j in 0..i {
					let b = Vec::<u64>::try_decode(&plaintexts.0[j], Encoding::poly_at_level(0));
					assert!(
						b.is_ok_and(|b| b == &a[j * params.degree()..(j + 1) * params.degree()])
					);
				}
			}
		})
	}
}
