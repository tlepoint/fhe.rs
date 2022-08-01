//! Ciphertext type in the BFV encryption scheme.

use crate::parameters::BfvParameters;
use math::rq::Poly;
use rand::SeedableRng;
use rand_chacha::ChaCha8Rng;
use std::rc::Rc;

/// A ciphertext encrypting a plaintext.
/// TODO: Can the members be private?
#[derive(Debug, Clone, PartialEq)]
pub struct Ciphertext {
	/// The parameters of the underlying BFV encryption scheme.
	pub par: Rc<BfvParameters>,

	/// The seed that generated the polynomial c0 in a fresh ciphertext.
	pub seed: Option<<ChaCha8Rng as SeedableRng>::Seed>,

	/// The ciphertext element c0.
	pub c0: Poly,

	/// The ciphertext element c1.
	pub c1: Poly,
}
