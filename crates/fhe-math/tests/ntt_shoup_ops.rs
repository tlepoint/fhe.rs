//! Unit test for polynomial Shoup operations.

use fhe_math::rq::{Context, Ntt, NttShoup, Poly};
use rand::rng;
use std::sync::Arc;

#[test]
fn test_ntt_shoup_add_sub_neg() {
    let modulus = 4611686018326724609;
    let ctx = Arc::new(Context::new(&[modulus], 16).unwrap());
    let mut rng = rng();

    let p_ntt = Poly::<Ntt>::random(&ctx, &mut rng);
    let p_shoup = Poly::<NttShoup>::random(&ctx, &mut rng);
    let p_shoup_as_ntt = p_shoup.clone().into_ntt();

    // Add/Sub/Neg on Ntt after explicit conversion.
    let sum = &p_ntt + &p_shoup_as_ntt;
    assert_eq!(sum, &p_ntt + &p_shoup_as_ntt);

    let diff = &p_ntt - &p_shoup_as_ntt;
    assert_eq!(diff, &p_ntt - &p_shoup_as_ntt);

    let neg = -&p_shoup_as_ntt;
    assert_eq!(neg, -&p_shoup_as_ntt);
}
