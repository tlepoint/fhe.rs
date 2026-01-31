//! Unit test for polynomial Shoup operations.

use fhe_math::rq::{Context, Poly, Representation};
use rand::rng;
use std::sync::Arc;

#[test]
fn test_ntt_shoup_add_sub_neg() {
    let modulus = 4611686018326724609;
    let ctx = Arc::new(Context::new(&[modulus], 16).unwrap());
    let mut rng = rng();

    let p_ntt = Poly::random(&ctx, Representation::Ntt, &mut rng);
    let p_shoup = Poly::random(&ctx, Representation::NttShoup, &mut rng);

    // Helper to get Ntt version of a poly
    let to_ntt = |p: &Poly| -> Poly {
        let mut q = p.clone();
        if *q.representation() == Representation::NttShoup {
            unsafe { q.override_representation(Representation::Ntt) };
            // Note: override_representation handles shoup cleanup if needed or just switch
            // enum But strict conversion:
            q.change_representation(Representation::Ntt);
        }
        q
    };

    let p_shoup_as_ntt = to_ntt(&p_shoup);

    // Case 1: Ntt + NttShoup
    let sum1 = &p_ntt + &p_shoup;
    assert_eq!(sum1.representation(), &Representation::Ntt);
    assert_eq!(sum1, &p_ntt + &p_shoup_as_ntt);

    // Case 2: NttShoup + Ntt
    let sum2 = &p_shoup + &p_ntt;
    assert_eq!(sum2.representation(), &Representation::Ntt);
    assert_eq!(sum2, &p_shoup_as_ntt + &p_ntt);

    // Case 3: NttShoup + NttShoup (should work if we relaxed AddAssign correctly)
    // Wait, AddAssign on LHS=NttShoup is forbidden.
    // But Add(&NttShoup, &NttShoup) -> converts LHS to Ntt, then adds RHS
    // (NttShoup). So LHS becomes Ntt. Ntt += NttShoup. This should work now.
    let p_shoup2 = Poly::random(&ctx, Representation::NttShoup, &mut rng);
    let p_shoup2_as_ntt = to_ntt(&p_shoup2);

    let sum3 = &p_shoup + &p_shoup2;
    assert_eq!(sum3.representation(), &Representation::Ntt);
    assert_eq!(sum3, &p_shoup_as_ntt + &p_shoup2_as_ntt);

    // Case 4: Neg NttShoup
    let neg = -&p_shoup;
    assert_eq!(neg.representation(), &Representation::Ntt);
    assert_eq!(neg, -&p_shoup_as_ntt);

    // Case 5: Sub Ntt - NttShoup
    let diff1 = &p_ntt - &p_shoup;
    assert_eq!(diff1.representation(), &Representation::Ntt);
    assert_eq!(diff1, &p_ntt - &p_shoup_as_ntt);

    // Case 6: Sub NttShoup - Ntt
    let diff2 = &p_shoup - &p_ntt;
    assert_eq!(diff2.representation(), &Representation::Ntt);
    assert_eq!(diff2, &p_shoup_as_ntt - &p_ntt);

    // Case 7: Sub NttShoup - NttShoup
    let diff3 = &p_shoup - &p_shoup2;
    assert_eq!(diff3.representation(), &Representation::Ntt);
    assert_eq!(diff3, &p_shoup_as_ntt - &p_shoup2_as_ntt);
}
