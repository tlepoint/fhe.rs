//! Ciphertext type in the BFV encryption scheme.

use crate::{parameters::BfvParameters, RelinearizationKey};
use itertools::izip;
use math::rq::{
	traits::{ContextSwitcher, TryConvertFrom},
	Poly, Representation,
};
use ndarray::{Array2, Axis};
use num_bigint::BigUint;
use rand::SeedableRng;
use rand_chacha::ChaCha8Rng;
use std::{
	ops::{Add, AddAssign, Neg, Sub, SubAssign},
	rc::Rc,
};

/// A ciphertext encrypting a plaintext.
#[derive(Debug, Clone, PartialEq)]
pub struct Ciphertext {
	/// The parameters of the underlying BFV encryption scheme.
	pub(crate) par: Rc<BfvParameters>,

	/// The seed that generated the polynomial c0 in a fresh ciphertext.
	pub(crate) seed: Option<<ChaCha8Rng as SeedableRng>::Seed>,

	/// The ciphertext element c0.
	pub(crate) c0: Poly,

	/// The ciphertext element c1.
	pub(crate) c1: Poly,
}

impl Add<&Ciphertext> for &Ciphertext {
	type Output = Ciphertext;

	fn add(self, rhs: &Ciphertext) -> Ciphertext {
		debug_assert_eq!(self.par, rhs.par);

		Ciphertext {
			par: self.par.clone(),
			seed: None,
			c0: &self.c0 + &rhs.c0,
			c1: &self.c1 + &rhs.c1,
		}
	}
}

impl AddAssign<&Ciphertext> for Ciphertext {
	fn add_assign(&mut self, rhs: &Ciphertext) {
		debug_assert_eq!(self.par, rhs.par);

		self.c0 += &rhs.c0;
		self.c1 += &rhs.c1;
		self.seed = None
	}
}

impl Sub<&Ciphertext> for &Ciphertext {
	type Output = Ciphertext;

	fn sub(self, rhs: &Ciphertext) -> Ciphertext {
		assert_eq!(self.par, rhs.par);

		Ciphertext {
			par: self.par.clone(),
			seed: None,
			c0: &self.c0 - &rhs.c0,
			c1: &self.c1 - &rhs.c1,
		}
	}
}

impl SubAssign<&Ciphertext> for Ciphertext {
	fn sub_assign(&mut self, rhs: &Ciphertext) {
		debug_assert_eq!(self.par, rhs.par);

		self.c0 -= &rhs.c0;
		self.c1 -= &rhs.c1;
		self.seed = None
	}
}

impl Neg for &Ciphertext {
	type Output = Ciphertext;

	fn neg(self) -> Ciphertext {
		Ciphertext {
			par: self.par.clone(),
			seed: None,
			c0: -&self.c0,
			c1: -&self.c1,
		}
	}
}

#[allow(dead_code)]
fn print_poly(s: &str, p: &Poly) {
	println!("{} = {:?}", s, Vec::<BigUint>::from(p))
}

/// Multiply two ciphertext and relinearize.
pub fn mul(
	ct0: &Ciphertext,
	ct1: &Ciphertext,
	rk: &RelinearizationKey,
) -> Result<Ciphertext, String> {
	// Extend
	let mut now = std::time::SystemTime::now();
	let c00 = ct0.par.extender.switch_context(&ct0.c0).unwrap();
	let c01 = ct0.par.extender.switch_context(&ct0.c1).unwrap();
	let c10 = ct1.par.extender.switch_context(&ct1.c0).unwrap();
	let c11 = ct1.par.extender.switch_context(&ct1.c1).unwrap();
	println!("Extend: {:?}", now.elapsed().unwrap());

	// Multiply
	now = std::time::SystemTime::now();
	let mut c0 = &c00 * &c10;
	let mut c1 = &c00 * &c11;
	c1 += &(&c01 * &c10);
	let mut c2 = &c01 * &c11;
	c0.change_representation(Representation::PowerBasis);
	c1.change_representation(Representation::PowerBasis);
	c2.change_representation(Representation::PowerBasis);
	println!("Multiply: {:?}", now.elapsed().unwrap());

	// Scale
	// TODO: This should be faster??
	now = std::time::SystemTime::now();
	let mut c0_scaled_coeffs = Array2::zeros((ct0.par.moduli().len(), ct0.par.degree()));
	let mut c1_scaled_coeffs = Array2::zeros((ct0.par.moduli().len(), ct0.par.degree()));
	let mut c2_scaled_coeffs = Array2::zeros((ct0.par.moduli().len(), ct0.par.degree()));

	let scaler = &ct0.par.extended_scaler;

	izip!(
		c0_scaled_coeffs.axis_iter_mut(Axis(1)),
		c0.coefficients().axis_iter(Axis(1)),
		c1_scaled_coeffs.axis_iter_mut(Axis(1)),
		c1.coefficients().axis_iter(Axis(1)),
		c2_scaled_coeffs.axis_iter_mut(Axis(1)),
		c2.coefficients().axis_iter(Axis(1))
	)
	.for_each(
		|(
			mut c0_scaled_column,
			c0_column,
			mut c1_scaled_column,
			c1_column,
			mut c2_scaled_column,
			c2_column,
		)| {
			scaler.scale(&c0_column, &mut c0_scaled_column, false);
			scaler.scale(&c1_column, &mut c1_scaled_column, false);
			scaler.scale(&c2_column, &mut c2_scaled_column, false);
		},
	);
	println!("Scale: {:?}", now.elapsed().unwrap());

	// Relinearize
	now = std::time::SystemTime::now();
	let mut c0 =
		Poly::try_convert_from(c0_scaled_coeffs, ct0.par.ctx(), Representation::PowerBasis)?;
	let mut c1 =
		Poly::try_convert_from(c1_scaled_coeffs, ct0.par.ctx(), Representation::PowerBasis)?;
	let c2 = Poly::try_convert_from(c2_scaled_coeffs, ct0.par.ctx(), Representation::PowerBasis)?;

	c0.change_representation(Representation::Ntt);
	c1.change_representation(Representation::Ntt);
	let out = rk.relinearize(&c0, &c1, &c2);
	println!("Relinearize: {:?}", now.elapsed().unwrap());

	out
}

#[cfg(test)]
mod tests {
	use super::mul;
	use crate::{
		traits::{Decoder, Decryptor, Encoder, Encryptor},
		BfvParameters, Encoding, Plaintext, RelinearizationKey, SecretKey,
	};
	use proptest::collection::vec as prop_vec;
	use proptest::prelude::any;
	use std::rc::Rc;

	proptest! {
		#[test]
		fn test_add(mut a in prop_vec(any::<u64>(), 8), mut b in prop_vec(any::<u64>(), 8)) {
			for params in [Rc::new(BfvParameters::default_one_modulus()),
						   Rc::new(BfvParameters::default_two_moduli())] {
				params.plaintext().reduce_vec(&mut a);
				params.plaintext().reduce_vec(&mut b);
				let mut c = a.clone();
				params.plaintext().add_vec(&mut c, &b);

				let sk = SecretKey::random(&params);

				for encoding in [Encoding::Poly, Encoding::Simd] {
					let pt_a = Plaintext::try_encode(&a, encoding.clone(), &params).unwrap();
					let pt_b = Plaintext::try_encode(&b, encoding.clone(), &params).unwrap();

					let mut ct_a = sk.encrypt(&pt_a).unwrap();
					let ct_b = sk.encrypt(&pt_b).unwrap();
					let ct_c = &ct_a + &ct_b;
					ct_a += &ct_b;

					let pt_c = sk.decrypt(&ct_c).unwrap();
					assert_eq!(Vec::<u64>::try_decode(&pt_c, encoding.clone()).unwrap(), c);
					let pt_c = sk.decrypt(&ct_a).unwrap();
					assert_eq!(Vec::<u64>::try_decode(&pt_c, encoding.clone()).unwrap(), c);
				}
			}
		}

		#[test]
		fn test_sub(mut a in prop_vec(any::<u64>(), 8), mut b in prop_vec(any::<u64>(), 8)) {
			for params in [Rc::new(BfvParameters::default_one_modulus()),
						   Rc::new(BfvParameters::default_two_moduli())] {
				params.plaintext().reduce_vec(&mut a);
				params.plaintext().reduce_vec(&mut b);
				let mut c = a.clone();
				params.plaintext().sub_vec(&mut c, &b);

				let sk = SecretKey::random(&params);

				for encoding in [Encoding::Poly, Encoding::Simd] {
					let pt_a = Plaintext::try_encode(&a, encoding.clone(), &params).unwrap();
					let pt_b = Plaintext::try_encode(&b, encoding.clone(), &params).unwrap();

					let mut ct_a = sk.encrypt(&pt_a).unwrap();
					let ct_b = sk.encrypt(&pt_b).unwrap();
					let ct_c = &ct_a - &ct_b;
					ct_a -= &ct_b;

					let pt_c = sk.decrypt(&ct_c).unwrap();
					assert_eq!(Vec::<u64>::try_decode(&pt_c, encoding.clone()).unwrap(), c);
					let pt_c = sk.decrypt(&ct_a).unwrap();
					assert_eq!(Vec::<u64>::try_decode(&pt_c, encoding.clone()).unwrap(), c);
				}
			}
		}

		#[test]
		fn test_neg(mut a in prop_vec(any::<u64>(), 8)) {
			for params in [Rc::new(BfvParameters::default_one_modulus()),
						   Rc::new(BfvParameters::default_two_moduli())] {
				params.plaintext().reduce_vec(&mut a);
				let mut c = a.clone();
				params.plaintext().neg_vec(&mut c);

				let sk = SecretKey::random(&params);
				for encoding in [Encoding::Poly, Encoding::Simd] {
					let pt_a = Plaintext::try_encode(&a, encoding.clone(), &params).unwrap();

					let ct_a = sk.encrypt(&pt_a).unwrap();
					let ct_c = -&ct_a;

					let pt_c = sk.decrypt(&ct_c).unwrap();
					assert_eq!(Vec::<u64>::try_decode(&pt_c, encoding.clone()).unwrap(), c);
				}
			}
		}
	}

	#[test]
	fn test_mul() {
		let par = Rc::new(BfvParameters::default_two_moduli());

		// We will encode `values` in an Simd format, and check that the product is computed correctly.
		let values = par.plaintext().random_vec(par.degree());
		let mut expected = values.clone();
		par.plaintext().mul_vec(&mut expected, &values);

		let sk = SecretKey::random(&par);
		let rk = RelinearizationKey::new(&sk).unwrap();
		let pt = Plaintext::try_encode(&values, Encoding::Simd, &par).unwrap();

		let ct1 = sk.encrypt(&pt).unwrap();
		let ct2 = sk.encrypt(&pt).unwrap();
		let ct3 = mul(&ct1, &ct2, &rk).unwrap();

		println!("Noise: {}", unsafe { sk.measure_noise(&ct3).unwrap() });
		let pt = sk.decrypt(&ct3).unwrap();
		assert_eq!(
			Vec::<u64>::try_decode(&pt, Encoding::Simd).unwrap(),
			expected
		);
	}
}

/*
q = [4611686018326724609, 4611686018309947393]
qp = [4611686018326724609, 4611686018309947393, 4611686018427387761, 4611686018427387617, 4611686018427387409]
Q = prod(q)
QP = prod(qp)
print(Q)
print(QP)
PP_QP.<XX> = PolynomialRing(Integers(QP))
R_QP.<X> = QuotientRing(PP_QP, XX^8+1)
PP_Q.<YY> = PolynomialRing(Integers(Q))
R_Q.<Y> = QuotientRing(PP_Q, YY^8+1)

def to_poly_QP(l, m):
	out = R_QP(0)
	for i in range(len(l)):
		if l[i] > (m >> 1):
			out += -R_QP(m-l[i]) * X^i
		else:
			out += R_QP(l[i]) * X^i
	return out
def to_poly_Q(l, m):
	out = R_Q(0)
	for i in range(len(l)):
		if l[i] > (m >> 1):
			out += -R_Q(m-l[i]) * Y^i
		else:
			out += R_Q(l[i]) * Y^i
	return out
m = [0, 0, 0, 0, 0, 0, 0, 0]
sk = [21267647931552827693776735476788494336, 21267647931552827693776735476788494336, 0, 0, 21267647931552827693776735476788494336, 1, 21267647931552827693776735476788494335, 21267647931552827693776735476788494336]
c00 = [17979034361998346255387772017037258607, 15796959577912698470831321809245575928, 19706298130934675706773674968457836373, 10794789768192925399934612065243809553, 10012727161876898756861647704265740002, 5284862896408540186353639085482949026, 14381082865694479880889794468306241614, 7486565450830229994418719784817455]
c01 = [3413170594790757764013135880673027487, 3637860477249825409694551869226806935, 11535330947900608776900086259459784662, 9039163981771718721030463516053331743, 15178198706057896137644415412341900430, 19353907151847390009113239015006699449, 17960534915281292880168378788075310746, 16389522930113674209846894009566306946]
c10 = [2691438789595652636719926324755836744, 8079062388787788416697448804292608149, 1800704143832233364185219383331763859, 19135365993747922589140962396833614682, 860691537301095266062501124003352508, 12473780818733416438750986758948770324, 1896097118849262389699419508619436387, 15650904598703575991783267165935128093]
c11 = [1316867955315271715782664797147936227, 19662303104517677336891774442339291682, 15662576057340609342231966321842203268, 1101549506482930334087417275893830153, 14580710721590246887703805376151447800, 9615144322328937513273994937716203987, 12876348390691566718717488290369209325, 19229404457624717279197585193264760965]
c0 = [2085924839667862180058226974748526733911397149918662699454371795143060706772543866027381743557, 38599293999824320220378037448738562411427152479335275074724018041111685534, 2085924839667862180092237815701085074455044268532214489028832797565162868373922197618218344671, 2085924839667862180113963035866082320089476094892795331630061633550649902264835939211397196095, 2085924839667862180163118808647573200523166825013012103559707596376709753030283478432731957032, 65104843177539539964535601032988537761454376876183506295805888099069235827, 105094199541069449575649547171711076540827599677013669603504961253895579315, 2085924839667862180117315273594990651773546960183000197908163565194025552126601846871921156235]
c1 = [268748494487531508157961182818553713256733034061925255716589824008507292012, 2085924839667862179941049253151636969085119985918606110010439956352258781612889887484051984653, 38234940861292873408882450393579489835096641753090649532347939659935995613, 36292029304116212545823905900241766578917845784625457664450455585533280665, 2085924839667862180166495801348263288727294095330942090584535381128315801883661561320830953162, 2085924839667862180157277086194903779398095006531909412979705568894551432959467193557545875453, 2085924839667862179956293567451983110182009874144067487245833361176815033565609498850609192301, 150882831992649484257652154449339230638952023273845718022728277377606824936]
c2 = [2085924839667862179950301375138779215942625177117230714666243110992040308089023006359411728170, 77379343305236200640718758558999887951197321460170050123573430030411430873, 2085924839667862180103149145377515818020811500590650448679655164480960553667989050650594364492, 2085924839667862180153555369486864086573438933041114001891450718501881374679251094658448408393, 28456216786823645283296879589218781455812238347805645679746415544985682046, 2085924839667862180078638054148823105189245296637446946976065553069145588255506498498352792318, 104280272707906811579044796117078484275491502392487473464291237716209497769, 2085924839667862179983349682444390532813775045648654034163797747998108160899270404859344181570]
c0p = [2214460209890395463727984965976617054, 8384887007909973531981371521512544877, 17053649048915092429236503645435701451, 3872094185515102171718337949668277628, 10337310740777919526968826656358002124, 20419003357633764889179799446549059102, 19093860865619131061260697155700404802, 15468454718714558785832338273058183964]
c1p = [1538152878873629169082294574666325794, 8591691619169993526245420662609165533, 9899566953590880788542899966899135222, 10905189077229487514080724357532826297, 2008084418300242521750712917802490394, 12650060149537536680063130231715385800, 5605620710711523362019543534655640558, 13154790265462723415946582529332142133]
c2p = [21028475206859524568828868495285318059, 5301557973171821646617013644232118731, 13104145912544932545075710812396797553, 2290488267658553355020637182045500119, 11448979270134516451679631842212261195, 2858933189728052268215444548041326094, 17503090618507706481762550052046110517, 4952705865954005767464953062310792713]

t = 1153

print(R_QP, X^9)
print((to_poly_QP(c00, Q) * to_poly_QP(c10, Q)).list() == [Integers(QP)(c) for c in c0])
print((to_poly_QP(c00, Q) * to_poly_QP(c11, Q)+to_poly_QP(c01, Q) * to_poly_QP(c10, Q)).list() == [Integers(QP)(c) for c in c1])
print((to_poly_QP(c01, Q) * to_poly_QP(c11, Q)).list() == [Integers(QP)(c) for c in c2])

def center(a, b):
	if a > b >> 1:
		return a - b
	else:
		return a
print(prod(qp[2:]))
print([(round(center(c0i, QP) * t / Q) % Q) - qii for c0i, qii in zip(c0, c0p)])
print([(round(center(c1i, QP) * t / Q) % Q) - qii for c1i, qii in zip(c1, c1p)])
#for c1i, qii in zip(c1, c1p):
#    print((round(center(c1i, Q) * t / Q) % Q) - qii, (round(center(c1i, Q) * t / Q) % Q), qii)
print([(round(center(c2i, QP) * t / Q) % Q) - qii for c2i, qii in zip(c2, c2p)])

delta = Q // t
print(delta)
print("ct0 dec:", (to_poly_Q(c00, Q) + to_poly_Q(c01, Q) * to_poly_Q(sk, Q) - (to_poly_Q(m, Q) * delta)).list())
print("ct1 dec:", (to_poly_Q(c10, Q) + to_poly_Q(c11, Q) * to_poly_Q(sk, Q) - (to_poly_Q(m, Q) * delta)).list())

m2 = [(mi * mi) % t for mi in m]
print(delta * delta)
print((to_poly_Q(c0p, Q) + to_poly_Q(c1p, Q) * to_poly_Q(sk, Q) + to_poly_Q(c2p, Q) * to_poly_Q(sk, Q) * to_poly_Q(sk, Q)).list())

print((to_poly_QP(c0, Q) + to_poly_QP(c1, Q) * to_poly_QP(sk, Q) + to_poly_QP(c2, Q) * to_poly_QP(sk, Q) * to_poly_QP(sk, Q)).list())

pc00 = to_poly_QP(c00, Q)
pc01 = to_poly_QP(c01, Q)
pc10 = to_poly_QP(c10, Q)
pc11 = to_poly_QP(c11, Q)
ps = to_poly_QP(sk, Q)
print()
print()
print((pc00 + pc01 * ps).list())
print((pc10 + pc11 * ps).list())
print(((pc00 + pc01 * ps) * (pc10 + pc11 * ps)).list())

for coeff in ((pc00 + pc01 * ps) * (pc10 + pc11 * ps)).list():
	coeff = coeff.lift()
	print(coeff, center(round(center(coeff, QP) * t / Q) % Q, Q))



*/
