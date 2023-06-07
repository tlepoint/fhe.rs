//! Secret keys for the BFV encryption scheme
use crate::bfv::{BfvParameters, Ciphertext, Plaintext};
use crate::{Error, Result};

use fhe_math::{
    rq::{traits::TryConvertFrom, Poly, Representation},
    zq::Modulus,
};
use fhe_traits::{Deserialize, FheDecrypter, FheEncrypter, FheParametrized, Serialize};
use fhe_util::sample_vec_cbd;
use itertools::Itertools;
use num_bigint::BigUint;
use rand::{thread_rng, CryptoRng, Rng, RngCore, SeedableRng};
use rand_chacha::ChaCha8Rng;
use std::fs::File;
use std::io::{Read, Write};
use std::sync::Arc;
use zeroize::{Zeroize, ZeroizeOnDrop, Zeroizing};

/// Secret key for the BFV encryption scheme.
#[derive(Debug, PartialEq, Eq, Clone)]
pub struct SecretKey {
    pub(crate) par: Arc<BfvParameters>,
    pub(crate) coeffs: Box<[i64]>,
}

impl Zeroize for SecretKey {
    fn zeroize(&mut self) {
        self.coeffs.zeroize();
    }
}

impl ZeroizeOnDrop for SecretKey {}

impl SecretKey {
    /// Generate a random [`SecretKey`] and write it to a file.
    pub fn random_write_key<R: RngCore + CryptoRng>(
        par: &Arc<BfvParameters>,
        rng: &mut R,
        file_name: &mut String,
    ) -> Self {
        let s_coefficients = sample_vec_cbd(par.degree(), par.variance, rng).unwrap();
        let secret_key = Self::new(s_coefficients, par);

        // write the secret key par to a file

        let par_file = file_name.clone() + ".par";

        let mut file = File::create(par_file).unwrap();

        let binding = secret_key.par.to_bytes();
        let par_bytes = binding.as_slice();
        file.write_all(par_bytes).unwrap();

        // write the coefficients to a file

        let coeff_file = file_name.clone() + ".coeff";

        let mut file = File::create(coeff_file).unwrap();

        let binding = secret_key.coeffs.to_vec();
        let coeff_bytes = binding.as_slice();

        let coeff_bytes_as_u8: &[u8] = unsafe {
            // Convert the &[i64] to a byte slice by transmuting the reference type
            std::slice::from_raw_parts(
                coeff_bytes.as_ptr() as *const u8,
                coeff_bytes.len() * std::mem::size_of::<i64>(),
            )
        };

        file.write_all(coeff_bytes_as_u8).unwrap();

        secret_key
    }

    /// Read a [`SecretKey`] from a file.
    pub fn read_key(file_name: String) -> Self {
        // read the secret key par from a file
        let par_file = file_name.clone() + ".par";

        let mut file = File::open(par_file).unwrap();
        let mut par_bytes = Vec::new();
        file.read_to_end(&mut par_bytes).unwrap();

        let par = Arc::new(BfvParameters::try_deserialize(&par_bytes).unwrap());

        let coeff_file = file_name.clone() + ".coeff";

        let mut file = File::open(coeff_file).unwrap();
        let mut coeff_bytes = Vec::new();
        file.read_to_end(&mut coeff_bytes).unwrap();

        let coeff_bytes_as_i64: &[i64] = unsafe {
            // Convert the &[u8] to a byte slice by transmuting the reference type
            std::slice::from_raw_parts(
                coeff_bytes.as_ptr() as *const i64,
                coeff_bytes.len() / std::mem::size_of::<i64>(),
            )
        };

        let coeffs = coeff_bytes_as_i64.to_vec();

        let secret_key = Self::new(coeffs, &par);

        secret_key
    }
    /// Generate a random [`SecretKey`].
    pub fn random<R: RngCore + CryptoRng>(par: &Arc<BfvParameters>, rng: &mut R) -> Self {
        let s_coefficients = sample_vec_cbd(par.degree(), par.variance, rng).unwrap();
        let secret_key = Self::new(s_coefficients, par);

        secret_key
    }

    /// Generate a [`SecretKey`] from its coefficients.
    pub(crate) fn new(coeffs: Vec<i64>, par: &Arc<BfvParameters>) -> Self {
        Self {
            par: par.clone(),
            coeffs: coeffs.into_boxed_slice(),
        }
    }

    /// Measure the noise in a [`Ciphertext`].
    ///
    /// # Safety
    ///
    /// This operations may run in a variable time depending on the value of the
    /// noise.
    pub unsafe fn measure_noise(&self, ct: &Ciphertext) -> Result<usize> {
        let plaintext = Zeroizing::new(self.try_decrypt(ct)?);
        let m = Zeroizing::new(plaintext.to_poly());

        // Let's create a secret key with the ciphertext context
        let mut s = Zeroizing::new(Poly::try_convert_from(
            self.coeffs.as_ref(),
            ct.c[0].ctx(),
            false,
            Representation::PowerBasis,
        )?);
        s.change_representation(Representation::Ntt);
        let mut si = s.clone();

        // Let's disable variable time computations
        let mut c = Zeroizing::new(ct.c[0].clone());
        c.disallow_variable_time_computations();

        for i in 1..ct.c.len() {
            let mut cis = Zeroizing::new(ct.c[i].clone());
            cis.disallow_variable_time_computations();
            *cis.as_mut() *= si.as_ref();
            *c.as_mut() += &cis;
            *si.as_mut() *= s.as_ref();
        }
        *c.as_mut() -= &m;
        c.change_representation(Representation::PowerBasis);

        let ciphertext_modulus = ct.c[0].ctx().modulus();
        let mut noise = 0usize;
        for coeff in Vec::<BigUint>::from(c.as_ref()) {
            noise = std::cmp::max(
                noise,
                std::cmp::min(coeff.bits(), (ciphertext_modulus - &coeff).bits()) as usize,
            )
        }

        Ok(noise)
    }

    pub(crate) fn encrypt_poly<R: RngCore + CryptoRng>(
        &self,
        p: &Poly,
        rng: &mut R,
    ) -> Result<Ciphertext> {
        assert_eq!(p.representation(), &Representation::Ntt);

        let level = self.par.level_of_ctx(p.ctx())?;

        let mut seed = <ChaCha8Rng as SeedableRng>::Seed::default();
        thread_rng().fill(&mut seed);

        // Let's create a secret key with the ciphertext context
        let mut s = Zeroizing::new(Poly::try_convert_from(
            self.coeffs.as_ref(),
            p.ctx(),
            false,
            Representation::PowerBasis,
        )?);
        s.change_representation(Representation::Ntt);

        let mut a = Poly::random_from_seed(p.ctx(), Representation::Ntt, seed);
        let a_s = Zeroizing::new(&a * s.as_ref());

        let mut b = Poly::small(p.ctx(), Representation::Ntt, self.par.variance, rng)
            .map_err(Error::MathError)?;
        b -= &a_s;
        b += p;

        // It is now safe to enable variable time computations.
        unsafe {
            a.allow_variable_time_computations();
            b.allow_variable_time_computations()
        }

        Ok(Ciphertext {
            par: self.par.clone(),
            seed: Some(seed),
            c: vec![b, a],
            level,
        })
    }
}

impl FheParametrized for SecretKey {
    type Parameters = BfvParameters;
}

impl FheEncrypter<Plaintext, Ciphertext> for SecretKey {
    type Error = Error;

    fn try_encrypt<R: RngCore + CryptoRng>(
        &self,
        pt: &Plaintext,
        rng: &mut R,
    ) -> Result<Ciphertext> {
        assert_eq!(self.par, pt.par);
        let m = Zeroizing::new(pt.to_poly());
        self.encrypt_poly(m.as_ref(), rng)
    }
}

impl FheDecrypter<Plaintext, Ciphertext> for SecretKey {
    type Error = Error;

    fn try_decrypt(&self, ct: &Ciphertext) -> Result<Plaintext> {
        if self.par != ct.par {
            Err(Error::DefaultError(
                "Incompatible BFV parameters".to_string(),
            ))
        } else {
            // Let's create a secret key with the ciphertext context
            let mut s = Zeroizing::new(Poly::try_convert_from(
                self.coeffs.as_ref(),
                ct.c[0].ctx(),
                false,
                Representation::PowerBasis,
            )?);
            s.change_representation(Representation::Ntt);
            let mut si = s.clone();

            let mut c = Zeroizing::new(ct.c[0].clone());
            c.disallow_variable_time_computations();

            for i in 1..ct.c.len() {
                let mut cis = Zeroizing::new(ct.c[i].clone());
                cis.disallow_variable_time_computations();
                *cis.as_mut() *= si.as_ref();
                *c.as_mut() += &cis;
                *si.as_mut() *= s.as_ref();
            }
            c.change_representation(Representation::PowerBasis);

            let d = Zeroizing::new(c.scale(&self.par.scalers[ct.level])?);

            // TODO: Can we handle plaintext moduli that are BigUint?
            let v = Zeroizing::new(
                Vec::<u64>::from(d.as_ref())
                    .iter_mut()
                    .map(|vi| *vi + self.par.plaintext.modulus())
                    .collect_vec(),
            );
            let mut w = v[..self.par.degree()].to_vec();
            let q = Modulus::new(self.par.moduli[0]).map_err(Error::MathError)?;
            q.reduce_vec(&mut w);
            self.par.plaintext.reduce_vec(&mut w);

            let mut poly =
                Poly::try_convert_from(&w, ct.c[0].ctx(), false, Representation::PowerBasis)?;
            poly.change_representation(Representation::Ntt);

            let pt = Plaintext {
                par: self.par.clone(),
                value: w.into_boxed_slice(),
                encoding: None,
                poly_ntt: poly,
                level: ct.level,
            };

            Ok(pt)
        }
    }
}

#[cfg(test)]
mod tests {
    use super::SecretKey;
    use crate::bfv::{parameters::BfvParameters, Encoding, Plaintext};
    use fhe_traits::{FheDecrypter, FheEncoder, FheEncrypter};
    use rand::thread_rng;
    use std::{error::Error, sync::Arc};

    #[test]
    fn keygen() {
        let mut rng = thread_rng();
        let params = Arc::new(BfvParameters::default(1, 8));
        let sk = SecretKey::random(&params, &mut rng);
        assert_eq!(sk.par, params);

        sk.coeffs.iter().for_each(|ci| {
            // Check that this is a small polynomial
            assert!((*ci).abs() <= 2 * sk.par.variance as i64)
        })
    }

    #[test]
    fn keygen_and_write() {
        let mut rng = thread_rng();
        let params = Arc::new(BfvParameters::default(1, 8));
        let sk: SecretKey = SecretKey::random_write_key(&params, &mut rng, &mut "key".to_string());
        assert_eq!(sk.par, params);
    }

    #[test]
    fn keygen_write_and_read() {
        let mut rng = thread_rng();
        let params = Arc::new(BfvParameters::default(1, 8));
        let sk: SecretKey = SecretKey::random_write_key(&params, &mut rng, &mut "key".to_string());

        let sk_read = SecretKey::read_key("key".to_string());

        assert_eq!(sk, sk_read);
    }

    #[test]
    fn encrypt_decrypt() -> Result<(), Box<dyn Error>> {
        let mut rng = thread_rng();
        for params in [
            Arc::new(BfvParameters::default(1, 8)),
            Arc::new(BfvParameters::default(6, 8)),
        ] {
            for level in 0..params.max_level() {
                for _ in 0..20 {
                    let sk = SecretKey::random(&params, &mut rng);

                    let pt = Plaintext::try_encode(
                        &params.plaintext.random_vec(params.degree(), &mut rng),
                        Encoding::poly_at_level(level),
                        &params,
                    )?;
                    let ct = sk.try_encrypt(&pt, &mut rng)?;
                    let pt2 = sk.try_decrypt(&ct)?;

                    println!("Noise: {}", unsafe { sk.measure_noise(&ct)? });
                    assert_eq!(pt2, pt);
                }
            }
        }

        Ok(())
    }
}
