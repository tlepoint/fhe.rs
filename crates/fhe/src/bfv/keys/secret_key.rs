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

#[derive(Debug, serde::Serialize, serde::Deserialize)]
struct Data {
    par: Vec<u8>,
    coeffs: Vec<u8>,
}

impl Zeroize for SecretKey {
    fn zeroize(&mut self) {
        self.coeffs.zeroize();
    }
}

impl ZeroizeOnDrop for SecretKey {}

impl SecretKey {
    /// Generate a random [`SecretKey`] and write it to a file.
    pub fn random_and_write_to_file<R: RngCore + CryptoRng>(
        par: &Arc<BfvParameters>,
        rng: &mut R,
        file_name: &mut String,
    ) -> Self {
        let secret_key = Self::random(par, rng);

        let key_file = file_name.clone() + ".key";

        // Serialize the secret key into a JSON string
        let json = Self::to_bytes(&secret_key);

        // Write the JSON to a file
        let mut file = File::create(key_file).expect("Failed to create file");
        file.write_all(&json).expect("Failed to write to file");

        secret_key
    }

    /// Derive a random [`SecretKey`] and write it to a file.
    pub fn derive_and_write_to_file<R: RngCore + CryptoRng>(
        par: &Arc<BfvParameters>,
        rng: &mut R,
        der_key: &String,
        file_name: &mut String,
    ) -> Self {
        let secret_key = Self::derive(par, rng, der_key);

        let key_file = file_name.clone() + ".key";

        // Serialize the secret key into a JSON string
        let json = Self::to_bytes(&secret_key);

        // Write the JSON to a file
        let mut file = File::create(key_file).expect("Failed to create file");
        file.write_all(&json).expect("Failed to write to file");

        secret_key
    }

    /// Read a [`SecretKey`] from a file.
    pub fn read_from_file(file_name: String) -> Self {
        // read the secret key par from a file
        let key_file = file_name.clone() + ".key";

        // Read the JSON from the file
        let mut file = File::open(key_file).expect("Failed to open file");
        let mut json = String::new();
        file.read_to_string(&mut json)
            .expect("Failed to read from file");

        let secret_key = Self::try_from_bytes(json.as_bytes()).unwrap();

        secret_key
    }

    /// Generate a random [`SecretKey`].
    pub fn random<R: RngCore + CryptoRng>(par: &Arc<BfvParameters>, rng: &mut R) -> Self {
        let s_coefficients = sample_vec_cbd(par.degree(), par.variance, Some(rng), None).unwrap();

        let secret_key = Self::new(s_coefficients, par);

        secret_key
    }

    /// Derives a [`SecretKey`] using an input.
    pub fn derive<R: RngCore + CryptoRng>(
        par: &Arc<BfvParameters>,
        rng: &mut R,
        der_key: &String,
    ) -> Self {
        let s_coefficients =
            sample_vec_cbd(par.degree(), par.variance, Some(rng), Some(der_key)).unwrap();

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

    /// Serializes the [`SecretKey`] into a byte vector.
    pub fn to_bytes(&self) -> Vec<u8> {
        let binding = self.coeffs.to_vec();
        let coeffs_bytes = binding.as_slice();

        let coeffs_bytes_as_u8: &[u8] = unsafe {
            std::slice::from_raw_parts(
                coeffs_bytes.as_ptr() as *const u8,
                coeffs_bytes.len() * std::mem::size_of::<i64>(),
            )
        };

        let serial_sk = Data {
            par: self.par.to_bytes(),
            coeffs: coeffs_bytes_as_u8.to_vec(),
        };

        let json = serde_json::to_string(&serial_sk).expect("Failed to serialize to JSON");

        json.as_bytes().to_vec()
    }

    /// Deserializes the [`SecretKey`] from a byte vector.
    pub fn try_from_bytes(bytes: &[u8]) -> Result<Self> {
        // Parse the JSON back into a Data object
        let parsed_data: Data =
            serde_json::from_slice(&bytes).expect("Failed to parse JSON into Data object");

        // Access the individual data sets
        let par_bytes = parsed_data.par;
        let coeffs_bytes = parsed_data.coeffs;

        let par = Arc::new(BfvParameters::try_deserialize(&par_bytes)?);
        let coeffs_bytes_as_i64: &[i64] = unsafe {
            // Convert the &[u8] to a byte slice by transmuting the reference type
            std::slice::from_raw_parts(
                coeffs_bytes.as_ptr() as *const i64,
                coeffs_bytes.len() / std::mem::size_of::<i64>(),
            )
        };

        let coeffs = coeffs_bytes_as_i64.to_vec();

        let secret_key = Self::new(coeffs, &par);

        Ok(secret_key)
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
        let sk: SecretKey =
            SecretKey::random_and_write_to_file(&params, &mut rng, &mut "alice".to_string());
        assert_eq!(sk.par, params);
    }

    #[test]
    fn keygen_random_write_and_read() {
        let mut rng = thread_rng();
        let params = Arc::new(BfvParameters::default(1, 8));
        let sk: SecretKey =
            SecretKey::random_and_write_to_file(&params, &mut rng, &mut "alice".to_string());

        let sk_read = SecretKey::read_from_file("alice".to_string());

        assert_eq!(sk, sk_read);
    }

    #[test]
    fn keygen_derive_write_and_read() {
        let mut rng = thread_rng();
        let der_key: String =
            "8da4ef21b864d2cc526dbdb2a120bd2874c36c9d0a1fb7f8c63d7f7a8b41de8f".to_string();

        let params = Arc::new(BfvParameters::default(1, 8));

        let sk = SecretKey::derive_and_write_to_file(
            &params,
            &mut rng,
            &der_key,
            &mut "alice".to_string(),
        );

        let sk_read = SecretKey::read_from_file("alice".to_string());

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
