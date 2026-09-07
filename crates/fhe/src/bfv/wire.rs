//! Traits used for the BFV homomorphic encryption scheme.

use crate::Result;
use crate::bfv::Parameters;

/// Internal conversion of validated protobuf DTOs with explicit parameter
/// binding.
pub(crate) trait FromProto<T>
where
    Self: Sized,
{
    /// Attempt to convert the `value` with a specific parameter.
    fn from_proto(value: T, par: &Parameters, limits: &crate::DecodeLimits) -> Result<Self>;
}

use crate::{
    DecodeLimits, Error,
    error::{SerializationError, SerializedObject},
};
use prost::encoding::{DecodeContext, WireType, decode_key, decode_varint, skip_field};

#[derive(Clone, Copy)]
enum Shape {
    Parameters,
    Secret,
    Ciphertext,
    PublicKey,
    KeySwitch,
    Relin,
    Galois,
    Evaluation,
    Rgsw,
}

/// Scan borrowed fields before prost allocates DTO vectors. Count all field
/// occurrences, including duplicates that protobuf would merge or overwrite.
pub(crate) fn preflight(
    bytes: &[u8],
    object: SerializedObject,
    par: Option<&Parameters>,
    limits: &DecodeLimits,
) -> crate::Result<()> {
    let (degree, moduli) = par
        .map(|p| (p.degree(), p.moduli().len()))
        .unwrap_or((0, 0));
    limits.check_context(bytes.len(), degree, moduli)?;
    let shape = match object {
        SerializedObject::Parameters => Shape::Parameters,
        SerializedObject::SecretKey => Shape::Secret,
        SerializedObject::Ciphertext => Shape::Ciphertext,
        SerializedObject::PublicKey => Shape::PublicKey,
        SerializedObject::RelinearizationKey => Shape::Relin,
        SerializedObject::EvaluationKey => Shape::Evaluation,
        SerializedObject::RgswCiphertext => Shape::Rgsw,
    };
    let mut budget = Budget {
        limits,
        degree,
        moduli,
        polynomials: 0,
        object,
    };
    budget.scan(bytes, shape)
}

struct Budget<'a> {
    limits: &'a DecodeLimits,
    degree: usize,
    moduli: usize,
    polynomials: usize,
    object: SerializedObject,
}
impl Budget<'_> {
    fn decode_error(&self, source: prost::DecodeError) -> Error {
        SerializationError::Decode {
            object: self.object,
            source,
        }
        .into()
    }

    fn charge(&mut self, count: usize) -> crate::Result<()> {
        self.polynomials = self.polynomials.saturating_add(count);
        self.limits
            .check_polynomials(self.polynomials, self.degree, self.moduli)?;
        Ok(())
    }

    fn scan(&mut self, mut input: &[u8], shape: Shape) -> crate::Result<()> {
        let mut scalars = 0usize;
        let mut first_components = 0usize;
        let mut seeded = false;
        if matches!(shape, Shape::Evaluation) && self.degree != 0 {
            self.charge(self.degree.ilog2() as usize)?;
        }
        while !input.is_empty() {
            let (tag, wire) = decode_key(&mut input).map_err(|e| self.decode_error(e))?;
            if wire == WireType::Varint {
                let value = decode_varint(&mut input).map_err(|e| self.decode_error(e))?;
                match (shape, tag) {
                    (Shape::Parameters, 1) => {
                        // Prost uint32 decoding casts the varint to u32.
                        self.limits.check(
                            "degree",
                            value as u32 as usize,
                            self.limits.max_degree,
                        )?;
                    }
                    (Shape::Parameters, 2) | (Shape::Secret, 1) => {
                        scalars = scalars.saturating_add(1);
                    }
                    _ => {}
                }
            } else if wire == WireType::LengthDelimited {
                // Let prost validate the complete field without allocating.
                let mut rest = input;
                skip_field(wire, tag, &mut rest, DecodeContext::default())
                    .map_err(|e| self.decode_error(e))?;
                let len = decode_varint(&mut input).map_err(|e| self.decode_error(e))? as usize;
                let (data, _) = input.split_at(len);
                input = rest;
                match (shape, tag) {
                    (Shape::Parameters, 2) | (Shape::Secret, 1) => {
                        let mut packed = data;
                        while !packed.is_empty() {
                            decode_varint(&mut packed).map_err(|e| self.decode_error(e))?;
                            scalars = scalars.saturating_add(1);
                        }
                    }
                    (Shape::Parameters, 5) => self.limits.check(
                        "plaintext modulus bytes",
                        data.len(),
                        self.limits.max_plaintext_bytes,
                    )?,
                    (Shape::Ciphertext, 1) | (Shape::KeySwitch, 1 | 2) => {
                        self.charge(1)?;
                        if tag == 1 {
                            first_components = first_components.saturating_add(1);
                        }
                    }
                    (Shape::Ciphertext, 2) | (Shape::KeySwitch, 3) => {
                        seeded |= !data.is_empty();
                    }
                    (Shape::PublicKey, 1) => self.scan(data, Shape::Ciphertext)?,
                    (Shape::Relin | Shape::Galois, 1) | (Shape::Rgsw, 1 | 2) => {
                        self.scan(data, Shape::KeySwitch)?
                    }
                    (Shape::Evaluation, 2) => {
                        // Also bounds the DTO vector when entries are empty.
                        self.charge(1)?;
                        self.scan(data, Shape::Galois)?;
                    }
                    _ => {}
                }
            } else {
                skip_field(wire, tag, &mut input, DecodeContext::default())
                    .map_err(|e| self.decode_error(e))?;
            }
            match shape {
                Shape::Parameters => {
                    self.limits
                        .check("moduli", scalars, self.limits.max_moduli)?
                }
                Shape::Secret => {
                    self.limits
                        .check("secret coefficients", scalars, self.limits.max_degree)?
                }
                Shape::Ciphertext
                | Shape::PublicKey
                | Shape::KeySwitch
                | Shape::Relin
                | Shape::Galois
                | Shape::Evaluation
                | Shape::Rgsw => {}
            }
        }
        if seeded {
            self.charge(if matches!(shape, Shape::KeySwitch) {
                first_components
            } else {
                1
            })?;
        }
        Ok(())
    }
}
