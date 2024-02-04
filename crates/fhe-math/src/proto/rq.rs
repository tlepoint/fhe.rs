#![allow(missing_docs)]

#[allow(clippy::derive_partial_eq_without_eq)]
#[derive(Clone, PartialEq, ::prost::Message)]
pub struct Rq {
    #[prost(enumeration = "Representation", tag = "1")]
    pub representation: i32,
    #[prost(uint32, tag = "2")]
    pub degree: u32,
    #[prost(bytes = "vec", tag = "3")]
    pub coefficients: ::prost::alloc::vec::Vec<u8>,
}
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, PartialOrd, Ord, ::prost::Enumeration)]
#[repr(i32)]
pub enum Representation {
    Unknown = 0,
    Powerbasis = 1,
    Ntt = 2,
    Nttshoup = 3,
}
impl Representation {
    /// String value of the enum field names used in the ProtoBuf definition.
    ///
    /// The values are not transformed in any way and thus are considered stable
    /// (if the ProtoBuf definition does not change) and safe for programmatic
    /// use.
    pub fn as_str_name(&self) -> &'static str {
        match self {
            Representation::Unknown => "UNKNOWN",
            Representation::Powerbasis => "POWERBASIS",
            Representation::Ntt => "NTT",
            Representation::Nttshoup => "NTTSHOUP",
        }
    }
    /// Creates an enum from field names used in the ProtoBuf definition.
    pub fn from_str_name(value: &str) -> ::core::option::Option<Self> {
        match value {
            "UNKNOWN" => Some(Self::Unknown),
            "POWERBASIS" => Some(Self::Powerbasis),
            "NTT" => Some(Self::Ntt),
            "NTTSHOUP" => Some(Self::Nttshoup),
            _ => None,
        }
    }
}
