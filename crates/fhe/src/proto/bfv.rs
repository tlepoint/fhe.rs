#![allow(missing_docs)]

#[allow(clippy::derive_partial_eq_without_eq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
#[derive(Clone, PartialEq, ::prost::Message)]
pub struct Ciphertext {
    #[prost(bytes = "vec", repeated, tag = "1")]
    pub c: ::prost::alloc::vec::Vec<::prost::alloc::vec::Vec<u8>>,
    #[prost(bytes = "vec", tag = "2")]
    pub seed: ::prost::alloc::vec::Vec<u8>,
    #[prost(uint32, tag = "3")]
    pub level: u32,
}
#[allow(clippy::derive_partial_eq_without_eq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
#[derive(Clone, PartialEq, ::prost::Message)]
pub struct RgswCiphertext {
    #[prost(message, optional, tag = "1")]
    pub ksk0: ::core::option::Option<KeySwitchingKey>,
    #[prost(message, optional, tag = "2")]
    pub ksk1: ::core::option::Option<KeySwitchingKey>,
}
#[allow(clippy::derive_partial_eq_without_eq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
#[derive(Clone, PartialEq, ::prost::Message)]
pub struct KeySwitchingKey {
    #[prost(bytes = "vec", repeated, tag = "1")]
    pub c0: ::prost::alloc::vec::Vec<::prost::alloc::vec::Vec<u8>>,
    #[prost(bytes = "vec", repeated, tag = "2")]
    pub c1: ::prost::alloc::vec::Vec<::prost::alloc::vec::Vec<u8>>,
    #[prost(bytes = "vec", tag = "3")]
    pub seed: ::prost::alloc::vec::Vec<u8>,
    #[prost(uint32, tag = "4")]
    pub ciphertext_level: u32,
    #[prost(uint32, tag = "5")]
    pub ksk_level: u32,
    #[prost(uint32, tag = "6")]
    pub log_base: u32,
}
#[allow(clippy::derive_partial_eq_without_eq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
#[derive(Clone, PartialEq, ::prost::Message)]
pub struct RelinearizationKey {
    #[prost(message, optional, tag = "1")]
    pub ksk: ::core::option::Option<KeySwitchingKey>,
}
#[allow(clippy::derive_partial_eq_without_eq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
#[derive(Clone, PartialEq, ::prost::Message)]
pub struct GaloisKey {
    #[prost(message, optional, tag = "1")]
    pub ksk: ::core::option::Option<KeySwitchingKey>,
    #[prost(uint32, tag = "2")]
    pub exponent: u32,
}
#[allow(clippy::derive_partial_eq_without_eq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
#[derive(Clone, PartialEq, ::prost::Message)]
pub struct EvaluationKey {
    #[prost(message, repeated, tag = "2")]
    pub gk: ::prost::alloc::vec::Vec<GaloisKey>,
    #[prost(uint32, tag = "3")]
    pub ciphertext_level: u32,
    #[prost(uint32, tag = "4")]
    pub evaluation_key_level: u32,
}
#[allow(clippy::derive_partial_eq_without_eq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
#[derive(Clone, PartialEq, ::prost::Message)]
pub struct Parameters {
    #[prost(uint32, tag = "1")]
    pub degree: u32,
    #[prost(uint64, repeated, tag = "2")]
    pub moduli: ::prost::alloc::vec::Vec<u64>,
    #[prost(uint64, tag = "3")]
    pub plaintext: u64,
    #[prost(uint32, tag = "4")]
    pub variance: u32,
}
#[allow(clippy::derive_partial_eq_without_eq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
#[derive(Clone, PartialEq, ::prost::Message)]
pub struct PublicKey {
    #[prost(message, optional, tag = "1")]
    pub c: ::core::option::Option<Ciphertext>,
}
