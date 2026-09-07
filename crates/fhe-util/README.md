# fhe-util [![crate version](https://img.shields.io/crates/v/fhe-util.svg)](https://crates.io/crates/fhe-util) [![documentation](https://docs.rs/fhe-util/badge.svg)](https://docs.rs/fhe-util)

Utility functions for the [`fhe.rs`](https://github.com/tlepoint/fhe.rs) ecosystem.

The crate contains helper routines such as primality testing, centered binomial sampling, modular arithmetic helpers and other small utilities relied upon by the [`fhe`](https://crates.io/crates/fhe) and `fhe-math` crates.

## Installation

```toml
[dependencies]
fhe-util = "0.1.2"
```

## Checked bit transcoding

Bit-transcoding helpers return `Result<_, TranscodeError>` for invalid widths,
length overflow, or words that would be truncated. Packing uses least-significant
bits first and zero-pads the final byte. `transcode_to_bytes_into` appends only
after validating all input, leaving its output unchanged on an error.

`transcode_from_bytes_exact(bytes, width, count)` rejects missing/extra bytes and
nonzero padding. Use it for a packed message with a known word count.
`transcode_from_bytes(bytes, width)` interprets every byte bit as data and
zero-extends a partial final word. `transcode_bidirectional` also zero-extends its
final word and does not retain the original bit count. No helper silently masks
out-of-range input words.

`sample_vec_cbd` returns the typed `InvalidVariance` error for variances outside
1..=32. The sampling algorithm and randomness consumption are unchanged.

## Decode limits

`DecodeLimits` and `DecodeLimitError` provide resource bounds used by BFV and
polynomial imports. They check encoded size, context dimensions, and expanded
polynomial storage before allocation. Limits can be tightened or raised through
explicit fields; arithmetic overflow is always rejected. These bounds are not a
quota for allocator overhead or total process RSS.

## Testing

```bash
cargo test -p fhe-util
```

## License

This project is licensed under the [MIT license](https://opensource.org/licenses/MIT).

## Security / Stability

The code in this crate has not undergone an independent security audit.
Use at your own risk.
