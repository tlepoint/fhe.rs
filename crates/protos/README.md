### Setup

Install the `protobuf-codegen` crate
```bash
cargo install protobuf-codegen
```

### Generate the Rust files from the proto

Run
```bash
protoc --rust_out src/protos src/protos/*.proto
```
