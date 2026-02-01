#![allow(missing_docs)]

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let proto_path = "src/proto/bfv.proto";
    let proto_dir = "src/proto";

    println!("cargo:rerun-if-changed={proto_path}");

    let mut config = prost_build::Config::new();
    config.compile_protos(&[proto_path], &[proto_dir])?;
    Ok(())
}
