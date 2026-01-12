//! CLI configuration for PIR example runs.
// Expect indexing in examples for simplicity

use clap::Parser;

/// Command-line arguments for configuring a PIR example run.
#[derive(Parser)]
pub struct Cli {
    #[arg(
        long,
        help = "The number of elements in the database",
        default_value = "65536"
    )]
    /// Number of elements in the database.
    pub database_size: usize,

    #[arg(
        long,
        help = "The size of each database element",
        default_value = "1024"
    )]
    /// Size in bytes of each database element.
    pub element_size: usize,
}

#[allow(dead_code)]
fn main() {}
