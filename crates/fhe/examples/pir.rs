// Allow indexing in examples for simplicity
#![allow(clippy::indexing_slicing)]
#![allow(missing_docs)]

use clap::Parser;

#[derive(Parser)]
pub struct Cli {
    #[arg(
        long,
        help = "The number of elements in the database",
        default_value = "65536"
    )]
    pub database_size: usize,

    #[arg(
        long,
        help = "The size of each database element",
        default_value = "1024"
    )]
    pub element_size: usize,
}

#[allow(dead_code)]
fn main() {}
