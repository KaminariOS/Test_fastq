pub use clap::Parser;

/// Simple program to greet a person
#[derive(Parser, Debug)]
#[command(author, version, about, long_about = None)]
pub struct Args {
    /// Name of the person to greet
    #[arg(short, long)]
    pub filename: String,

    /// Number of times to greet
    #[arg(short, long, default_value_t = 32)]
    pub k: u8,
}

