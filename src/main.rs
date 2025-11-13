// src/main.rs
use clap::Parser;
use vast::run;
use vast::Cli;

fn main() {
    // 1. Parse command-line arguments
    let cli = Cli::parse();

    // 2. Execute the core logic defined in src/lib.rs
    // Handle the VastResult<()> returned by the run function.
    if let Err(e) = run(cli) {
        // Print the error message to standard error (stderr)
        eprintln!("Error: {}", e);
        // Exit with a non-zero status code to signal failure
        std::process::exit(1);
    }
}