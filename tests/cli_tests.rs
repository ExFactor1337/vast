use clap::Parser;
use vast::{Cli, run}; // Import from the 'vast' crate (which is now src/lib.rs)

fn main() {
    let cli = Cli::parse();

    // Call run() and handle the final exit/error message.
    if let Err(e) = run(cli) {
        // Print the error message to standard error (stderr)
        eprintln!("Error: {}", e);
        // Exit with a non-zero status code to signal failure
        std::process::exit(1);
    }
}