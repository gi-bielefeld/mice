use anyhow::Result;
use clap::Parser;

mod cli;
mod compression;
mod collections;
mod io;

fn main() -> Result<()> {
    let args = cli::Cli::parse();
    compression::run_mice(&args)?;
    // compression::run_mice_hardcoded("/home/kml/phd/code/mice/example/graph.gff", "/home/kml/phd/code/mice/test_output",None);

    Ok(())
}
