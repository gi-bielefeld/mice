use anyhow::Result;
use clap::Parser;

mod cli;
mod compression;
mod collections;
mod io;

fn main() -> Result<()> {
    let args = cli::Cli::parse();
    compression::run_mice(&args)?;

    Ok(())
}
