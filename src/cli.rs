use clap::{builder::ValueParser, ArgAction, Parser};

#[derive(Parser, Debug)]
#[command(
    name = "graph-input",
    version,
    about = "Parse paths from a GFF/GFA file"
)]
pub struct Cli {
    /// Input graph file
    pub graph_input: String,

    /// Output directory
    #[arg(short = 'o', long = "out-dir", default_value = "mice_output")]
    pub out_dir: String,

    /// Remove an element if it occurs at least `remove-dup` times in any genome. Use 0 to disable removal.
    /// A value of 1 will result in an empty input; if you want this option enabled, use 2 or above.
    /// See `--no-group-by` for alternative way of defining "genome".
    #[arg(
        short = 'r',
        long = "remove-dup",
        default_value_t = 0,
        value_parser = ValueParser::new(|s: &str| -> Result<usize, String> {
            let v: usize = s.parse().map_err(|_| "Expected a positive integer".to_string())?;
            if v == 1 {
                Err("Value 1 is not allowed".to_string())
            } else {
                Ok(v)
            }
        })
    )]
    pub remove_duplicates: usize,

    /// Minimum number of genomes an element must appear in to be kept. Use 0 to disable filtering.
    /// See `--no-group-by` for alternative way of defining "genome".
    #[arg(
        short = 'q',
        long = "quorum",
        default_value_t = 0,
        value_parser = ValueParser::new(|s: &str| -> Result<usize, String> {
            let v: usize = s.parse().map_err(|_| "Expected a positive integer".to_string())?;
            Ok(v)
        })
    )]
    pub quorum: usize,

    /// Compression is first performed. Elements that remain unmerged and are shorter than this
    /// length (in basepairs) are then removed, and compression is performed again.
    #[arg(
        short = 'm',
        long = "min-size",
        default_value_t = 0,
        value_parser = ValueParser::new(|s: &str| -> Result<usize, String> {
            let v: usize = s.parse().map_err(|_| "Expected a positive integer".to_string())?;
            Ok(v)
        })
    )]
    pub min_size: usize,

    /// If true, every sequence region (chromosome, contig, FASTA header etc.) is treated as its own genome.
    #[arg(short = 's', long = "no-group-by", default_value_t = true, action = ArgAction::SetFalse)]
    pub group_by: bool,

    /// Avoid marking duplicated elements in each genome and merge them like they were unique.
    /// Setting this option activates the MICE (dup) mode described in the paper. 
    /// Keeping it unset instead keeps MICE (Bp bij) mode active.
    #[arg(long = "merge-with-dup", action = ArgAction::SetTrue)]
    pub merge_with_dup: bool,

    /// Add each block occurrence's (phrase's) elements to output.gff as a path attribute, 
    /// oriented so the block's canonical anchor reads '+' (filtered elements are omitted).
    #[arg(short = 'p', long = "phrase-paths", action = ArgAction::SetTrue)]
    pub write_phrases: bool,
}