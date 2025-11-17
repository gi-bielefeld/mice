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

    /// Remove an element if it occurs more than x times in any genome. Use 0 to disable removal.
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

    /// Minimum element length (in bp) to keep elements that were not merged after the first compression.
    ///
    /// Compression is first performed. Elements that remain unmerged and are shorter than this
    /// length are then removed, and compression is performed again.
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

    /// If set every path is treated as its own genome
    #[arg(short = 's', long = "no-group-by", default_value_t = true, action = ArgAction::SetFalse)]
    pub group_by: bool,

    /// Avoid marking duplicated elements in each genome and merge them like they were unique
    #[arg(long = "dirty", hide = true, action = ArgAction::SetTrue)]
    pub dirty: bool,
}
