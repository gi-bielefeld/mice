# mice 

`mice` computes **synteny blocks** from genomes expressed as sequences of genomic elements.
These elements can come from a genome graph (e.g., unitigs of a compacted de Bruijn graph), or from any other segmentation such as k-mers, genes, or MUMs/MEMs.

The input of `mice` is a GFF file in which each feature has an `ID` attribute (1-based index) specifying the element used in the path spelling the genome or chromosome.

## Installation

`mice` is written in rust, therefore you only need cargo to install it:

```bash
cargo install --path .
```

## Usage

```bash
mice [OPTIONS] <GRAPH_INPUT>
```

* `<GRAPH_INPUT>` – input graph file (GFF or GFA with path representing genomes)

### Options

* `-o, --out-dir <DIR>`
  Output directory (default: `mice_output`)

* `-r, --remove-dup <X>`
  Remove an element if it occurs more than *X* times in any genome (`0` = disable, default: `0`)

* `-m, --min-size <bp>`
  After first compression, drop unmerged elements shorter than `<bp>` base pairs, then recompress (default: `0`)

* `-s, --no-group-by`
  Treat every path as its own genome

* `-h, --help`, `-V, --version`

## Output

In `<OUT_DIR>` MICE writes:

* `output.gff`: block annotations (GFF)
* `paths.txt`: genomes rewritten as synteny blocks
* `partitions.txt`: each synteny block which element it contains
