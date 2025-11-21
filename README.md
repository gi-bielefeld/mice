# mice 

`mice` computes **synteny blocks** from genomes expressed as sequences of genomic elements.
These elements can come from a genome graph (e.g., unitigs of a compacted de Bruijn graph), or from any other segmentation such as k-mers, genes, or MUMs/MEMs.

The input of `mice` is a GFF file in which each feature has an `ID` attribute (1-based index) specifying the element used in the path spelling the genome or chromosome.

## Installation

`mice` is written in rust, therefore you only need cargo to install it:

```bash
cargo install --path .
```

## Quick start

We provide five *E. coli* genomes as an example dataset.

1. Use the provided graph  
   A precomputed `example/graph.gff.gz` is included.  
   Uncompress it and go directly to running `mice`.  

2. (*Optional*) Build the pangenome graph yourself

   Install `ggcat`:

   ```bash
   conda install -c conda-forge -c bioconda ggcat
   ```

   Build a compacted de Bruijn graph:

   ```bash
   ggcat build -k 31 -s 1 -l example/list.txt -o graph.gfa --gfa-v1
   ```

   Convert the graph to GFF:

   ```bash
   git clone https://github.com/lucaparmigiani/gfa2gff.git
   cd gfa2gff
   make
   gfa2gff 31 graph.gfa $(ls -1 example/*.fna.gz) > graph.gff
   ```

3. Run mice

   ```bash
   mice graph.gff
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
