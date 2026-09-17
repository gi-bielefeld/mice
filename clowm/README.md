# MICE on CloWM

Genomic rearrangements are major drivers of evolution and genetic disease, and a key focus of comparative genomics analyses.
Studying rearrangements requires segmenting genomes of interest into conserved regions, called **synteny blocks**, that highlight structural differences between genomes.
Synteny blocks are typically defined from annotated genes or derived as a by-product of whole-genome alignment.
However, as these procedures are heuristic and do not explicitly model rearrangements, they can obscure real variation, create false similarities, and affect phylogenetic inference.

MICE is an efficient, alignment-free program that instead computes synteny blocks according to an exact algorithm. This algorithm guarantees the resulting blocks have several advantageous properties, including no obscured intra-block rearrangements. It computes these blocks from genomes expressed as sequences of genomic elements, which can come from a genome graph (e.g., unitigs of a compacted de Bruijn graph), or from any other segmentation such as k-mers, genes, or MUMs/MEMs.

---

The input of MICE is a GFF file in which each feature has an `ID` attribute (1-based index) specifying the element used in the path spelling the genome or chromosome.

This workflow provides both the software implemention [`mice`](https://github.com/gi-bielefeld/mice) as well as two upstream tools [`ggcat`](https://github.com/algbio/ggcat) and [`gfa2gff`](https://github.com/lucaparmigiani/gfa2gff) which allow for starting the workflow with raw FASTA files, if you don't already have an appropriate GFF graph for direct input to `mice`.

For more info on input and output options, see the Usage tab.