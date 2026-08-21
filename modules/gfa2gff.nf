process gfa2gff {
    container "ghcr.io/gi-bielefeld/mice:f41a740_2026-06-10"

    input:
    path input_graph_gfa
    path input_fasta_dir
    val params

    output:
    path "mice_input_graph.gff", emit: out_gff_graph

    script:
    """
    gfa2gff ${params.k} ${input_graph_gfa} \$(find \$(realpath ${input_fasta_dir}) -type f -name '*.fa' -o -name '*.fa.gz' -o -name '*.fna' -o -name '*.fna.gz' -o -name '*.fasta' -o -name '*.fasta.gz') > mice_input_graph.gff
    """
}