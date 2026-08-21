process ggcat {
    container "ghcr.io/gi-bielefeld/mice:f41a740_2026-06-10"

    input:
    path input_fasta_dir
    val params

    output:
    path "ggcat_output.gfa", emit: out_gfa_graph

    // TODO threads
    script:
    """
    ggcat build --gfa-v1 -p -k ${params.k} -s 1 \$(find \$(realpath ${input_fasta_dir}) -type f -name '*.fa' -o -name '*.fa.gz' -o -name '*.fna' -o -name '*.fna.gz' -o -name '*.fasta' -o -name '*.fasta.gz') -o ggcat_output.gfa
    """
}