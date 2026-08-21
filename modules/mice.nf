process mice {
    container "ghcr.io/gi-bielefeld/mice:f41a740_2026-06-10"

    input:
    path input_graph_gff
    val params

    output:
    path "mice_output/output.gff", emit: out_graph
    path "mice_output/partitions.txt", emit: out_parts
    path "mice_output/paths.txt", emit: out_paths

    script:
    def remove_dup = "-r ${params.remove_dup}"
    def quorum = "-q ${params.quorum}"
    def min_size = "-m ${params.min_size}"

    def no_group_by = params.no_group_by ? "-s" : ""
    def merge_with_dup = params.no_group_by ? "--merge-with-dup" : ""

    """
    mice ${remove_dup} ${quorum} ${min_size} ${no_group_by} ${merge_with_dup} ${input_graph_gff}
    """
}