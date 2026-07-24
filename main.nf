#!/usr/bin/env nextflow
process mice {
    container "ghcr.io/gi-bielefeld/mice:f41a740_2026-06-10"

    input:
    path input_graph_gff
    val params

    output:
    path "./output.gff", emit: out_graph
    path "./partitions.txt", emit: out_parts
    path "./paths.txt", emit: out_paths

    script:
    def remove_dup = "-r ${params.remove_dup}"
    def quorum = "-q ${params.quorum}"
    def min_size = "-m ${params.min_size}"

    def no_group_by = params.no_group_by ? "-s" : ""
    def merge_with_dup = params.no_group_by ? "--merge-with-dup" : ""
    // def dirty = params.no_group_by ? "" : ""

    """
    mice ${remove_dup} ${quorum} ${min_size} ${no_group_by} ${merge_with_dup} -o ./ ${input_graph_gff}
    """
}

workflow {
    main:
    def input_channel = channel.empty()

    def input_filepath = file(params.input_graph_gff)

    if (input_filepath.exists() && !input_filepath.isDirectory()) {
        input_channel = channel.fromPath(params.input_graph_gff, checkIfExists: true)
    }
    else {
        error("[ERROR] Invalid input: '${params.input_graph_gff}' is not a valid file path for --input_graph_gff.")
    }

    mice(input_channel,params)

    publish:
    mice_out_graph = mice.out.out_graph
    mice_out_parts = mice.out.out_parts
    mice_out_paths = mice.out.out_paths
}

output {
    mice_out_graph {
        mode "copy"
    }
    mice_out_parts {
        mode "copy"
    }
    mice_out_paths {
        mode "copy"
    }
}