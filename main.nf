#!/usr/bin/env nextflow
process mice {
    container "mice"

    input:
    path input_graph
    val params

    output:
    path params.out_dir, emit: out_dir

    script:
    def remove_dup = "-r ${params.remove_dup}"
    def quorum = "-q ${params.quorum}"
    def min_size = "-m ${params.min_size}"

    def no_group_by = params.no_group_by ? "-s" : ""
    def merge_with_dup = params.no_group_by ? "--merge-with-dup" : ""
    // def dirty = params.no_group_by ? "" : ""

    """
    mice ${remove_dup} ${quorum} ${min_size} ${no_group_by} ${merge_with_dup} -o ${params.out_dir} ${input_graph}
    """
}

workflow {
    main:
    def input_channel = channel.empty()

    def input_filepath = file(params.input_graph)

    if (input_filepath.exists() && !input_filepath.isDirectory()) {
        input_channel = channel.fromPath(params.input_graph, checkIfExists: true)
    }
    else {
        error("[ERROR] Invalid input: '${params.input_graph}' is not a valid file path for --input_graph.")
    }

    mice(input_channel,params)

    publish:
    mice_output = mice.out.out_dir
}

output {
    mice_output {
        mode "copy"
    }
}