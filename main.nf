#!/usr/bin/env nextflow
include { mice } from './modules/mice.nf'

def parseS3Path(String path) {
    def s3_pattern = ~/^s3:\/\/([a-z0-9.-]{3,63})(?:\/(.*))?$/
    def matcher = path =~ s3_pattern

    if (matcher.matches()) {
        // matcher[0][1] is capture group 1 (bucket)
        // matcher[0][2] is capture group 2 (key)
        return [
            isValid: true, 
            bucket: matcher[0][1], 
            key: matcher[0][2] ?: "" // Handle empty key gracefully
        ]
    } else {
        return [isValid: false, bucket: null, key: null]
    }
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

    // output s3 dir - URI validation
    def s3_output_info = parseS3Path(params.out_dir)
    if (!s3_output_info.isValid) {
        error("Invalid output: '${params.out_dir}' is not a valid s3 URI for the output directory.")
    }

    mice(input_channel,params)

    publish:
    mice_out_graph = mice.out.out_graph
    mice_out_parts = mice.out.out_parts
    mice_out_paths = mice.out.out_paths

    onComplete:
    log.info("Pipeline finished. Cleaning up the artefact '.' folder in output dir ${params.output_dir}")
    try {
        to_delete = file("${params.output_dir}/.")
        to_delete.deleteDir()
        log.info("Successfully removed the '.' S3 artefact.")
    } catch (Exception e) {
        log.warn("Cleanup note: Could not remove the artefact '.' directory (${e.message}).")
    }
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