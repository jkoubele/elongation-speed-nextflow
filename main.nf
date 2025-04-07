#!/usr/bin/env nextflow

params.bam_dir = "/cellfile/datapublic/jkoubele/celegans_mutants/BAM_example"
params.out_dir = "/cellfile/datapublic/jkoubele/celegans_mutants/example_output"
params.genome_folder = "/cellfile/datapublic/jkoubele/reference_genomes/WBcel235/"
params.script_folder = "/cellfile/datapublic/jkoubele/elongation-speed-nextflow/scripts"
params.docker_image_path = "/cellfile/datapublic/jkoubele/pol-II-analysis/docker_images/bioinfo_tools.tar"

process runDockerScript {
    tag "$sample_name"

    input:
    tuple val(sample_name), path(sample_folder)

    output:
    path "${sample_name}"  // This should be the output folder created inside the script

    publishDir "${params.out_dir}/${sample_name}", mode: 'copy'  // Create sample-specific subfolders

    script:
    """
    bash ${params.script_folder}/extract_intronic_reads.sh \\
        -i ${sample_folder} \\
        -o ${params.out_dir}/${sample_name} \\
        -d ${params.docker_image_path} \\
        -g ${params.genome_folder} \\
        -s ${params.script_folder}
    """
}


workflow {
    Channel
        .fromPath("${params.bam_dir}/*/downsampled.bam")
        .map { bam_file ->
            def folder = bam_file.getParent()
            def name = folder.getName()
            println "DEBUG: sample_name = ${name}, sample_folder = ${folder}"
            tuple(name, folder)
        }
        | runDockerScript
}
