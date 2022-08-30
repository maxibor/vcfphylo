process GENOME2READS {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "conda-forge::openjdk=8.0.312" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/openjdk:8.0.121' :
        'quay.io/biocontainers/openjdk:8.0.121' }"

    input:
    tuple val(meta), path(input)

    output:
    tuple val(meta), path("*.fastq") , emit: fastq

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def read_length = 150
    def tiling_density = 2
    """
    java -jar Genome2Reads.jar \\
        $input \\
        $read_length \\
        $tiling_density \\
        > ${prefix}.fastq
    """
}
