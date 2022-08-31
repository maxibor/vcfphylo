process GENOME2READS {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "bioconda::pysam=0.19.1" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pysam:0.19.1--py310hd89ff4b_0' :
        'quay.io/biocontainers/pysam:0.19.1--py310hd89ff4b_0' }"

    input:
    tuple val(meta), path(input)

    output:
    tuple val(meta), path("*.fastq") , emit: fastq

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    genome2reads.py \\
        $input
    """
}
