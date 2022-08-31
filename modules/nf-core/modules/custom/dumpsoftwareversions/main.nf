process CUSTOM_DUMPSOFTWAREVERSIONS {
    label 'process_low'

    input:
    path versions

    output:
    path "software_versions.yml"    , emit: yml

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    cat $versions > software_versions.yml
    """
}
