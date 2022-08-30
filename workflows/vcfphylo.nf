/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// Validate input parameters
WorkflowVcfphylo.initialise(params, log)

// TODO nf-core: Add all file path parameters for the pipeline to the list below
// Check input path parameters to see if they exist
def checkPathParamList = [ params.input, params.multiqc_config, params.fasta ]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Check mandatory parameters
if (params.input) { ch_input = file(params.input) } else { exit 1, 'Input samplesheet not specified!' }

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    INPUT CHANNELS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

Channel
    .fromPath(params.ref_genomes)
    .map {
        it -> [['id':it.baseName, 'single_end': true], it]
    }
    .set { ref_genomes }
representative_genome = Channel.fromPath(params.representative_genome)
samp_vcf = Channel.fromPath(params.samp_vcf)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//
include { CUSTOM_DUMPSOFTWAREVERSIONS } from '../modules/nf-core/modules/custom/dumpsoftwareversions/main'
include { GUNZIP as GUNZIP_REF ; GUNZIP as GUNZIP_REP } from '../modules/nf-core/modules/gunzip/main'
include { SAMTOOLS_FAIDX } from '../modules/nf-core/modules/samtools/faidx/main'
include { SAMTOOLS_MPILEUP } from '../modules/nf-core/modules/samtools/mpileup/main'
include { BOWTIE2_INDEX } from '../modules/nf-core/modules/bowtie2/index/main'
include { BOWTIE2_ALIGN } from '../modules/nf-core/modules/bowtie2/align/main'
include { GENOME2READS } from '../modules/local/genome2reads'
include { BCFTOOLS_CALL } from '../modules/local/bcftools/call'
include { MULTIVCFANALYZER } from '../modules/nf-core/modules/multivcfanalyzer/main'
include { IQTREE } from '../modules/nf-core/modules/iqtree/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Info required for completion email and summary

workflow VCFPHYLO {

    ch_versions = Channel.empty()

    ref_genomes
        .branch {
            decompressed: it[1].toString().tokenize(".")[-1] != 'gz'
            compressed: it[1].toString().tokenize(".")[-1] == 'gz'
        }
        .set { genomes_fork }

    GUNZIP_REF (
        genomes_fork.compressed
    )

    GUNZIP_REF.out.gunzip
        .mix( genomes_fork.decompressed )
        .set { ref_genomes_pre_processed }

    GENOME2READS (
        ref_genomes_pre_processed
    )

    GUNZIP_REP (
        representative_genome
    )

    SAMTOOLS_FAIDX (
        GUNZIP_REP.out.gunzip
    )

    BOWTIE2_INDEX (
        GUNZIP_REP.out.gunzip
    )

    BOWTIE2_ALIGN (
        ref_genomes_pre_processed.out.fastq,
        BOWTIE2_INDEX.out.index.first(),
        false,
        true
    )

    SAMTOOLS_MPILEUP(
        BOWTIE2_ALIGN.out.bam
            .map {
                meta, bam -> [meta, bam, []]
            } ,
        GUNZIP_REP.out.gunzip
    )

    BCFTOOLS_CALL (
        SAMTOOLS_MPILEUP.out.mpileup
    )

    MULTIVCFANALYZER (
        BCFTOOLS_CALL.out.vcf.collect()
        .mix (samp_vcf.collect())
        .collect(), // vcfs
        GUNZIP_REP.out.gunzip, // fasta
        [], // snpeff_results
        [], // gff
        [], // allele_freqs
        30, // genotype_quality
        5, // coverage
        0.8, // homozygous_freq
        0.2, // heterozygous_freq
        [] //gff_exclude
    )

    IQTREE (
        MULTIVCFANALYZER.out.snp_alignment,
        []
    )

    CUSTOM_DUMPSOFTWAREVERSIONS (
        ch_versions.unique().collectFile(name: 'collated_versions.yml')
    )

    ch_versions    = ch_versions.mix(MULTIQC.out.versions)
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COMPLETION EMAIL AND SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow.onComplete {
    if (params.email || params.email_on_fail) {
        NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, multiqc_report)
    }
    NfcoreTemplate.summary(workflow, params, log)
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
