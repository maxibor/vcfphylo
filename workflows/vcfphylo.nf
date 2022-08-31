/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)


// TODO nf-core: Add all file path parameters for the pipeline to the list below
// Check input path parameters to see if they exist
// Check mandatory parameters
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    INPUT CHANNELS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

Channel
    .fromPath(params.ref_genomes)
    .map {
        it -> [['id':it.simpleName, 'single_end': true], it]
    }
    .set { ref_genomes }

Channel
    .fromPath(params.representative_genome)
    .map {
        it -> [['id':it.simpleName, 'single_end': true], it]
    }
    .set { representative_genome }

samp_vcf = params.samp_vcf ? Channel.fromPath(params.samp_vcf) : Channel.empty()

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
include { SAMTOOLS_INDEX } from '../modules/nf-core/modules/samtools/index/main'
include { BOWTIE2_BUILD } from '../modules/nf-core/modules/bowtie2/build/main'
include { BOWTIE2_ALIGN } from '../modules/nf-core/modules/bowtie2/align/main'
include { GENOME2READS } from '../modules/local/genome2reads'
include { PICARD_CREATESEQUENCEDICTIONARY } from '../modules/nf-core/modules/picard/createsequencedictionary/main'
include { GATK_UNIFIEDGENOTYPER } from '../modules/nf-core/modules/gatk/unifiedgenotyper/main'
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
        ref_genomes_pre_processed,
    )

    GUNZIP_REP (
        representative_genome
    )

    SAMTOOLS_FAIDX (
        GUNZIP_REP.out.gunzip
    )

    ch_versions = ch_versions.mix( SAMTOOLS_FAIDX.out.versions )

    BOWTIE2_BUILD (
        GUNZIP_REP.out.gunzip
    )

    BOWTIE2_ALIGN (
        GENOME2READS.out.fastq,
        BOWTIE2_BUILD.out.index.first(),
        false,
        true
    )

    ch_versions = ch_versions.mix( BOWTIE2_ALIGN.out.versions )

    SAMTOOLS_INDEX(
        BOWTIE2_ALIGN.out.bam
    )

    BOWTIE2_ALIGN.out.bam
            .join(SAMTOOLS_INDEX.out.bai)
            .view()

    PICARD_CREATESEQUENCEDICTIONARY (
        GUNZIP_REP.out.gunzip
    )

    ch_versions = ch_versions.mix( PICARD_CREATESEQUENCEDICTIONARY.out.versions )

    GATK_UNIFIEDGENOTYPER (
        BOWTIE2_ALIGN.out.bam
            .join(SAMTOOLS_INDEX.out.bai),
        GUNZIP_REP.out.gunzip
            .map {
            it -> it[1]
            }
            .first(),
        SAMTOOLS_FAIDX.out.fai
            .map {
            it -> it[1]
            }
            .first(),
        PICARD_CREATESEQUENCEDICTIONARY.out.reference_dict
            .map {
            it -> it[1]
            }
            .first(),
        [],
        [],
        [],
        []
    )

    ch_versions = ch_versions.mix( GATK_UNIFIEDGENOTYPER.out.versions )

    MULTIVCFANALYZER (
        GATK_UNIFIEDGENOTYPER.out.vcf
        .map {
            it -> [it[1]]
        }
        .collect()
        .mix (samp_vcf.collect().ifEmpty([]))
        .collect(), // vcfs
        GUNZIP_REP.out.gunzip
            .map {
                it -> [it[1]]
            }.first(), // fasta
        [], // snpeff_results
        [], // gff
        [], // allele_freqs
        30, // genotype_quality
        5, // coverage
        0.8, // homozygous_freq
        0.2, // heterozygous_freq
        [] //gff_exclude
    )

    ch_versions = ch_versions.mix( MULTIVCFANALYZER.out.versions )

    IQTREE (
        MULTIVCFANALYZER.out.snp_genome_alignment,
        []
    )

    ch_versions = ch_versions.mix( IQTREE.out.versions )

    CUSTOM_DUMPSOFTWAREVERSIONS (
        ch_versions.unique().collectFile(name: 'collated_versions.yml')
    )

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
