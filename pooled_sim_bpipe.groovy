// Bpipe pipeline
// Harriet Dashnow 22 March 2015
// Takes mapped exome bams. Simulates pooled exomes.

load 'pipeline_config.groovy'
load 'pipeline_stages.groovy'
EXOME_TARGET="/group/bioi1/shared/genomes/hg19/exome_targets/NIMBLEGENV2.bed"

// Settings
//CHR='chr22'
sample_pools = [2, 4, 6, 8, 10]

run {
//    intersect_targets +
    sample_pools * [
        set_pool + "%.bam" * [
            downsample_region
        ] + merge_bams + fix_header +
        dedup +
        call_variants
//        filter_vcf_qual
//        compress_vcf + index_vcf + 
//        cleanup
    ] + compare_sim
}
