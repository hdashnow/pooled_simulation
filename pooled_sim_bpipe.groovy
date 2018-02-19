// Bpipe pipeline
// Harriet Dashnow 22 March 2015
// Takes mapped exome bams. Simulates pooled exomes.

load 'pipeline_config.groovy'

combined_bed="combined_target.bed"

// Settings
CHR='chr22'
sample_pools = [2, 4, 6, 8, 10]

run {
    intersect_targets +
    sample_pools * [
        set_pool + "%.bam" * [
            downsample_region
        ] + merge_bams + fix_header +
        realignIntervals + realign + index_bam + 
        call_variants +
        compress_vcf + index_vcf //+ 
//        cleanup
    ]
}
