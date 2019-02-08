// Bpipe pipeline
// Harriet Dashnow 22 March 2015
// Takes mapped exome bams. Simulates pooled exomes.

load 'pipeline_config.groovy'
load 'pipeline_stages.groovy'


// Settings
EXOME_TARGET="/group/bioi1/shared/genomes/hg19/exome_targets/NIMBLEGENV2.bed"

sample_pools = [2, 4, 6, 8, 10, 20, 40]

// Required for set_joint stage
proband_gvcfs = args.grep(~/.*gvcf$/)


run {

    // Sample from individuals to create simulated pools
    sample_pools * [
        set_pool + "%.bam" * [
            downsample_region
        ] + merge_bams + fix_header +
        dedup +
        call_variants_gvcf +
        // Do joint calling on simulated pool with its constituant samples
        set_joint +
        combine_gvcfs +
        joint_calling +
        [compare_joint, compare_joint_readfilter]
    ]
}
