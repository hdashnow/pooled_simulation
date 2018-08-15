// Bpipe pipeline
// Harriet Dashnow 22 March 2015
// Takes mapped exome bams. Simulates pooled exomes.
// Simulates pools by simply taking all reads from the individual samples and pooling them together
// with no downsampling. e.g. a pool of 10 samples, each with 100x coverage would end up with 1000x coverage

load 'pipeline_config.groovy'
load 'simulation_stages.groovy'
combined_bed="/group/bioi1/shared/genomes/hg19/exome_targets/NIMBLEGENV2.bed"

// Settings
sample_pools = [2, 4, 6, 8, 10]

run {
//    intersect_targets +
    sample_pools * [
        set_pool +
        merge_bams + fix_header +
        call_variants +
        filter_vcf_qual
//        compress_vcf + index_vcf + 
//        cleanup
    ]
}
