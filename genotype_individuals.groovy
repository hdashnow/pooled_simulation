// Bpipe pipeline
// Harriet Dashnow 19 Feb 2018
// Takes mapped exome bams and calls variants

load 'pipeline_config.groovy'
load 'simulation_stages.groovy'

combined_bed="combined_target.bed"

set_ploidy = {
    branch.ploidy = 2
}

run {
    '%.bam' * [
        set_ploidy + 
        realignIntervals + realign + index_bam + 
        call_variants 
    ]
}
