// Bpipe pipeline
// Harriet Dashnow 19 Feb 2018
// Takes mapped exome bams and calls variants

load 'pipeline_config.groovy'
load 'pipeline_stages.groovy'

EXOME_TARGET="/group/bioi1/shared/genomes/hg19/exome_targets/NIMBLEGENV2.bed"

set_ploidy = {
    branch.ploidy = 2
}

run {
    '%_R*.fastq.gz' * [
        set_fastq_info +
        align_bwa + index_bam +
        dedup +
        set_ploidy + 
        call_variants + filter_vcf_qual
    ]
}
