load 'pipeline_config.groovy'
load 'analysis_stages.groovy'

// Settings
NSAMPLES=8 // Number of samples per capture (e.g. number in pools, or 1 if not pooled)
EXOME_TARGET="/group/bioi1/harrietd/pooled-parent/data/parentpool/SureSelectClinicalResearchExome_S06588914_Regions.bed"

run {
    //index_ref +
//    "%.fastq.gz" * [ fastqc ] +
    '%_R*.fastq.gz' * [
        set_sample_info +
        align + index_bam +
        realignIntervals + realign + index_bam +
        coverage +
        call_variants + 
        annotate_vcf +
        intersect_vcf
//        compress_vcf + index_vcf
    ]
}
