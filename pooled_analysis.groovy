load 'pipeline_config.groovy'
load 'pipeline_stages.groovy'

// Settings
EXOME_TARGET="/group/bioi1/harrietd/pooled-parent/data/parentpool/SureSelectClinicalResearchExome_S06588914_Regions.bed"
NSAMPLES=8 // Number of samples per capture (e.g. number in pools, or 1 if not pooled)
ploidy=NSAMPLES*2

run {
    //index_ref +
//    "%.fastq.gz" * [ fastqc ] +
    '%_R*.fastq.gz' * [
        set_fastq_info +
        align_bwa + index_bam +
        dedup +
        coverage +
        call_variants + 
        annotate_vcf +
        intersect_vcf
    ]
}
