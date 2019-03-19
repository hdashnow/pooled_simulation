load 'pipeline_config.groovy'
load 'pipeline_stages.groovy'

// Settings
EXOME_TARGET="/group/bioi1/shared/genomes/hg19/exome_targets/nextera_rapid_capture_exome/target_region.bed"
NSAMPLES=1 // Number of samples per capture (e.g. number in pools, or 1 if not pooled)
ploidy=NSAMPLES*2

run {
//    "%.fastq.gz" * [ fastqc ] +
    ~'(.*?)_.*_R[12].fastq.gz' * [
        set_sample_info +
        ~'(.*)_R[12].fastq.gz' * [
            set_fastq_info +
            align_bwa + index_bam
        ] + merge_lanes +
        dedup +
        //coverage +
        call_variants +
        annotate_vcf + 
        annotate_vep_table + annotate_vep_vcf +
        intersect_vcf +
        annotate_vep_table
    ]
}
