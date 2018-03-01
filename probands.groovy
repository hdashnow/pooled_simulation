load 'pipeline_config.groovy'
load 'analysis_stages.groovy'

// Settings
NSAMPLES=1 // Number of samples per capture (e.g. number in pools, or 1 if not pooled)
EXOME_TARGET="/group/bioi1/shared/genomes/hg19/exome_targets/nextera_rapid_capture_exome/target_region.bed"

run {
//    "%.fastq.gz" * [ fastqc ] +
    ~"(.*)_R[0-9][_.].*fastq.gz" * [ 
        set_sample_info + align + index_bam 
    ] +
    ~"(.*)_AGRF_.*bam" * [ 
        set_sample_info + merge_bams +
         //+ index_bam +
        realignIntervals + realign + index_bam +
        coverage +
        call_variants +
        annotate_vcf + 
        intersect_vcf
//        compress_vcf + index_vcf
    ]
}
