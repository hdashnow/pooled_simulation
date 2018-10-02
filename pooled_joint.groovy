// perform GATK4 joint calling on probands and parent pool

load 'pipeline_config.groovy'
load 'pipeline_stages.groovy'

// Settings
// Using the intersection of the exome target regions for the pool and probands if different captures
INT_BED="/group/bioi1/harrietd/pooled-parent/filter/intersected/intersection_Nextera_SureSelectClinical.bed"

def get_fastq(directory) {

    def fastqs = []
    new File(directory).eachFileMatch(~/.*.fastq.gz/) { file ->  
    fastqs << file  
    }  
    //def fastqs = new groovy.util.FileNameFinder().getFileNames(directory, '*.fastq.gz')
    //println(fastqs)
    return fastqs
}

set_pool_info = {

    NSAMPLES=8 // Number of samples per capture in the pool (i.e. 1 if not pooled)
    branch.ploidy=NSAMPLES*2
    branch.EXOME_TARGET="/group/bioi1/harrietd/pooled-parent/data/parentpool/SureSelectClinicalResearchExome_S06588914_Regions.bed"
    
    //forward inputs
}

set_proband_info = {

    NSAMPLES=1 // Number of samples per capture in the pool (i.e. 1 if not pooled)
    branch.ploidy=NSAMPLES*2
    branch.EXOME_TARGET="/group/bioi1/shared/genomes/hg19/exome_targets/nextera_rapid_capture_exome/target_region.bed"

    //forward inputs
}

//variant_discovery = {
//        println inputs
//       ~'(.*?)_.*_R[12].fastq.gz' * [
//           set_sample_info +
//           ~'(.*)_R[12].fastq.gz' * [
//               set_fastq_info +
//               align_bwa + index_bam
//           ] + merge_lanes +
//        dedup +
//        coverage +
//        call_variants_gvcf
//        ]
//}

run {
//    "%.fastq.gz" * [ fastqc ] +

    [
        [ set_pool_info + 
            ~'(Parent.*)_.*_R[12].fastq.gz' * [
                set_sample_info + set_fastq_info +
                align_bwa + index_bam +
                dedup +
                coverage +
                call_variants_gvcf ]
        ],

        [ set_proband_info +
            ~'(0.*?)_AGRF_.*_R[12].fastq.gz' * [
                set_sample_info +
                ~'(0.*_AGRF_.*)_R[12].fastq.gz' * [
                    set_fastq_info +
                    align_bwa + index_bam
                ] + 
                merge_lanes +
                dedup +
                coverage +
                call_variants_gvcf 
            ]
        ]
    ] + 
    combine_gvcfs + 
    joint_calling +
    annotate_vcf +
    intersect_vcf
}
