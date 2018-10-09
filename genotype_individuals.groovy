// Bpipe pipeline
// Harriet Dashnow 19 Feb 2018
// Takes exome fastqs, aligns them, then does gvcf calling (to be later used for joint calling)

load 'pipeline_config.groovy'
load 'pipeline_stages.groovy'

// Settings
ploidy = 2
EXOME_TARGET="/group/bioi1/shared/genomes/hg19/exome_targets/NIMBLEGENV2.bed"

@preserve("*.bam")
align_bwa = {
    doc "Align reads with bwa mem algorithm."

    output.dir="align"

    def fastaname = get_fname(REF)
            //set -o pipefail
    from('fastq.gz', 'fastq.gz') produce(branch.sample + '.bam') {
        exec """
            bwa mem -M -t $threads
            -R "@RG\\tID:${flowcell}.${lane}\\tPL:$PLATFORM\\tPU:${flowcell}.${lane}.${library}\\tLB:${library}\\tSM:${sample}"
            $REF
            $inputs |
            samtools view -bSuh - | samtools sort -o $output.bam -T $output.bam.prefix
        """, "bwa"
    }
}

run {
    '%_R*.fastq.gz' * [
        set_fastq_info + set_sample_info +
        align_bwa + index_bam +
        dedup +
        coverage +
        call_variants_gvcf ]
}
