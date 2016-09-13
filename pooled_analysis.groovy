BASE="/vlsci/VR0320/shared/hdashnow/MelbGenomicsPipelineRepo"

// Ref files
// Note, dbSNP is an old version, may want to update
REFBASE="/group/bioi1/shared/genomes/hg19/gatk"
REF="$REFBASE/gatk.ucsc.hg19.fasta"
DBSNP="$REFBASE/dbsnp_132.hg19.vcf "
GOLD_STANDARD_INDELS="$REFBASE/Mills_and_1000G_gold_standard.indels.b37.chr.vcf"

// Scripts
MEDIAN_AWK="/group/bioi1/harrietd/git/micro-genotyper-long/repeat_genotyper_bpipes/median.awk"

// Settings
NSAMPLES=4
EXOME_TARGET="/group/bioi1/harrietd/pooled-parent/data/parentpool/SureSelectClinicalResearchExome_S06588914_Regions.bed"
PLATFORM='illumina'
threads=8 //for BWA

fastqc = {
    doc "Run FASTQC to generate QC metrics for raw reads"
    output.dir = "fastqc"
    from('.fastq.gz')  produce(output.prefix.prefix.prefix + '_fastqc.zip') {
        exec "fastqc -o ${output.dir} $input.gz"
    }
}

def get_info(filename) {
    return(filename.split("/")[-1].split("\\.")[0].split("_"))
}

set_sample_info = {

    doc "Validate and set information about the sample to be processed"

    def info = get_info(input)
    branch.sample = info[0]
    branch.lane = info[1]
    }

align = {
    doc "Extract reads from bam then align with bwa mem algorithm."
    output.dir="align"
    from('fastq.gz', 'fastq.gz') produce(output.prefix.prefix.prefix + '.bam') {
        exec """
            bwa mem -M
            -t $threads
            -R "@RG\\tID:${sample}_${lane}\\tPL:$PLATFORM\\tPU:NA\\tLB:${lane}\\tSM:${sample}"
            $REF $inputs |
            samtools view -bSuh - |
            samtools sort -T $output.prefix -o $output.prefix -
        """, "bwa"
    }
}

@transform('coverage')
coverage = {
    exec """
        bedtools coverage
            -sorted -d
            -g ${REF}.genome
            -a $EXOME_TARGET
            -b $input.bam |
            cut -f 5 |
            sort -n |
            awk -f $MEDIAN_AWK > $output.coverage
            """
}

realignIntervals = {
    doc "Discover candidate regions for realignment in an alignment with GATK"
    output.dir="align"
    exec """
        java -Xmx4g -jar $GATKDIR/GenomeAnalysisTK.jar
            -T RealignerTargetCreator
            -R $REF
            -I $input.bam
            --known $GOLD_STANDARD_INDELS
            -o $output.intervals
            -L $combined_bed
    """, "realign_target_creator"
}

realign = {
    doc "Apply GATK local realignment to specified intervals in an alignment"
    output.dir="align"
    exec """
        java -Xmx5g -jar $GATKDIR/GenomeAnalysisTK.jar
             -T IndelRealigner
             -R $REF
             -I $input.bam
             -targetIntervals $input.intervals
             -o $output.bam
             -L $combined_bed
    ""","local_realign"
}

index_bam = {

    doc "Create an index for a BAM file"

    transform("bam") to ("bam.bai") {
        exec "samtools index $input.bam","index_bam"
    }
    forward input
}

call_variants = {

    output.dir="variants"

    var call_conf:5.0,
    emit_conf:5.0
    branch.ploidy = NSAMPLES*2

     transform("bam","bam") to("metrics","vcf") {
         exec """
                 java -Xmx8g -jar $GATKDIR/GenomeAnalysisTK.jar -T UnifiedGenotyper
                    -R $REF
                    -I $input.bam
                    -nt 4
                    --dbsnp $DBSNP
                    -stand_call_conf $call_conf -stand_emit_conf $emit_conf
                    -dcov 1600
                    -l INFO
                    -L $combined_bed
                    -A AlleleBalance -A Coverage -A FisherStrand
                    -glm BOTH
                    -metrics $output.metrics
                    -o $output.vcf
                    -ploidy $ploidy
             ""","unified_genotyper"
     }
 }

compress_vcf = {
    transform("vcf") to ("vcf.gz") {
        output.dir="variants"
        exec """
            bgzip -f -c $input > $output
        """
    }
}

index_vcf = {
    transform("vcf.gz") to ("vcf.gz.tbi") {
        output.dir="variants"
        exec """
            tabix -f -p vcf $input
        """
        forward input
    }
}

run {
    "%.fastq.gz" * [ fastqc ] +
    '%_R*.fastq.gz' * [
        set_sample_info +
        align + index_bam +
        realignIntervals + realign + index_bam +
        coverage +
        call_variants +
        compress_vcf + index_vcf
    ]
}
