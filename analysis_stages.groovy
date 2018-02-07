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

index_ref = {
    exec "bwa index $REF"
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

        msg "Aligning $inputs size=${inputs.size()}"

        exec """
            bwa mem -M
            -t $threads
            -R "@RG\\tID:${sample}_${lane}\\tPL:$PLATFORM\\tPU:NA\\tLB:${lane}\\tSM:${sample}"
            $REF $inputs |
            samtools view -bSuh - |
            samtools sort -T $output.prefix -o $output.bam -
        """, "bwa"
    }
}

merge_bams = {
    doc """
        Merge the BAM files from multiple lanes together.
        """

    output.dir="align"

    produce(sample + ".merge.bam") {
            msg "Merging $inputs.bam size=${inputs.bam.size()}"
            exec """
                java -Xmx2g -jar $PICARD MergeSamFiles
                    ${inputs.bam.withFlag("INPUT=")}
                    VALIDATION_STRINGENCY=LENIENT
                    ASSUME_SORTED=true
                    CREATE_INDEX=true
                    OUTPUT=$output.bam
             """, "merge"
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
            -L $EXOME_TARGET
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
             -L $EXOME_TARGET
    ""","local_realign"
}

index_bam = {

    doc "Create an index for a BAM file"
    output.dir="align"
    transform("bam") to ("bam.bai") {
        exec "samtools index $input.bam","index_bam"
    }
    forward input
}

call_variants = {

    output.dir="variants"

    branch.ploidy = NSAMPLES*2

     transform("bam","bam") to("metrics","vcf") {
         exec """
                 java -Xmx8g -jar $GATKDIR/GenomeAnalysisTK.jar -T UnifiedGenotyper
                    -R $REF
                    -I $input.bam
                    -nt 4
                    --dbsnp $DBSNP
                    -dcov 1600
                    -l INFO
                    -L $EXOME_TARGET
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


