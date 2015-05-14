// Bpipe pipeline
// Harriet Dashnow 22 March 2015
// Takes mapped exome bams. Simulates pooled exomes.


BASE="/vlsci/VR0320/shared/hdashnow/MelbGenomicsPipelineRepo"

// Tools
TOOLS="$BASE/tools"
BEDTOOLS="$TOOLS/bedtools/2.18.2"
PICARD="$TOOLS/picard/picard-tools-1.65/lib"
SAMTOOLS="$TOOLS/samtools/0.1.19"
GATK="$TOOLS/gatk/2.8-1-g932cd3a"

// Ref files
EXOME_TARGET="$BASE/designs/nextera_rapid_capture_exome_1.2/target_regions.bed"
EXCLUDE='/vlsci/VR0320/shared/hdashnow/pooled_simulation/CS_excluded.bed'
REFBASE="$BASE/hg19" 
REF="$REFBASE/ucsc.hg19.fasta"
DBSNP="$REFBASE/dbsnp_138.hg19.vcf"
GOLD_STANDARD_INDELS="$REFBASE/Mills_and_1000G_gold_standard.indels.hg19.vcf"

combined_bed="combined_target.bed"

// Settings
//SAMPLES=4
//DOWNSAMPLE=7 + 1.0/SAMPLES
CHR='chr22'
sample_pools = [2, 4, 6, 8, 10]
set_pool = {
//    println(inputs)
    branch.num_samples = branch.name.toInteger()
    branch.ploidy = branch.num_samples*2
    def all_inputs = "$inputs".split(" ").toList()
//    println(branch.num_samples)
    branch.pool_samples = all_inputs[0..branch.num_samples-1]
    exec """
        echo $num_samples $pool_samples > ${num_samples}.txt
    """
    // need to check if there are enough input sequences to create pool of this size
    forward(pool_samples)
}

//  -v  Only report those entries in A that have _no overlaps_ with B.
//        - Similar to "grep -v" (an homage).
intersect_targets = {
    produce('combined_target.bed') {
        exec """
           $BEDTOOLS/bin/bedtools intersect -v -a $EXOME_TARGET -b $EXCLUDE  > $output.bed
        """
//    forward inputs
    }    
}

//@filter('downsampled')
downsample_region = {
    produce(output.prefix+'.downsampled'+branch.num_samples+'.bam') {
        output.dir="align"
        DOWNSAMPLE=7 + 1.0/branch.num_samples
        exec """
            $SAMTOOLS/samtools view -b -s $DOWNSAMPLE $input.bam > $output 
        """
    }
}

//@filter('merged')
//merge_bams = {
//    produce(CHR+'.merge.'+branch.num_samples+'.bam') {
//        exec """
//            samtools merge -h $input.bam $output.bam $inputs.bam
//        """
//    }
//}

merge_bams = {
    produce('merge.'+branch.num_samples+'.bam') {
        output.dir="align"
        exec """
            java -Xmx4g -jar $PICARD/MergeSamFiles.jar
                ${inputs.bam.withFlag('I=')}
                O=$output.bam
        """
    }
}

@filter('RGfixed')
fix_header = {
    output.dir="align"
    exec """
        java -Xmx4g -jar $PICARD/AddOrReplaceReadGroups.jar
            I=$input.bam
            O=$output.bam
            SORT_ORDER=coordinate
            RGPU=na
            RGID=1
            RGLB=input
            RGPL=Illumina
            RGSM=Company
            CREATE_INDEX=True
    """
}

realignIntervals = {
    doc "Discover candidate regions for realignment in an alignment with GATK"
    output.dir="align"
    exec """
        java -Xmx4g -jar $GATK/GenomeAnalysisTK.jar 
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
        java -Xmx5g -jar $GATK/GenomeAnalysisTK.jar 
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

    // A bit of a hack to ensure the index appears in the
    // same directory as the input bam, no matter where it is
    // nb: fixed in new version of Bpipe
    output.dir=file(input.bam).absoluteFile.parentFile.absolutePath
    transform("bam") to ("bam.bai") {
        exec "$SAMTOOLS/samtools index $input.bam"
    }
    forward input
}

call_variants = {
 
    output.dir="variants"
     var call_conf:5.0,
         emit_conf:5.0
 
     transform("bam","bam") to("metrics","vcf") {
         exec """
                 java -Xmx8g -jar $GATK/GenomeAnalysisTK.jar -T UnifiedGenotyper 
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
             """
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
//    intersect_targets +
//    '%.bam' * [
    sample_pools * [
        set_pool + "%.bam" * [
            downsample_region
        ] + merge_bams + fix_header +
        realignIntervals + realign + index_bam + 
        call_variants +
        compress_vcf + index_vcf
    ]
}
