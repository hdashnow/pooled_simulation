// Bpipe pipeline
// Harriet Dashnow 22 March 2015
// Takes mapped exome bams. Simulates pooled exomes.


BASE="/group/bioi1/harrietd/pooled-parent/pooled_simulation"

// Tools
TOOLS="$BASE/tools"
PICARD="/usr/local/installed/picard/2.0.1/picard.jar"
GATK="/usr/local/installed/gatk/3.8"

// Ref files
EXOME_TARGET="/group/bioi1/shared/genomes/hg19/exome_targets/nextera_rapid_capture_exome/target_region.bed"
EXCLUDE="$BASE/CS_excluded.bed"
REFBASE="/group/bioi1/shared/genomes/hg19/gatk" 
REF="$REFBASE/gatk.ucsc.hg19.fasta"
DBSNP="$REFBASE/dbsnp_138.hg19.vcf"
GOLD_STANDARD_INDELS="$REFBASE/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf"

combined_bed="combined_target.bed"

// Settings
CHR='chr22'
sample_pools = [2, 4, 6, 8, 10]
set_pool = {
    branch.num_samples = branch.name.toInteger()
    branch.ploidy = branch.num_samples*2
    def all_inputs = "$inputs".split(" ").toList()

    Collections.shuffle(all_inputs)
    branch.pool_samples = all_inputs[0..branch.num_samples-1]

    produce(branch.num_samples + ".txt") {
    exec """
        echo "$pool_samples" > $output.txt
    """
    // need to check if there are enough input sequences to create pool of this size
    forward(pool_samples)
}
}

//  -v  Only report those entries in A that have _no overlaps_ with B.
//        - Similar to "grep -v" (an homage).
intersect_targets = {
    produce(combined_bed) {
        exec """
           bedtools intersect -v -a $EXOME_TARGET -b $EXCLUDE  > $output.bed
        """
    forward inputs
    }    
}

//@filter('downsampled')
downsample_region = {
    produce(output.prefix.prefix+'.downsampled.'+branch.num_samples+'.bam') {
        output.dir="align"
        DOWNSAMPLE=7 + 1.0/branch.num_samples
        exec """
            samtools view -b -s $DOWNSAMPLE $input.bam > $output 
        """
    }
}

merge_bams = {
    produce('merge.'+branch.num_samples+'.bam') {
        output.dir="align"
        exec """
            java -Xmx4g -jar $PICARD MergeSamFiles
                ${inputs.bam.withFlag('I=')}
                O=$output.bam
        """
    }
}

@filter('RGfixed')
fix_header = {
    output.dir="align"
    exec """
        java -Xmx4g -jar $PICARD AddOrReplaceReadGroups
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

@preserve("*.bam")
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
        exec "samtools index $input.bam"
    }
    forward input
}

call_variants = {
 
    output.dir="variants"
 
     transform("bam","bam") to("metrics","vcf") {
         exec """
                 java -Xmx8g -jar $GATK/GenomeAnalysisTK.jar -T UnifiedGenotyper 
                    -R $REF 
                    -I $input.bam 
                    -nt 4
                    --dbsnp $DBSNP 
                    -dcov 1600 
                    -l INFO 
                    -L $combined_bed
                    -A AlleleBalance -A Coverage -A FisherStrand 
                    -glm BOTH
                    -metrics $output.metrics
                    -o $output.vcf
                    -ploidy $ploidy
             """, "unified_genotyper"
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

cleanup = {
    cleanup "*.bam", "*.vcf", "*.intervals"
}

run {
    intersect_targets +
    sample_pools * [
        set_pool + "%.bam" * [
            downsample_region
        ] + merge_bams + fix_header +
        realignIntervals + realign + index_bam + 
        call_variants +
        compress_vcf + index_vcf //+ 
//        cleanup
    ]
}
