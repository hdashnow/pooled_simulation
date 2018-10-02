// Bpipe pipeline
// Harriet Dashnow 22 March 2015
// Takes mapped exome bams. Simulates pooled exomes.

///////////////////
// Helper functions

def get_fname(path) {
    x = path.split('/')[-1]
    return(x)
}

def get_info(filename) {
    return(filename.split('/')[-1].split('\\.')[0].split('_'))
}

set_sample_info = {

    doc "Extract sample name from file (up to the first underscore)"

    def info = get_info(input)
    branch.sample = info[0]
}

set_fastq_info = {
    def info = get_info(input)

    branch.sample = info[0]

    if (info.length >= 2) {
        branch.lane = info[-2]
    } else {
        branch.lane = 'L001'
    }

    if (info.length >= 5) {
        branch.flowcell = info[-4]
        branch.library = info[-3]
        branch.lane = info[-2]
    } else {
        branch.flowcell = 'NA'
        branch.library = 'NA'
        branch.lane = 'L001'
    }

}

fastqc = {
    doc "Run FASTQC to generate QC metrics for raw reads"
    output.dir = "fastqc"
    from('.fastq.gz')  produce(output.prefix.prefix.prefix + '_fastqc.zip') {
        exec "fastqc -o ${output.dir} $input.gz"
    }
}

@preserve("*.bam")
align_bwa = {
    doc "Align reads with bwa mem algorithm."

    output.dir="align"

    def fastaname = get_fname(REF)
            //set -o pipefail
    from('fastq.gz', 'fastq.gz') produce(branch.name + '.bam') {
        exec """
            bwa mem -M -t $threads
            -R "@RG\\tID:${flowcell}.${lane}\\tPL:$PLATFORM\\tPU:${flowcell}.${lane}.${library}\\tLB:${library}\\tSM:${sample}"
            $REF
            $inputs |
            samtools view -bSuh - | samtools sort -o $output.bam -T $output.bam.prefix
        """, "bwa"
    }
}

@preserve("*.bai")
index_bam = {
    
    output.dir="align"
    
    transform("bam") to("bam.bai") {
        exec "samtools index $input.bam"
    }
    forward input
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

set_pool = {
    branch.num_samples = branch.name.toInteger()
    branch.ploidy = branch.num_samples*2
    def all_inputs = "$inputs".split(" ").toList()
    
    var seed : false
    if(!seed)
        seed=0

    branch.randomseed = branch.num_samples + seed

    Collections.shuffle(all_inputs, new Random(branch.randomseed))
    branch.pool_samples = all_inputs[0..branch.num_samples-1]

    produce(branch.num_samples + ".txt") {
    exec """
        echo "$pool_samples" > $output.txt
    ""","small"
    // need to check if there are enough input sequences to create pool of this size
    forward(pool_samples)
}
}


//@filter('downsampled')
downsample_region = {

    doc """Downsample bam files. Calculate the proportion of reads to select
        by taking 2/pool size so that the final pool will have on average 2x
        the number of reads in a single sample. If only pooling 2 samples,
        the two bam files will simply be combined."""

    output.dir="align"

    produce(output.prefix.prefix+'.downsampled.'+branch.num_samples+'.bam') {

        if(branch.num_samples > 2) {
            // the integer part (7) is the random seed
            // fractional part is the proprtion of reads to sample
            DOWNSAMPLE=7 + 2.0/branch.num_samples
            exec """
                samtools view -b -s $DOWNSAMPLE $input.bam > $output.bam
            """
        } else if (branch.num_samples == 2) {
            exec """
                cp $input.bam $output.bam
            """
        } else {
            fail "Tried to pool $num_samples samples. This isn't supported."
        }
    }
}

merge_bams = {

    output.dir="align"

    produce('merge.'+branch.num_samples+'.bam') {
        msg "Merging $inputs.bam size=${inputs.bam.size()}"
        exec """
            java -Xmx8g -jar $PICARD MergeSamFiles
                ${inputs.bam.withFlag("INPUT=")}
                OUTPUT=$output.bam
                VALIDATION_STRINGENCY=LENIENT
                ASSUME_SORTED=true
                CREATE_INDEX=true
        ""","merge"
    }
}

merge_lanes = {
    doc """
        Merge the BAM files from multiple lanes together.
        """

    output.dir="align"

    produce(branch.sample + ".merge.bam") {
        msg "Merging $inputs.bam size=${inputs.bam.size()}"
        exec """
            java -Xmx8g -jar $PICARD MergeSamFiles
                ${inputs.bam.withFlag("INPUT=")}
                OUTPUT=$output.bam
                VALIDATION_STRINGENCY=LENIENT
                ASSUME_SORTED=true
                CREATE_INDEX=true
         """, "merge"
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

dedup = {
    doc "Remove PCR duplicates from reads"
    output.dir="align"

    def safe_tmp_dir = [TMPDIR, UUID.randomUUID().toString()].join( File.separator )
    exec """
        mkdir -p "$safe_tmp_dir"

        java -Xmx4g -Djava.io.tmpdir=$safe_tmp_dir -jar $PICARD MarkDuplicates
             INPUT=$input.bam 
             REMOVE_DUPLICATES=true 
             VALIDATION_STRINGENCY=LENIENT 
             AS=true 
             METRICS_FILE=$output.metrics
             CREATE_INDEX=true
             OUTPUT=$output.bam

        rm -r "$safe_tmp_dir"
    """
}

@transform('vcf')
call_variants = {
 
    output.dir="variants"
 
    exec """
        $GATK --java-options "-Xmx48g" HaplotypeCaller 
            -R $REF 
            -I $input.bam 
            --dbsnp $DBSNP 
            -L $EXOME_TARGET
            -O $output.vcf
            -ploidy $ploidy
    """, "callvariants"
}

compress_vcf = {

    output.dir="variants"

    transform("vcf") to ("vcf.gz") {
        exec """
            bgzip -f -c $input > $output
        """
    }
}

index_vcf = {

    output.dir="variants"

    transform("vcf.gz") to ("vcf.gz.tbi") {
        exec """
            tabix -f -p vcf $input
        """
        forward input
    }
}

cleanup = {
    cleanup "*.bam", "*.vcf", "*.intervals"
}

@filter('highqual')
filter_vcf_qual = {

    output.dir="variants"

    exec """
        /group/bioi1/harrietd/src/freebayes/vcflib/bin/vcffilter -f "QUAL > 20" $input.vcf > $output.vcf
    """
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

@filter('annotated')
annotate_vcf = {
    output.dir="variants"
    exec """
        set -o pipefail

        vcfanno -base-path $GNOMAD $VCFANNO_CONFIG $input.vcf > $output.vcf
    """
}

@filter('intersect')
intersect_vcf = {
    doc "Filter vcfs to the intersection of the target regions of the two exome captures used"
    output.dir="variants"
    exec """
        bedtools intersect -header -a $input.vcf -b $INT_BED > $output.vcf
    """
}

compare_sim = {

    output.dir="variants"

    from('*.RGfixed.vcf') produce('pooled_sim_compare.csv', 'pooled_sim_compare_falsepos.csv') {

        exec """
            /group/bioi1/harrietd/git/STRetch/tools/bin/python /group/bioi1/harrietd/git/pooled_simulation/compare_sim_vcf.py --individual_vcfs /group/bioi1/harrietd/pooled-parent/pooled_simulation2/simplex/individuals/variants/SRR???????.vcf --pool_vcfs $inputs.vcf --pool_specs $inputs.txt --output $output1.csv --falsepos $output2.csv
    """
    }
}
