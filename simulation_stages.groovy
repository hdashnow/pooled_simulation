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
    produce(output.prefix.prefix+'.downsampled.'+branch.num_samples+'.bam') {
        output.dir="align"
        DOWNSAMPLE=7 + 1.0/branch.num_samples
        exec """
            samtools view -b -s $DOWNSAMPLE $input.bam > $output.bam 
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

@transform('vcf')
call_variants = {
 
    output.dir="variants"
 
    exec """
        $GATK --java-options "-Xmx48g" HaplotypeCaller 
            -R $REF 
            -I $input.bam 
            --dbsnp $DBSNP 
            -L $combined_bed
            -O $output.vcf
            -ploidy $ploidy
    """, "callvariants"
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

@filter('highqual')
filter_vcf_qual = {

    output.dir="variants"

    exec """
        /group/bioi1/harrietd/src/freebayes/vcflib/bin/vcffilter -f "QUAL > 20" $input.vcf > $output.vcf
    """
}
