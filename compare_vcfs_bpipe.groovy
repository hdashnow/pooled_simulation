// Bpipe pipeline
// Harriet Dashnow 23 March 2015
// Compares 


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

combined_bed="/vlsci/VR0320/shared/hdashnow/pooled_simulation/analysis/combined_target.bed"

pool4 = ["../vcfs/020430101.merge.dedup.realign.recal.vcf", "../vcfs/020430501.merge.dedup.realign.recal.vcf", "../vcfs/020430901.merge.dedup.realign.recal.vcf", "../vcfs/020431001.merge.dedup.realign.recal.vcf"]
pool6 = ["../vcfs/020430101.merge.dedup.realign.recal.vcf", "../vcfs/020430501.merge.dedup.realign.recal.vcf", "../vcfs/020430901.merge.dedup.realign.recal.vcf", "../vcfs/020431001.merge.dedup.realign.recal.vcf", "../vcfs/020431101.merge.dedup.realign.recal.vcf", "../vcfs/020431201.merge.dedup.realign.recal.vcf"]
pool8 = ["../vcfs/020430101.merge.dedup.realign.recal.vcf", "../vcfs/020430501.merge.dedup.realign.recal.vcf", "../vcfs/020430901.merge.dedup.realign.recal.vcf", "../vcfs/020431001.merge.dedup.realign.recal.vcf", "../vcfs/020431101.merge.dedup.realign.recal.vcf", "../vcfs/020431201.merge.dedup.realign.recal.vcf", "../vcfs/020431401.merge.dedup.realign.recal.vcf", "../vcfs/020431701.merge.dedup.realign.recal.vcf"]
pool10 = ["../vcfs/020430101.merge.dedup.realign.recal.vcf", "../vcfs/020430501.merge.dedup.realign.recal.vcf", "../vcfs/020430901.merge.dedup.realign.recal.vcf", "../vcfs/020431001.merge.dedup.realign.recal.vcf", "../vcfs/020431101.merge.dedup.realign.recal.vcf", "../vcfs/020431201.merge.dedup.realign.recal.vcf", "../vcfs/020431401.merge.dedup.realign.recal.vcf", "../vcfs/020431701.merge.dedup.realign.recal.vcf", "../vcfs/020431801.merge.dedup.realign.recal.vcf", "../vcfs/020431901.merge.dedup.realign.recal.vcf"]

//pool_individuals = [pool4, pool6, pool8, pool10]
//pool_individuals = [pool4]
sample_pools = [2, 4, 6, 8, 10]
//sample_pools = [4]

set_pool = {
    branch.num_samples = branch.name.toInteger()
    def all_inputs = "$inputs".split(" ").toList()
    branch.pool_samples = all_inputs[0..branch.num_samples-1]
    exec """
        echo $num_samples $pool_samples > ${num_samples}.txt
    """
    forward(pool_samples)
}

//@filter('chr22')
filter_vcf = {
    println(output)
    produce(output.prefix.prefix+'.recode.vcf') {
        println(output)
        println("$output.prefix.prefix")
        exec """
            vcftools --vcf $input.vcf --out $output.prefix.prefix --bed $combined_bed --recode
        """
    }
}

compress_vcf = {
    transform("vcf") to ("vcf.gz") {
        exec """
            bgzip -f -c $input > $output
        """
    }
}

index_vcf = {
    transform("vcf.gz") to ("vcf.gz.tbi") {
        exec """
            tabix -f -p vcf $input
        """
        forward input
    }
}

merge_vcfs = {
//    branch.all_inputs = "$inputs".split(" ").toList()
//    branch.num_samples = branch.all_inputs.length()
//    println()
//    println("$inputs.vcf")
    produce('merge.'+branch.num_samples+'.vcf') {
        exec """
            vcf-merge $inputs > $output.vcf
        """
    }
}

run {
    sample_pools * [
        set_pool + 
        "%.vcf" * [
            filter_vcf + compress_vcf + index_vcf 
        ] + merge_vcfs + compress_vcf + index_vcf
    ]
}
