load 'pipeline_stages.groovy'

set_pool = {
    branch.num_samples = branch.name.split("\\.")[1].toInteger()
    branch.ploidy = branch.num_samples*2
}


@transform("vcf")
call_variants_freebayes = {
    exec """
        freebayes 
            --ploidy $ploidy 
            --pooled-discrete 
            --use-best-n-alleles 4 
            --fasta-reference /group/bioi1/shared/genomes/hg19/gatk/gatk.ucsc.hg19.fasta 
            $input.bam > $output.vcf
    ""","freebayes"
}

run {
    "%.bam" * [ set_pool + call_variants_freebayes + filter_vcf_qual ] + compare_sim_freebayes
}
