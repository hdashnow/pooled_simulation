load 'pipeline_stages.groovy'

@transform("vcf")
call_variants_freebayes = {
    exec """
        freebayes --fasta-reference /group/bioi1/shared/genomes/hg19/gatk/gatk.ucsc.hg19.fasta $input.bam > $output.vcf
    """, "freebayes"
}

run {
    "%.bam" * [ call_variants_freebayes + filter_vcf_qual ]
}
