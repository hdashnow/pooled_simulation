// Bpipe pipeline
// Harriet Dashnow 23 March 2015
// Compares 

load 'pipeline_config.groovy'

// ***Work in progress, doesn't yet run***

def get_sampleid(filename) {
    String fileContents = new File(filename).text
    println(fileContents.split(', '))
}

// Input is a text file containing input bam files for that pool. 
// Use to get sample IDs to match up with vcfs 
set_pool = {
    branch.num_samples = branch.name.toInteger()
    def all_vcfs = "$inputs.vcf".split(" ").toList()
    branch.pool_sample_ids = get_sampleid(input.txt) //extract read ids from file
    branch.pool_samples = //vcf inputs starting with pool_sample_ids
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
    "%.txt" * [
        set_pool //+ 
//        "%.vcf" * [
//            filter_vcf + compress_vcf + index_vcf 
//        ] + merge_vcfs + compress_vcf + index_vcf
    ]
}
