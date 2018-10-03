// perform GATK4 joint calling on probands and parent pool

load 'pipeline_config.groovy'
load 'pipeline_stages.groovy'

// Settings
// Using the intersection of the exome target regions for the pool and probands if different captures
INT_BED="/group/bioi1/harrietd/pooled-parent/filter/intersected/intersection_Nextera_SureSelectClinical.bed"

/*
sample_file (e.g. samples.txt) describes the probands and pooled samples
type: pool or proband. Expecting a single pool which will be used to filter variants in all probands
ploidy: number of samples in files * 2. For example 2 for probands, 8 for a pool of 4 parents.
sample_id: sample ID or pool name (must be unique within this file)
target: the exome capture target bed file for that sample
files: a comma separated list of all fastq.gz files for that sample, including multiple lanes if any
This is a whitespace separated file. Do not include any whitespace within any of the fields.
*/

sample_file = args[0]

println "Processing samples from $sample_file"

// Parse input file specifying which files belong to which samples and which are pools

// Get the headers of our file
headers = new File(sample_file).readLines()[0].tokenize()

samples =
    new File(sample_file)
        .readLines()[1..-1] // read all the lines except first
        .collect { it.tokenize() } // split each line by whitespace
        .collect { fields ->
            // transpose is like 'zip' in Python, this creates a Map with info about sample
            return [ headers, fields ].transpose().collectEntries()
        }

// Parse fastq files to a list
samples.each {  samplemap ->
    samplemap.files = samplemap.files.tokenize(',')
}

// The above is a List of samples (each sample being a Map / dictionary of info about the sample)
// However, it's nice to be able to look up any sample by its id, so let's index them that way
sample_index = Collections.synchronizedMap(samples.collectEntries { [ it.sample_id, it ] })

// Extract the maps for the probands and pool
probands = samples.grep { it.type == 'proband' }*.sample_id
pool = samples.grep { it.type == 'pool' }*.sample_id
all_samples = samples*.sample_id

println "Probands: " + probands.join(',')
println "Pool: " + pool.join(',')

set_sample_info = {

    // Set a branch variable so that downstream stages could see which sample we are processing
    branch.sample = branch.name
    branch.sample_map = sample_index[branch.sample]
    branch.ploidy = sample_map.ploidy
    branch.EXOME_TARGET = sample_map.target

    // Forward the input file for this sample to later stages
    forward sample_map.files

    println "Processing sample: $branch.sample\nPloidy: $branch.ploidy\nTarget region: $branch.EXOME_TARGET\n"
}


run {
     all_samples * [
        set_sample_info +
        ~'(.*)_R[12].fastq.gz' * [
            set_fastq_info +
            align_bwa + index_bam
        ] +
        merge_lanes +
        dedup +
        coverage +
        call_variants_gvcf
    ] +
    combine_gvcfs +
    joint_calling +
    annotate_vcf +
    intersect_vcf +
    compare_joint
}
