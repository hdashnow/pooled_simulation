
/**
 * This example shows how to construct a pipeline that reads the data
 * to process from a file.
 *
 * To run this pipeline, use:
 * 
 *     bpipe run pipeline.groovy  samples.txt 
 */


/*
This file describes the probands and pooled samples
type: pool or proband. Expecting a single pool which will be used to filter variants in all probands
ploidy: number of samples in files * 2. For example 2 for probands, 8 for a pool of 4 parents.
sample_id: sample ID or pool name (must be unique within this file)
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

probands = samples.grep { it.type == 'proband' }*.sample_id
pool = samples.grep { it.type == 'pool' }*.sample_id

println "Probands: " + probands.join(',')
println "Pool: " + pool.join(',')

set_proband_info = {

    // Set a branch variable so that downstream stages could see which sample we are processing
    branch.sample = sample_index[branch.name]
    branch.ploidy = sample.ploidy

    // Forward the input file for this sample to later stages
    forward sample.files

    println "Processing proband: $branch.name with ploidy: $branch.ploidy"
    println "$sample.files"
}

run {
     probands * [ 
        set_proband_info
    ]
}
