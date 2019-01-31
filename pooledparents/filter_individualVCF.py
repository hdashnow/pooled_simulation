import argparse
import sys
import vcf
import os
from collections import Counter

__author__ = "Harriet Dashnow"
__credits__ = ["Harriet Dashnow"]
__license__ = "MIT"
__version__ = "0.1.0"
__email__ = "h.dashnow@gmail.com"

def parse_args():
    """Parse the input arguments, use '-h' for help"""
    parser = argparse.ArgumentParser(description='Filter proband vcfs based on pooled partent vcf')
    parser.add_argument(
        '--individual_vcfs', type=str, required=True, nargs='+',
        help='One or more vcfs containing variant calls from individuals used to simulate the pools')
    parser.add_argument(
        '--pool_vcfs', type=str, required=True, nargs='+',
        help='One or more vcfs containing variant calls for the simulated pools')
    parser.add_argument(
        '--pool_specs', type=str, required=True, nargs='+',
        help='One or more text files specifying the bam files used for the simulated pools. Names must correspond to the vcf names in --pool_vcfs')
    parser.add_argument(
        '--output', type=str, required=False,
        help='Output file name. Defaults to stdout.')
    parser.add_argument(
        '--falsepos', type=str, required=False, default='pooled_sim_variants_falsepos.csv',
        help='Output file name for assumed false positives. Defaults to pooled_sim_variants_falsepos.csv.')
    return parser.parse_args()

def variant_id(record):
    ALTstr = '/'.join([str(x) for x in record.ALT]) # join ALT loci
    POSstr = '{0:09d}'.format(record.POS) # add leading 0s
    return '_'.join([str(x) for x in [record.CHROM, POSstr, record.REF, ALTstr]])

def sample_id_from_fname(fname):
    sample_id = os.path.basename(fname).split('.')[0]
    if sample_id == 'merge':
        sample_id = os.path.basename(fname).split('.')[1]
    return(sample_id)

def parse_pool_specs(spec_files):
    """ Expecting file contents in the form:
    [/my/dir/sample1.bam, /my/dir/sample2.bam]
    """
    pool_specs = {}
    for spec_file in spec_files:
        with open(spec_file) as f:
            samples_txt = f.read()
            pool_id = sample_id_from_fname(spec_file)
            samples = []
            for sample_txt in samples_txt.split(', '):
                sample_bam = sample_txt.lstrip().rstrip().lstrip('[').rstrip(']')
                sample_id = sample_id_from_fname(sample_bam)

                samples.append(sample_id)
            pool_specs[pool_id] = samples
    return(pool_specs)

def count_nonref_alleles(GT_string):
    alleles = GT_string.split('/')
    nonref_allele_count = len(alleles) - sum([allele == "0" for allele in alleles])
    return(nonref_allele_count)

def main():
    # Parse command line arguments
    args = parse_args()
    individual_vcf_files = args.individual_vcfs
    pool_vcf_files = args.pool_vcfs
    pool_spec_files = args.pool_specs
    outfile = args.output
    falsepos_file = args.falsepos

    if outfile:
        outstream = open(outfile, 'w')
    else:
        outstream = sys.stdout

    pool_specs = parse_pool_specs(pool_spec_files)
#    print(pool_specs)

    # Parse vcfs on individuals
    individual_vars = {}
    pooled_individual_vars = {} # {pool : { variant_id: non-ref count } }
    for pool in pool_specs:
        pooled_individual_vars[pool] = {}

    for vcf_file in individual_vcf_files:
        sample_id = sample_id_from_fname(vcf_file)
        individual_vars[sample_id] = set()

        pools_sample_is_in = [ pool for pool in pool_specs if sample_id in pool_specs[pool] ]
        
        with open(vcf_file, 'r') as this_vcf:
            vcf_reader = vcf.Reader(this_vcf)
            vcf_writer = vcf.Writer(open(sample_id+".filtered.vcf", 'w'), vcf_reader)
            for record in vcf_reader:

                var_id = variant_id(record)
                individual_vars[sample_id].add(var_id)

                nonref_allele_count = count_nonref_alleles(record.samples[0]['GT'])

                for pool in pools_sample_is_in:
                    if var_id not in pooled_individual_vars[pool]:
                        pooled_individual_vars[pool][var_id] = nonref_allele_count
                    else:
                        pooled_individual_vars[pool][var_id] += nonref_allele_count

                #if var_id not in pool_vars:
                #    vcf_writer.write_record(record)

    # Parse vcfs for pools
    pool_vars = {}
    pool_var_counts = {}
    nonref_allele_counts = {}
    for pool_vcf_file in pool_vcf_files:

        pool = sample_id_from_fname(pool_vcf_file)
        pool_vars[pool] = set()
        pool_var_counts[pool] = {}

        with open(pool_vcf_file, 'r') as this_vcf:
            for record in vcf.Reader(this_vcf):
                var_id = variant_id(record)
                pool_vars[pool].add(var_id)

                nonref_allele_count = count_nonref_alleles(record.samples[0]['GT'])
                nonref_allele_counts[var_id] = nonref_allele_count
   
                if var_id not in pool_var_counts[pool]:
                    pool_var_counts[pool][var_id] = nonref_allele_count
                else:
                    pool_var_counts[pool][var_id] += nonref_allele_count


    header = 'pool,variant,nonref_allele_count_truth,nonref_allele_count_obs,recovered\n'
    outstream.write(header)
    for pool in pooled_individual_vars:
        for variant in pooled_individual_vars[pool]:
            if variant in pool_vars[pool]:
                recovered = "TRUE"
            else:
                recovered = "FALSE"
            nonref_allele_count_truth = pooled_individual_vars[pool][variant]
            try:
                nonref_allele_count_obs = pool_var_counts[pool][variant] 
            except KeyError:
                nonref_allele_count_obs = 'NA'
            outstream.write('{},{},{},{},{}\n'.format(pool,
                variant,
                nonref_allele_count_truth,
                nonref_allele_count_obs,
                recovered))
        
    outstream.close()

    outstream = open(falsepos_file, "w")
    header = 'pool,variant,nonref_allele_count_truth,nonref_allele_count_obs,false_positive\n'
    outstream.write(header)
    for pool in pooled_individual_vars:
        for variant in pool_vars[pool]:
            if variant in pooled_individual_vars[pool]:
                falsepos = "FALSE"
            else:
                falsepos = "TRUE"
            nonref_allele_count_obs = pool_var_counts[pool][variant] 
            try:
                nonref_allele_count_truth = pooled_individual_vars[pool][variant]
            except KeyError:
                nonref_allele_count_truth = 'NA'
            outstream.write('{},{},{},{},{}\n'.format(pool,
                variant,
                nonref_allele_count_truth,
                nonref_allele_count_obs,
                falsepos))
    outstream.close
    # For each pool, check how many variants were found and record their nonref_allele_count
#    for pool in pooled_individual_vars:
#        recovered_vars = [var for var in pooled_individual_vars[pool] if var in pool_vars[pool]]
#        recovered_allele_counts = [pooled_individual_vars[pool][var] for var in recovered_vars]
#        total_allele_counts = [pool_var_counts[pool][var] for var in pool_var_counts[pool]]
#        n_total = len(pool_vars[pool])
#        n_recovered = len(recovered_vars)
#        print()
#        print('Pool:', pool)
#        print(n_total, 'variants')
#        print(n_recovered, 'recovered')
#        print('{:05.2f}% recovered'.format(100*n_recovered/n_total))
#        print('Non-ref allele count for all variants in pool')
#        print(Counter(total_allele_counts))
#        print('Non-ref allele count for recovered variants in pool')
#        print(Counter(recovered_allele_counts))
#
#    for sample_id in individual_vars:
#        print()
#        print('Individual:', sample_id)
#        n_unfiltered = len(individual_vars[sample_id])
#        print(n_unfiltered, 'variants before filtering')
#                   

if __name__ == '__main__':
    main()
