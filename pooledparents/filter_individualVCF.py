import argparse
import sys
import vcf
import os
from collections import Counter
# shared functions
from filtervcf import *

__author__ = "Harriet Dashnow"
__credits__ = ["Harriet Dashnow"]
__license__ = "MIT"
__version__ = "0.1.0"
__email__ = "h.dashnow@gmail.com"

def parse_args():
    """Parse the input arguments, use '-h' for help"""
    parser = argparse.ArgumentParser(description='Filter proband vcfs based on a pooled parent vcf and report recall per variant')
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
        '--out_csv', type=str, required=False, default='pools_probands.compare.csv',
        help='Output filename for csv (default: %(default)s)')
    parser.add_argument(
        '--falsepos', action='store_true',
        help='Report false positives as additional lines in the output csv.')

    return parser.parse_args()

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

def main():
    # Parse command line arguments
    args = parse_args()
    individual_vcf_files = args.individual_vcfs
    pool_vcf_files = args.pool_vcfs
    pool_spec_files = args.pool_specs
    outfile = args.out_csv

    outstream = open(outfile, 'w')

    # outstream.write(('variant,nonref_alleles_pool,total_alleles_pool,'
    #                 'nonref_alleles_probands,total_alleles_probands,'
    #                 'nonref_reads_pool,total_reads_pool,nonref_reads_probands,'
    #                 'recovered_all,falsepos,QD,AF_EXOMESgnomad,AF_GENOMESgnomad,'
    #                 'proband,recovered_in_proband,GT_pool\n'))

    header = ('pool,variant,nonref_allele_count_truth,nonref_allele_count_obs,'
            'recovered,falsepos\n')
    outstream.write(header)

    pool_specs = parse_pool_specs(pool_spec_files)

    # Parse vcfs of individuals
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

                variant = variant_id(record)
                individual_vars[sample_id].add(variant)

                nonref_allele_count, _total_alleles = count_nonref_alleles(record.samples[0]['GT'])

                for pool in pools_sample_is_in:
                    if variant not in pooled_individual_vars[pool]:
                        pooled_individual_vars[pool][variant] = nonref_allele_count
                    else:
                        pooled_individual_vars[pool][variant] += nonref_allele_count

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
                variant = variant_id(record)
                pool_vars[pool].add(variant)

                nonref_allele_count, _total_alleles = count_nonref_alleles(record.samples[0]['GT'])
                nonref_allele_counts[variant] = nonref_allele_count

                if variant not in pool_var_counts[pool]:
                    pool_var_counts[pool][variant] = nonref_allele_count
                else:
                    pool_var_counts[pool][variant] += nonref_allele_count

                if args.falsepos: # if reporting false positives
                    # Check if false positive (and only report those)
                    if variant not in pooled_individual_vars[pool]:
                        recovered = 'FALSE'
                        falsepos = 'TRUE'
                        nonref_allele_count_obs = pool_var_counts[pool][variant]
                        try:
                            nonref_allele_count_truth = pooled_individual_vars[pool][variant]
                        except KeyError:
                            nonref_allele_count_truth = 'NA'
                        outstream.write('{},{},{},{},{},{}\n'.format(
                        pool,
                        variant,
                        nonref_allele_count_truth,
                        nonref_allele_count_obs,
                        recovered,
                        falsepos))

    for pool in pooled_individual_vars:
        for variant in pooled_individual_vars[pool]:
            falsepos = 'FALSE'
            if variant in pool_vars[pool]:
                recovered = "TRUE"
            else:
                recovered = "FALSE"
            nonref_allele_count_truth = pooled_individual_vars[pool][variant]
            try:
                nonref_allele_count_obs = pool_var_counts[pool][variant]
            except KeyError:
                nonref_allele_count_obs = 'NA'
            outstream.write('{},{},{},{},{},{}\n'.format(
                pool,
                variant,
                nonref_allele_count_truth,
                nonref_allele_count_obs,
                recovered,
                falsepos))

    outstream.close

if __name__ == '__main__':
    main()
