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

    # Write header
    outstream.write(('pool,proband,variant,recovered_proband,recovered,falsepos,'
                    'nonref_alleles_pool,total_alleles_pool,'
                    'nonref_alleles_probands,total_alleles_probands,'
                    'nonref_reads_pool,total_reads_pool,nonref_reads_probands,'
                    'QD,GT_pool'
                    '\n'))

    pool_specs = parse_pool_specs(pool_spec_files)

    # Parse vcfs of individuals
    #individual_vars = {}
    # Variants from individual probands, allocated to the pool to which they belong
    proband_vars_by_pool = {} # {pool : { variant_id: {key:values} } }
    for pool in pool_specs:
        proband_vars_by_pool[pool] = {}

    for vcf_file in individual_vcf_files:
        sample_id = sample_id_from_fname(vcf_file)
        #individual_vars[sample_id] = set()

        pools_sample_is_in = [ pool for pool in pool_specs if sample_id in pool_specs[pool] ]

        with open(vcf_file, 'r') as this_vcf:
            vcf_reader = vcf.Reader(this_vcf)
            vcf_writer = vcf.Writer(open(sample_id+".filtered.vcf", 'w'), vcf_reader)
            for record in vcf_reader:

                variant = variant_id(record)
                #individual_vars[sample_id].add(variant)

                nonref_alleles_probands, total_alleles_probands = count_nonref_alleles(record.samples[0]['GT'])

                for pool in pools_sample_is_in:
                    if variant not in proband_vars_by_pool[pool]:
                        proband_vars_by_pool[pool][variant] = {
                                                            'nonref_alleles_probands': 0,
                                                            'total_alleles_probands': 0,
                                                            }
                    proband_vars_by_pool[pool][variant]['nonref_alleles_probands'] += nonref_alleles_probands
                    proband_vars_by_pool[pool][variant]['total_alleles_probands'] += total_alleles_probands

    # Parse vcfs for pools
    pool_vars = {}
    nonref_allele_counts = {}
    pool_vars_details = {} # {pool : { variant_id: {key:values} } }
    for pool_vcf_file in pool_vcf_files:

        pool = sample_id_from_fname(pool_vcf_file)
        pool_vars[pool] = set()
        pool_vars_details[pool] = {}

        with open(pool_vcf_file, 'r') as this_vcf:
            for record in vcf.Reader(this_vcf):
                variant = variant_id(record)
                pool_vars[pool].add(variant)

                nonref_alleles_pool, total_alleles_pool = count_nonref_alleles(record.samples[0]['GT'])
                #nonref_allele_counts[variant] = nonref_allele_pool

                if variant not in pool_vars_details[pool]: #XXX not sure about this logic, should variant occur multiple times?
                    pool_vars_details[pool][variant] = {
                                                        'nonref_alleles_pool': nonref_alleles_pool,
                                                        'total_alleles_pool': total_alleles_pool,
                                                        }
                else:
                    raise Exception


                if args.falsepos: # if reporting false positives

                    # Assign some placeholder values
                    proband = None
                    recovered_proband = None
                    total_alleles_pool = None
                    #total_alleles_probands = None
                    nonref_reads_pool = None
                    total_reads_pool = None
                    nonref_reads_probands = None
                    QD = None
                    GT_pool = None

                    # Check if false positive (and only report those)
                    if variant not in proband_vars_by_pool[pool]:
                        recovered = 'FALSE'
                        falsepos = 'TRUE'
                        nonref_alleles_pool = pool_vars_details[pool][variant]
                        if variant in proband_vars_by_pool[pool]:
                            nonref_alleles_probands = proband_vars_by_pool[pool][variant]['nonref_alleles_probands']
                            total_alleles_probands = proband_vars_by_pool[pool][variant]['total_alleles_probands']

                        else:
                            nonref_alleles_probands = 'NA'
                            total_alleles_probands = 'NA'

                        outstream.write(','.join([str(x) for x in [
                            pool, proband, variant, recovered, recovered_proband, falsepos,
                            nonref_alleles_pool, total_alleles_pool,
                            nonref_alleles_probands, total_alleles_probands,
                            nonref_reads_pool, total_reads_pool, nonref_reads_probands,
                            QD,GT_pool
                            ]]) + '\n')

    for pool in proband_vars_by_pool:
        for variant in proband_vars_by_pool[pool]:

            # Assign some placeholder values
            proband = None
            recovered_proband = None
            total_alleles_pool = None
            total_alleles_probands = None
            nonref_reads_pool = None
            total_reads_pool = None
            nonref_reads_probands = None
            QD = None
            GT_pool = None

            falsepos = 'FALSE'
            if variant in pool_vars[pool]:
                recovered = "TRUE"
            else:
                recovered = "FALSE"
            nonref_alleles_probands = proband_vars_by_pool[pool][variant]['nonref_alleles_probands']
            total_alleles_probands = proband_vars_by_pool[pool][variant]['total_alleles_probands']

            if variant in pool_vars_details[pool]:
                nonref_alleles_pool = pool_vars_details[pool][variant]['nonref_alleles_pool']
                total_alleles_pool = pool_vars_details[pool][variant]['total_alleles_pool']
            else:
                nonref_alleles_pool = 'NA'
                total_alleles_pool = 'NA'

            outstream.write(','.join([str(x) for x in [
                pool, proband, variant, recovered, recovered_proband, falsepos,
                nonref_alleles_pool, total_alleles_pool,
                nonref_alleles_probands, total_alleles_probands,
                nonref_reads_pool, total_reads_pool, nonref_reads_probands,
                QD,GT_pool
                ]]) + '\n')

    outstream.close

if __name__ == '__main__':
    main()
