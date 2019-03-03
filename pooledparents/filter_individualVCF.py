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
                proband_id = sample_id_from_fname(sample_bam)

                samples.append(proband_id)
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
                    'QD_pool,GT_pool,QD_proband,AF_EXOMESgnomad'
                    '\n'))

    pool_specs = parse_pool_specs(pool_spec_files)

    # Parse vcfs of individuals
    # Variants from individual probands, allocated to the pool to which they belong
    proband_vars_by_pool = {} # {pool : { variant_id: {key:values} } }
    for pool in pool_specs:
        proband_vars_by_pool[pool] = {}

    for vcf_file in individual_vcf_files:
        proband_id = sample_id_from_fname(vcf_file)

        pools_sample_is_in = [ pool for pool in pool_specs if proband_id in pool_specs[pool] ]

        with open(vcf_file, 'r') as this_vcf:
            vcf_reader = vcf.Reader(this_vcf)
            #vcf_writer = vcf.Writer(open(proband_id+".filtered.vcf", 'w'), vcf_reader)
            #XXX add functionality to output filtered vcfs
            for record in vcf_reader:

                variant = variant_id(record)

                nonref_alleles_this_proband, total_alleles_this_proband = count_nonref_alleles(record.samples[0]['GT'])
                nonref_reads_this_proband = count_nonref_reads(record.samples[0])

                AF_EXOMESgnomad = extract_record_info(record, 'AF_EXOMESgnomad')
                qual = record.QUAL
                QD = qual/record.INFO['DP']

                for pool in pools_sample_is_in:
                    if variant not in proband_vars_by_pool[pool]:
                        proband_vars_by_pool[pool][variant] = {
                                                            'nonref_alleles_probands': 0,
                                                            'total_alleles_probands': 0,
                                                            'nonref_reads_probands': 0,
                                                            'probands_with_variant': set(),
                                                            }
                    proband_vars_by_pool[pool][variant]['nonref_alleles_probands'] += nonref_alleles_this_proband
                    proband_vars_by_pool[pool][variant]['total_alleles_probands'] += total_alleles_this_proband
                    proband_vars_by_pool[pool][variant]['nonref_reads_probands'] += nonref_reads_this_proband
                    proband_vars_by_pool[pool][variant]['probands_with_variant'].add(proband_id)
                    proband_vars_by_pool[pool][variant]['AF_EXOMESgnomad'] = AF_EXOMESgnomad
                    proband_vars_by_pool[pool][variant]['QD_proband'] = QD

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

                GT_pool = record.samples[0]['GT']
                nonref_alleles_pool, total_alleles_pool = count_nonref_alleles(GT_pool)
                #nonref_allele_counts[variant] = nonref_allele_pool
                nonref_reads_pool = count_nonref_reads(record.samples[0])
                total_reads_pool = record.samples[0]['DP']
                qual = record.QUAL
                try:
                    QD_pool = qual/record.INFO['DP']
                except KeyError:
                    QD_pool = 'NA'

                if variant not in pool_vars_details[pool]:
                    pool_vars_details[pool][variant] = {
                                                        'nonref_alleles_pool': nonref_alleles_pool,
                                                        'total_alleles_pool': total_alleles_pool,
                                                        'nonref_reads_pool': nonref_reads_pool,
                                                        'total_reads_pool': total_reads_pool,
                                                        'GT_pool': GT_pool,
                                                        'QD_pool': QD_pool
                                                        }
                else:
                    raise Exception('This variant already reported in this pool. Logic error?')

                if args.falsepos: # if reporting false positives
                    # Check if false positive (and only report those)
                    if variant not in proband_vars_by_pool[pool]:
                        recovered = 'FALSE'
                        falsepos = 'TRUE'
                        # Assign some NA values
                        proband = 'NA'
                        recovered_proband = 'NA'
                        QD_proband = 'NA'
                        AF_EXOMESgnomad = 'NA'

                        nonref_alleles_pool = pool_vars_details[pool][variant]['nonref_alleles_pool']
                        nonref_reads_pool = pool_vars_details[pool][variant]['nonref_reads_pool']
                        total_reads_pool = pool_vars_details[pool][variant]['total_reads_pool']
                        GT_pool = pool_vars_details[pool][variant]['GT_pool']
                        QD_pool = pool_vars_details[pool][variant]['QD_pool']

                        if variant in proband_vars_by_pool[pool]:
                            raise Error(("False positive variant from pool found "
                                        "in one of the probands. This doesn't "
                                        "make sense. Logic error in code?"))
                        else:
                            nonref_alleles_probands = 'NA'
                            total_alleles_probands = 'NA'
                            nonref_reads_probands = 'NA'

                        outstream.write(','.join([str(x) for x in [
                            pool, proband, variant, recovered, recovered_proband, falsepos,
                            nonref_alleles_pool, total_alleles_pool,
                            nonref_alleles_probands, total_alleles_probands,
                            nonref_reads_pool, total_reads_pool, nonref_reads_probands,
                            QD_pool,GT_pool,QD_proband,AF_EXOMESgnomad
                            ]]) + '\n')

    for pool in proband_vars_by_pool:
        for variant in proband_vars_by_pool[pool]:

            falsepos = 'FALSE'
            if variant in pool_vars[pool]:
                recovered = "TRUE"
                nonref_alleles_pool = pool_vars_details[pool][variant]['nonref_alleles_pool']
                total_alleles_pool = pool_vars_details[pool][variant]['total_alleles_pool']
                nonref_reads_pool = pool_vars_details[pool][variant]['nonref_reads_pool']
                total_reads_pool = pool_vars_details[pool][variant]['total_reads_pool']
                GT_pool = pool_vars_details[pool][variant]['GT_pool']
                QD_pool = pool_vars_details[pool][variant]['QD_pool']
                recovery_by_proband = which_recovered(proband_vars_by_pool[pool][variant]['probands_with_variant'],
                                        pool_specs[pool])
            else:
                recovered = "FALSE"
                nonref_alleles_pool = 'NA'
                total_alleles_pool = 'NA'
                nonref_reads_pool = 'NA'
                total_reads_pool = 'NA'
                GT_pool = 'NA'
                QD_pool = 'NA'
                # report recovered_proband FALSE for each proband
                recovery_by_proband = which_recovered([],proband_vars_by_pool[pool][variant]['probands_with_variant'])

            nonref_alleles_probands = proband_vars_by_pool[pool][variant]['nonref_alleles_probands']
            total_alleles_probands = proband_vars_by_pool[pool][variant]['total_alleles_probands']
            nonref_reads_probands = proband_vars_by_pool[pool][variant]['nonref_reads_probands']
            QD_proband = proband_vars_by_pool[pool][variant]['QD_proband']
            AF_EXOMESgnomad = proband_vars_by_pool[pool][variant]['AF_EXOMESgnomad']

            for proband in recovery_by_proband:
                recovered_proband = recovery_by_proband[proband]

                outstream.write(','.join([str(x) for x in [
                    pool, proband, variant, recovered, recovered_proband, falsepos,
                    nonref_alleles_pool, total_alleles_pool,
                    nonref_alleles_probands, total_alleles_probands,
                    nonref_reads_pool, total_reads_pool, nonref_reads_probands,
                    QD_pool,GT_pool,QD_proband,AF_EXOMESgnomad
                    ]]) + '\n')

    outstream.close

if __name__ == '__main__':
    main()
