import argparse
import sys
import vcf
from vcf.parser import _Filter
import os
from collections import Counter
import copy
# shared functions
from filtervcf import *

__author__ = "Harriet Dashnow"
__credits__ = ["Harriet Dashnow"]
__license__ = "MIT"
__version__ = "0.1.0"
__email__ = "h.dashnow@gmail.com"

def parse_args():
    """Parse the input arguments, use '-h' for help"""
    parser = argparse.ArgumentParser(description='Filter a multiple-sample vcf of probands and a parent pool based on alleles found in the pool, to find potential de novo variants.')
    parser.add_argument(
        '--in_vcf', type=str, required=True,
        help='A single multi-sample VCF including all probands and the pool, joint called with GATK GenotypeGVCFs')
    parser.add_argument(
        '--pool', type=str, required=False,
        help='Sample name used in the VCF for the pool (assumed to be the last sample if not given)')
    parser.add_argument(
        '--probands', type=str, required=False, nargs='+',
        help='Sample name used in the VCF for the probands (assumed to be all but the last sample if not given)')
    parser.add_argument(
        '--out_csv', type=str, required=False, default='pools_probands.compare.csv',
        help='Output filename for csv (default: %(default)s)')
    parser.add_argument(
        '--out_vcf', type=str, required=False, default='pools_probands.compare.vcf',
        help='Output filename for vcf (default: %(default)s)')
    parser.add_argument(
        '--out_proband_vcf', type=str, required=False, default='.denovo.vcf',
        help='Suffix for output filtered vcfs for each proband, prefix is sample name (default: %(default)s)')
    parser.add_argument(
        '--filter_reads', type=int, required=False,
        help='For the purposes of filtering, consider a variant allele called in the pool if it is supported by this many reads. If not set variants will be filtered only if the variant allele was called in the parent pool genotype.')
    parser.add_argument(
        '--ploidy_filter', type=int, required=False,
        help='The ploidy of the pooled sample (i.e. 2*number of samples in the pool). If provided, allow variants to called in the pool only if they have enough supporting reads based on the ploidy: locus depth * (1/ploidy) * 0.5.')
    parser.add_argument(
        '--falsepos', action='store_true',
        help='Report false positives as additional lines in the output csv.')
    parser.add_argument(
        '--tech_variation', type=float, default=0.5,
        help='Amount of technical variation to allow when choosing an allele frequency threshold based on ploidy. I.e. if 0.5, allow half as many reads to call a variant in the pool.')

    return parser.parse_args()

def main():

    # Parse command line arguments
    args = parse_args()
    vcf_file = args.in_vcf
    outfile = args.out_csv
    out_vcf = args.out_vcf
    report_FPs = args.falsepos
    tech_variation = args.tech_variation

    outstream = open(outfile, 'w')

    # Write header
    outstream.write(('variant,nonref_alleles_pool,total_alleles_pool,'
                    'nonref_alleles_probands,total_alleles_probands,'
                    'nonref_reads_pool,total_reads_pool,nonref_reads_probands,'
                    'recovered_all,falsepos,QD,AF_EXOMESgnomad,AF_GENOMESgnomad,'
                    'proband,recovered_in_proband,GT_pool\n'))

    with open(vcf_file, 'r') as this_vcf:
        vcf_reader = vcf.Reader(this_vcf)
        # Add an aditional filter that will be inherited by the vcf writer
        vcf_reader.filters['InPool'] = _Filter('InPool',
            'All alleles found in the probands are also found in the pool.')
        # Create vcf writer based on the header from the input vcf
        vcf_writer = vcf.Writer(open(out_vcf, 'w'), vcf_reader)

        # Fetch sample names
        all_vcf_samples = vcf_reader.samples
        pool_name = check_pool_name(args.pool, all_vcf_samples)
        proband_names = check_proband_names(args.probands, all_vcf_samples)
        sys.stderr.write('Pool name: {0}\nProband names: {1}\n'.format(pool_name,
                        ', '.join(proband_names)))
        pool_size = len(proband_names)
        pool_pos = all_vcf_samples.index(pool_name)
        probands_pos = [all_vcf_samples.index(proband) for proband in proband_names]

        # Create a vcf writer for each proband
        probandVCF_dict = {}
        for proband in proband_names:
            # Can't deepcopy vcf reader object, so editing it and returning it to previous state
            vcf_reader.samples= [proband]
            proband_out_vcf = proband + args.out_proband_vcf
            probandVCF_dict[proband] = vcf.Writer(open(proband_out_vcf, 'w'), vcf_reader)
        vcf_reader.samples = all_vcf_samples

        for record in vcf_reader:
            # Extract gnomad allele frequency data if available
            AF_EXOMESgnomad = extract_record_info(record, 'AF_EXOMESgnomad')
            AF_GENOMESgnomad = extract_record_info(record, 'AF_GENOMESgnomad')

            var_id = variant_id(record)
            nonref_alleles_pool, total_alleles_pool = count_nonref_alleles(record.samples[pool_pos]['GT'])
            qual = record.QUAL
            QD = qual/record.INFO['DP']

            nonref_reads_pool = count_nonref_reads(record.samples[pool_pos])
            total_reads_pool = record.samples[pool_pos]['DP']

            nonref_reads_probands = 0
            for proband_pos in probands_pos:
                nonref_reads_probands += count_nonref_reads(record.samples[proband_pos])

            GT_pool = record.samples[pool_pos]['GT']
            alleles_in_pool = get_nonref_alleles(record.samples[pool_pos]['GT'])
            alleles_in_probands = set.union(*[get_nonref_alleles(record.samples[pos]['GT']) for pos in probands_pos])

            filtered = 'FALSE'
            falsepos = 'FALSE'
            # Calculate a mininum read filter based on the filter_reads or ploidy_filter
            # arguments, if given
            if args.filter_reads or args.ploidy_filter:
                min_read_filter = set_read_filter(total_reads_pool, args.filter_reads,
                    args.ploidy_filter, tech_variation)
                alleles_in_pool_by_reads = set(alleles_supported(record, pool_pos,
                    min_read_filter, include_ref = False))
                if is_recovered(alleles_in_probands, alleles_in_pool_by_reads):
                    filtered = 'TRUE'
                # likely false positive if found in the pool but not in any of the probands
                if len(alleles_in_pool_by_reads - alleles_in_probands) > 0:
                    falsepos = 'TRUE'

            else:
                # Filter if all the variants found in the probands are also found in the pool
                if is_recovered(alleles_in_probands, alleles_in_pool):
                    filtered = 'TRUE'
                # likely false positive if found in the pool but not in any of the probands
                if len(alleles_in_pool - alleles_in_probands) > 0:
                    falsepos = 'TRUE'

            if filtered == 'TRUE':
                record.FILTER = 'InPool'

            # Count nonref alleles and total alleles in probands
            # Write a filtered vcf for each proband
            nonref_alleles_probands = 0
            total_alleles_probands = 0
            for proband_pos in probands_pos:
                proband = all_vcf_samples[proband_pos]
                nonref, total = count_nonref_alleles(record.samples[proband_pos]['GT'])
                nonref_alleles_probands += nonref
                total_alleles_probands += total

            for proband_pos in probands_pos:
                proband = all_vcf_samples[proband_pos]
                nonref, total = count_nonref_alleles(record.samples[proband_pos]['GT'])

                # Skip variant if this individual has no non-ref alleles e.g. GT is  ./. or 0/0
                if nonref == 0:
                    continue

                # Check if variant is recovered for this proband specifically
                alleles_this_proband = get_nonref_alleles(record.samples[proband_pos]['GT'])

                # Write out the variant (GT for this sample only) to the vcf file for that proband
                # only if the variant is not found in the parent pool
                recovered_proband = 'FALSE'

                if args.filter_reads or args.ploidy_filter:
                    if is_recovered(alleles_this_proband, alleles_in_pool_by_reads):
                        recovered_proband = 'TRUE'
                else:
                    if is_recovered(alleles_this_proband, alleles_in_pool):
                        recovered_proband = 'TRUE'
                if recovered_proband == 'FALSE':
                    tmp_record = copy.deepcopy(record)
                    tmp_record.samples = [record.samples[proband_pos]]
                    probandVCF_dict[proband].write_record(tmp_record)

                outstream.write(','.join([str(x) for x in [var_id,nonref_alleles_pool,
                    total_alleles_pool,nonref_alleles_probands,total_alleles_probands,
                    nonref_reads_pool,total_reads_pool,nonref_reads_probands,filtered,falsepos,QD,
                    AF_EXOMESgnomad, AF_GENOMESgnomad, proband, recovered_proband, GT_pool]]) + '\n')

            # If none of the probands have any non-ref alleles at this locus
            # Still report it in the csv for false positives counts
            if report_FPs and nonref_alleles_probands == 0:
                outstream.write(','.join([str(x) for x in [var_id,nonref_alleles_pool,
                    total_alleles_pool,nonref_alleles_probands,total_alleles_probands,
                    nonref_reads_pool,total_reads_pool,nonref_reads_probands,filtered,falsepos,QD,
                    AF_EXOMESgnomad, AF_GENOMESgnomad, 'NA', 'NA', 'NA']]) + '\n')



            # Write all samples from all variants to VCF
            # includes FILTER InPool for variants where all alleles in probands are recovered in pool
            vcf_writer.write_record(record)

if __name__ == '__main__':
    main()
