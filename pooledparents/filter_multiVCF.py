import argparse
import sys
import vcf
from vcf.parser import _Filter
import os
from collections import Counter
import copy
import math

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

    return parser.parse_args()

def variant_id(record):
    """Create a unique ID for each variant so they can be compared"""
    ALTstr = '/'.join([str(x) for x in record.ALT]) # join ALT loci
    POSstr = '{0:09d}'.format(record.POS) # add leading 0s
    return '_'.join([str(x) for x in [record.CHROM, POSstr, record.REF, ALTstr]])

def sample_id_from_fname(fname):
    """Extract same id from filename"""
    sample_id = os.path.basename(fname).split('.')[0]
    if sample_id == 'merge':
        sample_id = os.path.basename(fname).split('.')[1]
    return(sample_id)

def count_nonref_alleles(GT_string):
    """Count the number of non-reference alleles from a VCF genotype (GT) string"""
    alleles = GT_string.split('/')
    alleles = [allele for allele in alleles if allele != '.'] # remove missing genotypes
    total_alleles = len(alleles)
    nonref_allele_count = total_alleles - sum([allele == "0" for allele in alleles])
    return(nonref_allele_count, total_alleles)

def count_nonref_reads(record_sample):
    """Count the number of reads supporting all non-reference alleles"""
    allelic_depths = record_sample['AD']
    return(sum(allelic_depths[1:]))

def alleles_supported(record, sample_pos, n, include_ref = True):
    """Return alleles supported by at least n reads
    Args:
        record (pyvcf record object)
        sample_pos (int): position of the desired sample in the vcf (0-based)
        n (int): return allele if it is supported by at least n reads
        include_ref (bool): if True, return the reference allele if it is supported by n reads
    Returns:
        list of alleles that are supported by at least n reads (strings)
    """
    # Get numeric representation of the alleles i.e. ref = 0, first alt = 1
    all_alleles = [str(x) for x in range(len(record.alleles))]
    allelic_depths = record.samples[sample_pos]['AD']
    if not include_ref:
        all_alleles = all_alleles[1:]
        allelic_depths = allelic_depths[1:]
    supported_alleles = [all_alleles[i] for i in range(len(all_alleles)) if allelic_depths[i] >= n]
    return(supported_alleles)

def get_nonref_alleles(GT_string):
    """Take a VCF genotype (GT) string and return a set containing all non-reference alleles"""
    alleles = set(GT_string.split('/'))
    try:
        alleles.remove('.') # remove missing genotypes
    except KeyError:
        pass
    try:
        alleles.remove('0') # remove reference alleles
    except KeyError:
        pass
    return(alleles)

def is_recovered(alleles_in_probands, alleles_in_pool):
    """Filter if all the variants found in the proband(s) are also found in the pool.
    This tends to result in multi-alleleic sites not getting filtered in many cases.

    alleles_in_probands, alleles_in_pool: iterable consisting of items that can be compared in a set
    """
    if len(set(alleles_in_probands) - set(alleles_in_pool)) == 0:
        return True
    else:
        return False

def set_read_filter(total_reads_pool, filter_reads = None, ploidy = None,
                    tech_variation = 0.5, combine_filters = 'strict'):
    """Set a minimum read depth filter based on either an input value, the ploidy,
        or a combination of the two.
    Args:
        total_reads_pool (int): read depth at this position
        filter_reads (int): filter variants with less than this many reads
        ploidy (int): perform variant frequency filter based on this ploidy
        tech_variation (float): 0 > tech_variation >= 1, ploidy filter will be multiplied
            by this values
        combine_filters (str): 'strict': take the maximum if two filters are used, 
            'lenient': take the minimum if two filters are used
    Returns:
        int
    """
    if filter_reads or ploidy:
        min_read_filter_reads = min_read_filter_ploidy = 0
        if filter_reads:
            min_read_filter_reads = filter_reads
        if ploidy:
            min_read_filter_ploidy = math.trunc(total_reads_pool * (1/ploidy) * tech_variation)
    if filter_reads and ploidy:
        if combine_filters == 'strict':
            min_read_filter = max(min_read_filter_reads, min_read_filter_ploidy)
        elif combine_filters == 'lenient':
            min_read_filter = min(min_read_filter_reads, min_read_filter_ploidy)
        else:
            raise ValueError("Valid values of combine_filters are 'strict' or 'lenient'. The value '{}' was given.".format(combine_filters))
    elif filter_reads:
        min_read_filter = min_read_filter_reads
    elif ploidy:
        min_read_filter = min_read_filter_ploidy
    if min_read_filter < 1:
        min_read_filter = 1

    return(min_read_filter)

def main():
    # Parse command line arguments
    args = parse_args()
    vcf_file = args.in_vcf
    outfile = args.out_csv
    out_vcf = args.out_vcf

    outstream = open(outfile, 'w')

    # Amount of technical variation to allow when choosing an allele frequency threshold
    # based on ploidy. I.e. if 0.5, allow half as many reads to call a variant in the pool.
    tech_variation = 0.5

    outstream.write('variant,nonref_alleles_pool,total_alleles_pool,nonref_alleles_probands,total_alleles_probands,nonref_reads_pool,total_reads_pool,recovered_all,falsepos,QD,AF_EXOMESgnomad,AF_GENOMESgnomad,proband,recovered_in_proband,GT_pool\n')

    with open(vcf_file, 'r') as this_vcf:
        vcf_reader = vcf.Reader(this_vcf)
        # Add an aditional filter that will be inherited by the vcf writer
        vcf_reader.filters['InPool'] = _Filter('InPool',
            'All alleles found in the probands are also found in the pool.')
        # Create vcf writer based on the header from the input vcf
        vcf_writer = vcf.Writer(open(out_vcf, 'w'), vcf_reader)

        all_vcf_samples = vcf_reader.samples
        # Fetch sample names
        if args.pool:
            pool_name = args.pool
            if pool_name not in all_vcf_samples:
                sys.exit('ERROR: Specified pool name "{}" was not found in the VCF. Samples in the VCF include: {}'.format(pool_name, ', '.join(all_vcf_samples)))
        else:
            sys.stderr.write('WARNING: Pool name not given, inferred to be last sample in VCF.\n')
            pool_name = all_vcf_samples[-1]
        if args.probands:
            proband_names = args.probands
            for proband in proband_names:
                if proband not in all_vcf_samples:
                    sys.exit('ERROR: Specified proband name "{}" was not found in the VCF. Samples in the VCF include: {}'.format(proband, ', '.join(all_vcf_samples)))
        else:
            sys.stderr.write('WARNING: Proband names not given, inferred to be all but the last sample in VCF.\n')
            proband_names = all_vcf_samples[:-1]
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
            try:
                AF_EXOMESgnomad = record.INFO['AF_EXOMESgnomad']
                if len(AF_EXOMESgnomad) == 1:
                    AF_EXOMESgnomad = AF_EXOMESgnomad[0]
                else:
                    AF_EXOMESgnomad = 'NA' # throw away gnomad data at multiallelic sites
            except KeyError:
                AF_EXOMESgnomad = 'NA'
            try:
                AF_GENOMESgnomad = record.INFO['AF_GENOMESgnomad']
                if len(AF_GENOMESgnomad) == 1:
                    AF_GENOMESgnomad = AF_GENOMESgnomad[0]
                else:
                    AF_GENOMESgnomad = 'NA'
            except KeyError:
                AF_GENOMESgnomad = 'NA'

            var_id = variant_id(record)
            nonref_alleles_pool, total_alleles_pool = count_nonref_alleles(record.samples[pool_pos]['GT'])
            qual = record.QUAL
            QD = qual/record.INFO['DP']

            nonref_reads_pool = count_nonref_reads(record.samples[pool_pos])
            total_reads_pool = record.samples[pool_pos]['DP']

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
                
                if args.filter_reads:
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
                    nonref_reads_pool,total_reads_pool,filtered,falsepos,QD,
                    AF_EXOMESgnomad, AF_GENOMESgnomad, proband, recovered_proband, GT_pool]]) + '\n')

            # If none of the probands have any non-ref alleles at this locus
            # Still report it in the csv for false positives counts
            if nonref_alleles_probands == 0:
                outstream.write(','.join([str(x) for x in [var_id,nonref_alleles_pool,
                    total_alleles_pool,nonref_alleles_probands,total_alleles_probands,
                    nonref_reads_pool,total_reads_pool,filtered,falsepos,QD,
                    AF_EXOMESgnomad, AF_GENOMESgnomad, 'NA', 'NA', 'NA']]) + '\n')



            # Write all samples from all variants to VCF
            # includes FILTER InPool for variants where all alleles in probands are recovered in pool
            vcf_writer.write_record(record)

#            outstream.write(','.join([str(x) for x in [var_id,nonref_alleles_pool,
#                total_alleles_pool,nonref_alleles_probands,total_alleles_probands,
#                nonref_reads_pool,total_reads_pool,filtered,falsepos,QD,
#                AF_EXOMESgnomad, AF_GENOMESgnomad]]) + '\n')

if __name__ == '__main__':
    main()
