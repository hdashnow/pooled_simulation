import argparse
import sys
import vcf
from vcf.parser import _Filter
import os
from collections import Counter

__author__ = "Harriet Dashnow"
__credits__ = ["Harriet Dashnow"]
__license__ = "MIT"
__version__ = "0.1.0"
__email__ = "h.dashnow@gmail.com"

def parse_args():
    """Parse the input arguments, use '-h' for help"""
    parser = argparse.ArgumentParser(description='Filter a multiple-sample vcf of probands and a parent pool based on alleles found in the pool.')
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
        '--filter_reads', type=int, required=False,
        help='Filter variants where this number of variant reads is observed in the parent pool. If not set variants will be filtered only if the variant allele was called in the parent pool.')
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

def main():
    # Parse command line arguments
    args = parse_args()
    vcf_file = args.in_vcf
    outfile = args.out_csv
    out_vcf = args.out_vcf

    outstream = open(outfile, 'w')

    # Write header
    outstream.write('variant,nonref_alleles_pool,total_alleles_pool,nonref_alleles_probands,total_alleles_probands,nonref_reads_pool,total_reads_pool,recovered,falsepos,QUAL\n')

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
        
        for record in vcf_reader:
            var_id = variant_id(record)
            nonref_alleles_pool, total_alleles_pool = count_nonref_alleles(record.samples[pool_pos]['GT'])
            qual = record.QUAL

            nonref_alleles_probands = 0
            total_alleles_probands = 0
            for proband_pos in probands_pos:
                nonref, total = count_nonref_alleles(record.samples[proband_pos]['GT'])
                nonref_alleles_probands += nonref
                total_alleles_probands += total
            
            nonref_reads_pool = count_nonref_reads(record.samples[pool_pos])
            total_reads_pool = record.samples[pool_pos]['DP']

            alleles_in_pool = get_nonref_alleles(record.samples[pool_pos]['GT'])
            alleles_in_probands = set.union(*[get_nonref_alleles(record.samples[pos]['GT']) for pos in probands_pos])

            filtered = 'FALSE'
            if args.filter_reads:
                # Filter if any non-reference reads are observed
                min_read_filter = args.filter_reads
                if nonref_reads_pool >= min_read_filter:
                    filtered = 'TRUE'
            else:
                # Filter if all the variants found in the probands are also found in the pool
                # This tends to result in multi-alleleic sites not getting filtered in many cases
                # Would need to split sites/individuals to fix this?
                if len(alleles_in_probands - alleles_in_pool) == 0:
                    filtered = 'TRUE'

            if filtered == 'TRUE':
                record.FILTER = 'InPool'

            falsepos = 'FALSE'
            # likely false positive if found in the pool but not in any of the probands
            if len(alleles_in_pool - alleles_in_probands) > 0:
                falsepos = 'TRUE'

            vcf_writer.write_record(record)

            outstream.write(','.join([str(x) for x in [var_id,nonref_alleles_pool,
                total_alleles_pool,nonref_alleles_probands,total_alleles_probands,
                nonref_reads_pool,total_reads_pool,filtered,falsepos,qual]]) + '\n')

if __name__ == '__main__':
    main()
