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
        '--probands', type=str, required=True, nargs='+',
        help='One or more vcfs containing proband variant calls')
    parser.add_argument(
        '--parentpool', type=str, required=True,
        help='Single vcf containing variant calls for the parent pool')
    parser.add_argument(
        '--output', type=str, required=False,
        help='Output file name. Defaults to stdout.')
    return parser.parse_args()

def variant_id(record):
    ALTstr = '/'.join([str(x) for x in record.ALT]) # join ALT loci
    POSstr = '{0:09d}'.format(record.POS) # add leading 0s
    return '_'.join([str(x) for x in [record.CHROM, POSstr, record.REF, ALTstr]])

def main():
    # Parse command line arguments
    args = parse_args()
    proband_vcfs_files = args.probands
    parent_vcf_file = args.parentpool
    outfile = args.output

    if outfile:
        outstream = open(outfile, 'w')
    else:
        outstream = sys.stdout

    parent_vars = set()
    minor_allele_counts = []
    
    with open(parent_vcf_file, 'r') as this_vcf:
        for record in vcf.Reader(this_vcf):
            var_id = variant_id(record)
            parent_vars.add(var_id)

            GT = record.samples[0]['GT']
            alleles = GT.split('/')
            minor_allele_count = sum([allele == "1" for allele in alleles])
            minor_allele_counts.append(minor_allele_count)

    print('Parent pool')
    print(len(parent_vars), 'variants')
    print()
    print(Counter(minor_allele_counts))

    proband_vars = {}
    print()

    for vcf_file in proband_vcfs_files:
        proband_id = os.path.basename(vcf_file).split('.')[0]
        proband_vars[proband_id] = set()
        
        with open(vcf_file, 'r') as this_vcf:
            vcf_reader = vcf.Reader(this_vcf)
            vcf_writer = vcf.Writer(open(proband_id+".filtered.vcf", 'w'), vcf_reader)
            for record in vcf_reader:

                var_id = variant_id(record)
                proband_vars[proband_id].add(var_id)
                if var_id not in parent_vars:
                    vcf_writer.write_record(record)

    for key in proband_vars:
        print('Proband:', key)
        print(len(proband_vars[key]), 'variants before filtering')
        print(len(proband_vars[key] - parent_vars), 'variants after filtering')
        print()

if __name__ == '__main__':
    main()
