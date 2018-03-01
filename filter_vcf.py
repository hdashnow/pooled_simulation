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

def sample_id_from_fname(fname):
    sample_id = os.path.basename(fname).split('.')[0]
    if sample_id == 'merge':
        sample_id = os.path.basename(fname).split('.')[1]
    return(sample_id)

def count_nonref_alleles(GT_string):
    alleles = GT_string.split('/')
    nonref_allele_count = len(alleles) - sum([allele == "0" for allele in alleles])
    return(nonref_allele_count)

def get_gnomad_annotations(record):
    try:
        AF_EXOMESgnomad = record.INFO['AF_EXOMESgnomad'][0]
    except KeyError:
        AF_EXOMESgnomad = 'NA'
    try:
        Hom_Indiv_Exomes = record.INFO['Hom_Indiv_Exomes'][0]
    except KeyError:
        Hom_Indiv_Exomes = 'NA'
    try:
        AF_GENOMESgnomad = record.INFO['AF_GENOMESgnomad'][0]
    except KeyError:
        AF_GENOMESgnomad = 'NA'
    try:
        Hom_genomes = record.INFO['Hom_genomes'][0]
    except KeyError:
        Hom_genomes = 'NA'

    #print(AF_EXOMESgnomad, Hom_Indiv_Exomes, AF_GENOMESgnomad, Hom_genomes)
    #sys.exit()
    return {'AF_exome': AF_EXOMESgnomad, 'Hom_exome': Hom_Indiv_Exomes, 'AF_genome': AF_GENOMESgnomad, 'Hom_genome': Hom_genomes}

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
    var_dict = {}
    #pool_minor_allele_counts = []
    
    with open(parent_vcf_file, 'r') as this_vcf:
        for record in vcf.Reader(this_vcf):
            var_id = variant_id(record)
            parent_vars.add(var_id)

            pool_minor_allele_count = count_nonref_alleles(record.samples[0]['GT'])
            #pool_minor_allele_counts.append(pool_minor_allele_count)
            record_dict = get_gnomad_annotations(record)
            record_dict['pool_minor_allele_count'] = pool_minor_allele_count
            record_dict['proband_minor_allele_count'] = 'NA'
            record_dict['in_pool'] = 'TRUE'
            record_dict['in_any_proband'] = 'FALSE'
            var_dict[var_id] = record_dict
            
#    print('Parent pool')
#    print(len(parent_vars), 'variants')
#    print()
#    print(Counter(pool_minor_allele_counts))
#    print()

    #proband_vars = set()
    for vcf_file in proband_vcfs_files:
        proband_id = sample_id_from_fname(vcf_file)
        #proband_vars[proband_id] = set()
        
        with open(vcf_file, 'r') as this_vcf:
            vcf_reader = vcf.Reader(this_vcf)
            vcf_writer = vcf.Writer(open(proband_id+".filtered.vcf", 'w'), vcf_reader)
            for record in vcf_reader:

                var_id = variant_id(record)
                #proband_vars[proband_id].add(var_id)

                if var_id not in parent_vars:
                    vcf_writer.write_record(record)
                if var_id not in var_dict:
                    var_dict[var_id] = get_gnomad_annotations(record)
                    var_dict[var_id]['in_pool'] = 'FALSE'
                    var_dict[var_id]['pool_minor_allele_count'] = 'NA'
                    var_dict[var_id]['proband_minor_allele_count'] = 0
                if var_dict[var_id]['proband_minor_allele_count'] == 'NA':
                    var_dict[var_id]['proband_minor_allele_count'] = 0
                proband_minor_allele_count = count_nonref_alleles(record.samples[0]['GT'])
                var_dict[var_id]['proband_minor_allele_count'] += proband_minor_allele_count
                var_dict[var_id]['in_any_proband'] = 'TRUE'
                #proband_vars.add(var_id)

    header = '\t'.join(['variant', 'in_pool', 'in_any_proband', 'pool_minor_allele_count', 
        'proband_minor_allele_count', 'AF_exome', 'Hom_exome', 'AF_genome', 'Hom_genome']) + '\n'
    outstream.write(header)
    for var_id in var_dict:
        var_line = '{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(var_id, 
            var_dict[var_id]['in_pool'], var_dict[var_id]['in_any_proband'], 
            var_dict[var_id]['pool_minor_allele_count'], 
            var_dict[var_id]['proband_minor_allele_count'],
            var_dict[var_id]['AF_exome'], 
            var_dict[var_id]['Hom_exome'], var_dict[var_id]['AF_genome'], 
            var_dict[var_id]['Hom_genome'])
        outstream.write(var_line)
    outstream.close()        

#    for key in proband_vars:
#        print('Proband:', key)
#        print(len(proband_vars[key]), 'variants before filtering')
#        print(len(proband_vars[key] - parent_vars), 'variants after filtering')
#        print()

if __name__ == '__main__':
    main()
