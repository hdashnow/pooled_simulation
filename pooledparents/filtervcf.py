#import sys
#import vcf
#from vcf.parser import _Filter
import os
#from collections import Counter
#import copy
#import math

__author__ = "Harriet Dashnow"
__credits__ = ["Harriet Dashnow"]
__license__ = "MIT"
__version__ = "0.1.0"
__email__ = "h.dashnow@gmail.com"

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
    """Count the number of non-reference alleles from a VCF genotype (GT) string
    Returns the number of non-ref alleles and the total number of alleles"""
    alleles = GT_string.split('/')
    alleles = [allele for allele in alleles if allele != '.'] # remove missing genotypes
    total_alleles = len(alleles)
    nonref_allele_count = total_alleles - sum([allele == "0" for allele in alleles])
    return(nonref_allele_count, total_alleles)

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

