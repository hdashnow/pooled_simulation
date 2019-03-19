import os
import math
import sys

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

def count_nonref_reads(record_sample):
    """Count the number of reads supporting all non-reference alleles"""
    allelic_depths = record_sample['AD']
    try:
        nonref_reads = sum(allelic_depths[1:])
    except TypeError: # Occurs when AD is a single value, not a list
        nonref_reads = 0
    return(nonref_reads)

def alleles_supported(record, sample_pos, n, include_ref = True):
    """Return alleles supported by at least n reads
    Args:
        record (pyvcf record object)
        sample_pos (int): position of the desired sample in the vcf (0-based)
        n (int): return allele if it is supported by at least n reads
        include_ref (bool): if True, return the reference allele if it is supported by n reads
    Returns:
        list of numeric positions of alleles that are supported by at least n reads (strings)
        i.e. ref = '0', first nonref = '1'
    """
    # Get numeric representation of the alleles i.e. ref = 0, first alt = 1
    all_alleles = [str(x) for x in range(len(record.alleles))]
    allelic_depths = record.samples[sample_pos]['AD']
    if not include_ref:
        all_alleles = all_alleles[1:]
        allelic_depths = allelic_depths[1:]
    supported_alleles = [all_alleles[i] for i in range(len(all_alleles)) if allelic_depths[i] >= n]
    return(supported_alleles)

def is_recovered(alleles_in_probands, alleles_in_pool):
    """True if all the variants found in the proband(s) are also found in the pool.
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
    if total_reads_pool is None:
        total_reads_pool = 0

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

def check_pool_name(pool_name, all_vcf_samples):
    if pool_name:
        if pool_name not in all_vcf_samples:
            raise KeyError('ERROR: Specified pool name "{}" was not found in the VCF. Samples in the VCF include: {}'.format(pool_name, ', '.join(all_vcf_samples)))
    else:
        sys.stderr.write('WARNING: Pool name not given, inferred to be last sample in VCF.\n')
        pool_name = all_vcf_samples[-1]
    return pool_name

def check_proband_names(proband_names, all_vcf_samples):
    if proband_names:
        for proband in proband_names:
            if proband not in all_vcf_samples:
                raise KeyError('ERROR: Specified proband name "{}" was not found in the VCF. Samples in the VCF include: {}'.format(proband, ', '.join(all_vcf_samples)))
    else:
        sys.stderr.write('WARNING: Proband names not given, inferred to be all but the last sample in VCF.\n')
        proband_names = all_vcf_samples[:-1]
    return proband_names

def extract_record_info(record, field):
    try:
        info = record.INFO[field]
        if len(info) == 1:
            info = info[0]
        else:
            info = 'NA' # throw away gnomad data at multiallelic sites
        return info
    except KeyError:
        return 'NA'

def which_recovered(recovered_names, all_names):
    """Produce dictionary of recovered TRUE/FALSE"""
    recovery_dict = {}
    for name in all_names:
        if name in recovered_names:
            recovery_dict[name] = 'TRUE'
        else:
            recovery_dict[name] = 'FALSE'
    return(recovery_dict)
