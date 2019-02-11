from filtervcf import *
import pytest
import vcf

@pytest.mark.parametrize("position, expected", [
    (69270, 'chr1_000069270_A_G'),
    (894573, 'chr1_000894573_G_A'),
    (26083991, 'chr1_026083991_CGGGG_C/CG'),
])
def test_variant_id(position, expected):
    dir_this_file = os.path.dirname(os.path.realpath(__file__)) + '/'
    datadir = dir_this_file + "test_data/"
    vcf_file = datadir+"SRR1301861.vcf"
    with open(vcf_file, 'r') as this_vcf:
        vcf_reader = vcf.Reader(this_vcf)
        for record in vcf_reader:
            if record.POS == position:
                var_id = variant_id(record)
                assert var_id == expected

@pytest.mark.parametrize("fname, expected", [
    ('test.bam', 'test'),
    ('merge.10.bam', '10'),
])
def test_sample_id_from_fname(fname, expected):
    assert sample_id_from_fname(fname) == expected

@pytest.mark.parametrize("GT_string, expected", [
    ('0/0', (0,2)), # returns (nonref_allele_count, total_alleles)
    ('1/1', (2,2)),
    ('./.', (0,0)),
    ('./1', (1,1)),
    ('0/1/2/2', (3,4)),
])
def test_count_nonref_alleles(GT_string, expected):
    assert count_nonref_alleles(GT_string) == expected

@pytest.mark.parametrize("GT_string, expected", [
    ('0/0', set() ), # returns (nonref_allele_count, total_alleles)
    ('1/1', set('1') ),
    ('./.', set() ),
    ('./1', set('1') ),
    ('0/1/2/2', set(['1','2']) ),
])
def test_get_nonref_alleles(GT_string, expected):
    assert get_nonref_alleles(GT_string) == expected

@pytest.mark.parametrize("position, expected", [
    (69270, 7),
    (26083991, 89 + 27),
    (1560964, 0),
])
def test_count_nonref_reads(position, expected):
    dir_this_file = os.path.dirname(os.path.realpath(__file__)) + '/'
    datadir = dir_this_file + "test_data/"
    vcf_file = datadir + "merge.2.pools_probands.vcf"
    with open(vcf_file, 'r') as this_vcf:
        vcf_reader = vcf.Reader(this_vcf)
        for record in vcf_reader:
            if record.POS == position:
                nonref_reads = count_nonref_reads(record.samples[0])
                assert nonref_reads == expected

@pytest.mark.parametrize("position, min_read_filter, include_ref, expected", [
    (69270, 1, True, ['1']),
    (69270, 1, False, ['1']),
    (69270, 10, True, []),
    (26083991, 1, True, ['0', '1', '2']),
    (26083991, 1, False, ['1', '2']),
])
def test_alleles_supported(position, min_read_filter, include_ref, expected):
    dir_this_file = os.path.dirname(os.path.realpath(__file__)) + '/'
    datadir = dir_this_file + "test_data/"
    vcf_file = datadir + "merge.2.pools_probands.vcf"
    with open(vcf_file, 'r') as this_vcf:
        vcf_reader = vcf.Reader(this_vcf)
        for record in vcf_reader:
            if record.POS == position:
                alleles_by_reads = alleles_supported(record, 0,
                                    min_read_filter, include_ref)
                assert alleles_by_reads == expected

@pytest.mark.parametrize("alleles_in_probands, alleles_in_pool, expected", [
    ([1], [1], True),
    ([1, 2], [1], False),
])
def test_is_recovered(alleles_in_probands, alleles_in_pool, expected):
    assert is_recovered(alleles_in_probands, alleles_in_pool) == expected

@pytest.mark.parametrize("total_reads_pool,filter_reads,expected", [
    (100, 10, 10),
    (1, 1, 1),
    (1, 5, 5), #should probably throw error
    #(100, 0, 1), #should throw error
])
def test_set_read_filter_reads(total_reads_pool, filter_reads, expected):
    assert set_read_filter(total_reads_pool, filter_reads = filter_reads) == expected

@pytest.mark.parametrize("total_reads_pool, ploidy, tech_variation, expected", [
    (100, 10, 1, 10),
    (100, 10, 0.5, 5),
    (10, 1000, 1, 1),
    (None, 1000, 1, 1),
])
def test_set_read_filter_ploidy(total_reads_pool, ploidy, tech_variation, expected):
    assert set_read_filter(total_reads_pool, filter_reads = None, ploidy = ploidy,
        tech_variation = tech_variation) == expected

@pytest.mark.parametrize("total_reads_pool, filter_reads, ploidy, tech_variation, combine_filters, expected", [
    (100, 1, 10, 1, "strict", 10),
    (100, 1, 10, 1, "lenient", 1),
])
def test_set_read_filter_both(total_reads_pool, filter_reads, ploidy, tech_variation,
                            combine_filters, expected):
    assert set_read_filter(total_reads_pool, filter_reads = filter_reads, ploidy = ploidy,
        tech_variation = tech_variation, combine_filters = combine_filters) == expected

def test_set_read_filter_wrongargs():
    with pytest.raises(ValueError):
        set_read_filter(10, filter_reads = 1, ploidy = 2, tech_variation = 0.5, combine_filters = 'max')

@pytest.mark.parametrize("pool_name, all_vcf_samples, expected", [
    (None, ['proband1', 'proband2', 'pool'], 'pool'),
    ('pool', ['proband1', 'pool', 'proband2'], 'pool'),
])
def test_check_pool_name(pool_name, all_vcf_samples, expected):
    assert check_pool_name(pool_name, all_vcf_samples) == expected

def test_check_pool_name_error():
    with pytest.raises(KeyError):
        assert check_pool_name('other', ['proband1', 'pool', 'proband2'])

@pytest.mark.parametrize("proband_names, all_vcf_samples, expected", [
    (None, ['proband1', 'proband2', 'pool'], ['proband1', 'proband2']),
    (['proband1', 'proband2'], ['proband1', 'pool', 'proband2'], ['proband1', 'proband2']),
])
def test_check_proband_names(proband_names, all_vcf_samples, expected):
    assert check_proband_names(proband_names, all_vcf_samples) == expected

def test_check_proband_names_error():
    with pytest.raises(KeyError):
        assert check_proband_names(['some', 'other'], ['proband1', 'pool', 'proband2'])

@pytest.mark.parametrize("position, field, expected", [
    (69270, 'AF_EXOMESgnomad', 'NA'),
    (26083991, 'AF_GENOMESgnomad', 'NA'),
    (1560964, 'other', 'NA'),
])
def test_extract_record_info(position, field, expected):
    dir_this_file = os.path.dirname(os.path.realpath(__file__)) + '/'
    datadir = dir_this_file + "test_data/"
    vcf_file = datadir + "merge.2.pools_probands.vcf" #XXX this test file isn't annotated, replace
    with open(vcf_file, 'r') as this_vcf:
        vcf_reader = vcf.Reader(this_vcf)
        for record in vcf_reader:
            if record.POS == position:
                assert extract_record_info(record, field) == expected

@pytest.mark.parametrize("recovered_names, all_names, expected", [
    (['A'], ['A', 'B'], {'A': 'TRUE', 'B': 'FALSE'}),
    ([], ['A', 'B'], {'A': 'FALSE', 'B': 'FALSE'}),
])
def test_which_recovered(recovered_names, all_names, expected):
    assert which_recovered(recovered_names, all_names) == expected
