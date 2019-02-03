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
