from filtervcf import *
import pytest

def test_variant_id():
    pass

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
