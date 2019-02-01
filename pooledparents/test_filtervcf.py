from filtervcf import *
import pytest

@pytest.mark.parametrize("fname, expected", [
    ('test.bam', 'test'),
])
def test_sample_id_from_fname(fname, expected):
    assert sample_id_from_fname(fname) == expected

