from filter_multiVCF import *
import pytest

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
