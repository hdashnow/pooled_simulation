from filter_multiVCF import *
import pytest
import glob

pytest_plugins = ["pytester"]

def test_count_nonref_reads():
    pass
    
def test_alleles_supported():
    pass
    
def test_is_recovered():
    pass

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

@pytest.fixture
def my_run(testdir, scriptdir = '.'):
    def do_run(*args, scriptdir):
        args = ['python', scriptdir+"filter_multiVCF.py"] + list(args)
        return testdir.run(*args)
    return do_run

def test_end2end(tmpdir, my_run):
    dir_this_file = os.path.dirname(os.path.realpath(__file__)) + '/'
    datadir = dir_this_file + "test_data/"
    output = "pooled_sim_compare.csv"
    result = my_run("--in_vcf", datadir+"merge.2.pools_probands.vcf",
                    "--pool", '2',
                    "--probands", 'SRR1301698', 'SRR1301861',
                    "--out_csv", output,
                    scriptdir = dir_this_file)
    assert result.ret == 0
    with open(output, "r") as f:
        newcontent = f.read()

