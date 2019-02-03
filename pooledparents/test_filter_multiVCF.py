from filter_multiVCF import *
import pytest
import glob
import vcf

pytest_plugins = ["pytester"]

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

