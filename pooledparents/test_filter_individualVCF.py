from filter_individualVCF import *
import pytest
import os
import shutil
import glob

pytest_plugins = ["pytester"]

def test_parse_pool_specs(tmpdir):
    filename1 = tmpdir + '/2.txt'
    contents1 = '[/my/dir/sample1.bam, /my/dir/sample2.bam]'
    filename2 = tmpdir + '/4.txt'
    contents2 = '[/my/dir/sample1.bam, /my/dir/sample2.bam, /my/dir/sample3.bam, /my/dir/sample4.bam]'
    spec_files = [filename1, filename2]

    with open(filename1, 'w') as f1, open(filename2, 'w') as f2:
        f1.write(contents1)
        f2.write(contents2)

    expected = "{'2': ['sample1', 'sample2'], '4': ['sample1', 'sample2', 'sample3', 'sample4']}"
    assert str(parse_pool_specs(spec_files)) == expected

@pytest.fixture
def my_run(testdir, scriptdir = '.'):
    def do_run(*args, scriptdir):
        args = ['python', scriptdir+"filter_individualVCF.py"] + list(args)
        return testdir.run(*args)
    return do_run

@pytest.mark.parametrize("additional_args, exit_status", [
    ([], 0),
    (['--falsepos'], 0),
])
def test_end2end(tmpdir, my_run, additional_args, exit_status):
    dir_this_file = os.path.dirname(os.path.realpath(__file__)) + '/'
    datadir = dir_this_file + "test_data/"
    output = "pooled_sim_compare.csv"
    result = my_run("--individual_vcfs", *glob.glob(datadir+"SRR???????.vcf"),
                    "--pool_vcfs", datadir+"merge.2.RGfixed.dedup.vcf",
                    "--pool_specs", datadir+"2.txt",
                    "--out_csv", output,
                    *additional_args,
                    scriptdir = dir_this_file)
    assert result.ret == exit_status
    with open(output, "r") as f:
        newcontent = f.read()
