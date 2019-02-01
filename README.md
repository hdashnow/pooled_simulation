[![Build Status](https://travis-ci.org/hdashnow/pooled_simulation.svg?branch=master)](https://travis-ci.org/hdashnow/pooled_simulation)

[![codecov](https://codecov.io/github/hdashnow/pooled_simulation/branch/master/graphs/badge.svg)](https://codecov.io/github/hdashnow/pooled_simulation)

Pooled parent code in development. Not yet ready for external use.

Usage:

Individual calling:

`bpipe run -p seed=1 pooled_sim_bpipe.groovy samples/*.bam`

`bpipe run -p seed=1 pooled_depth_sim_bpipe.groovy samples/*.bam`

Joint calling:

`bpipe run -p seed=1 pooled_sim_joint.groovy samples/*.bam samples/*.gvcf`

`bpipe run -p seed=1 pooled_depth_joint.groovy samples/*.bam samples/*.gvcf`

Install:

`git clone git@github.com:hdashnow/pooled_simulation.git`
`cd pooled_simulation`
`- pip install -U .`

Testing:

`python -m pytest pooledparents/test_*.py`
