#!/usr/bin/env python

from distutils.core import setup

LONG_DESCRIPTION = \
'''Pooled parent project code'''

with open('requirements.txt') as f:
    requirements = f.read().splitlines()

setup(
    name='pooledparent',
    version='0.1.0.0',
    author='Harriet Dashnow',
    author_email='h.dashnow@gmail.com',
    packages=['pooledparent'],
    package_dir={'pooledparent': 'pooledparent'},
    entry_points={
        'console_scripts': ['pooledparent = pooledparent.pooledparent:main']
    },
    url='https://github.com/hdashnow/pooled_simulation',
    license='LICENSE',
    description=('Pooled parent project'),
    long_description=(LONG_DESCRIPTION),
    install_requires=requirements,
)
