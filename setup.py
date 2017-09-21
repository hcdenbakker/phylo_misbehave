import os
import shutil
import sys
import glob
from setuptools import setup, find_packages

def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()

version = 'x.y.z'
if os.path.exists('VERSION'):
  version = open('VERSION').read().strip()

setup(
    name='phylomisbehave',
    version=version,
    description='PhyloMisbehave: infer SNP-dense regions and the influence of different filtering strategies',
	long_description=read('README.md'),
    packages = find_packages(),
	package_data={'phylomisbehave': ['example_data/*', 'databases/*']},
    author='Henk den Bakker',
    author_email='xxxxxxxxxxxxx',
    url='https://github.com/sanger-pathogens/saffrontree',
    scripts=glob.glob('scripts/*'),
    test_suite='nose.collector',
    tests_require=['nose >= 1.3'],
    install_requires=[
			'ete3', 
			'cython',
			'seaborn',
			'matplotlib',
			'biopython >= 1.68'
       ],
    license='MIT',
    classifiers=[
        'Development Status :: 4 - Beta',
		'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Programming Language :: Python :: 3 :: Only',
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)'
    ],
)
