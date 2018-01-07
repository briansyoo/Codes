#!/usr/bin/env python
import sys
from setuptools import setup

with open('README.md') as f:
	readme = f.read()

if sys.argv[-1] == 'readme':
    print(readme)
    sys.exit()


setup(name='CassandraAnalysisLibrary',
	version='0.1',
	description='Cassandra Analysis Library (CAL) for post analysis of CASSANDRA output files.',
	url='https://github.com/briansyoo/Codes',
	author='Brian Yoo',
	author_email='briansyoo@gmail.com',
	license='MIT',
	scripts=['bin/cal_angle.py','bin/cal_rdf.py','bin/cal_angle.py','bin/cal_dihedral.py','bin/cal_plot.py'],
	classifiers=[
		'Programming Language : Python :: 2.6',
		'Programming Language : Python :: 2.7',
	],
	install_requires=[
	'numpy','matplotlib',
	  ],
	  zip_safe=False)


