#!/usr/bin/env python

if sys.argv[-1] == 'readme':
    print(readme)
    sys.exit()

m setuptools import setup

setup(name='funniest',
      version='0.1',
      description='Cassandra Analysis Library (CAL) for post analysis of CASSANDRA output files.',
      url='http://github.com/storborg/funniest',
      author='Brian Yoo',
      author_email='briansyoo@gmail.com',
      license='MIT',
      packages=['funniest'],
      install_requires=[
          'numpy',
      ],
      zip_safe=False)


