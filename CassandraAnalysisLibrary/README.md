<!---
==============================================================
   ____    _    ____ ____    _    _   _ ____  ____      _    
  / ___|  / \  / ___/ ___|  / \  | \ | |  _ \|  _ \    / \   
 | |     / _ \ \___ \___ \ / _ \ |  \| | | | | |_) |  / _ \  
 | |___ / ___ \ ___) |__) / ___ \| |\  | |_| |  _ <  / ___ \ 
  \____/_/   \_\____/____/_/   \_\_| \_|____/|_| \_\/_/   \_\
    / \   _ __   __ _| |_   _ ___(_)___                      
   / _ \ | '_ \ / _` | | | | / __| / __|                     
  / ___ \| | | | (_| | | |_| \__ \ \__ \                     
 /_/   \_\_| |_|\__,_|_|\__, |___/_|___/                     
 | |   (_) |__  _ __ __ |___/_ _   _                         
 | |   | | '_ \| '__/ _` | '__| | | |                        
 | |___| | |_) | | | (_| | |  | |_| |                        
 |_____|_|_.__/|_|  \__,_|_|   \__, |                        
                               |___/             
==============================================================
-->
# Welcome to Cassandra Analysis Library version 0.1!


Last Updated 1/6/2018: Restructured library to include Python egg installation. Code is currently broken.

Addtional tools in progress:

	cal_cluster - Stillinger custer analysis

Currently working on adding to the library and improving the capabilities of the scripts.
Future implementation will implement machine learning algorithms from scikitlearn.

REQUIREMENTS: gcc, gfortran, python2.7 (recommended), f2py, matplotlib, numpy


## To install with root access:

	1. Compile fortran subroutines

	>> ./install.sh

	2. Installation with setup.py

	>> python setup.py install

## To install without root access:

	1. Compile fortran subroutines

	>> ./install.sh

	2. Add following lines to .bashrc file (setup of local env)

	PYTHONPATH="${PYTHONPATH}:path/to/usr/.local/lib/python/"
	export PATH=path/to/usr/.local/bin:$PATH

	3. Source .bashrc file

	>> source .bashrc

	4. Exceute installation with setup.py

	>> python setup.py install --home=~/.local


## Example Usage

NOTE: Make sure the mcf files are in the correct order. Otherwise these scripts will not
work properly. As of now, there is no check for correctness in the order of the mcf files.

For more detail on these scripts and on additional optional flags execute the script with the
help flag (-h) in the command line


### To Execute:

	Density Profile Analysis

	>>cal_density.py -f 'xyzfilename'.xyz -m 'species1filename'.mcf 'species2filename'.mcf ... 

	Radial Distribution Analysis

	>>cal_rdf.py -f 'xyzfilename'.xyz -m 'species1filename'.mcf 'species2filename'mcf ...

	Angle Distribution Analysis

	>>cal_angle.py -f 'xyzfilename'.xyz -m 'species1filename'.mcf 'species2filename'mcf ...
	
	Dihedral Angle Analysis

	>>cal_dihedral.py -f 'xyzfilename'.xyz -m 'species1filename'.mcf 'species2filename'mcf ...

	Time Series Plot Analysis

	>>cal_plot.py -f 'prpfilename'.prp

