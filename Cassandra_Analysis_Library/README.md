<>   ____    _    ____ ____    _    _   _ ____  ____      _    
<>  / ___|  / \  / ___/ ___|  / \  | \ | |  _ \|  _ \    / \   
<> | |     / _ \ \___ \___ \ / _ \ |  \| | | | | |_) |  / _ \  
<> | |___ / ___ \ ___) |__) / ___ \| |\  | |_| |  _ <  / ___ \ 
<>  \____/_/   \_\____/____/_/   \_\_| \_|____/|_| \_\/_/   \_\
<>    / \   _ __   __ _| |_   _ ___(_)___                      
<>   / _ \ | '_ \ / _` | | | | / __| / __|                     
<>  / ___ \| | | | (_| | | |_| \__ \ \__ \                     
<> /_/   \_\_| |_|\__,_|_|\__, |___/_|___/                     
<> | |   (_) |__  _ __ __ |___/_ _   _                         
<> | |   | | '_ \| '__/ _` | '__| | | |                        
<> | |___| | |_) | | | (_| | |  | |_| |                        
<> |_____|_|_.__/|_|  \__,_|_|   \__, |                        
                               |___/             

Cassandra Analysis Library
version 0.1

Last Updated (BY) 1/6/2018: restructured library.

Currently working on adding to the library and improving the capabilities of the scripts.

REQUIREMENTS: f2py, matplotlib, numpy

Addtional tools in progress:
	cas_cluster - Stillinger custer analysis
	cas_combine_boxes - creates single xyz file from GEMC xyz files. 
		(Creating this since there are some issues with topotools.)

To install:
	>> ./install.sh



To execute:

	-cas_density.py
		>> python cas_density.py -f 'xyzfilename'.xyz -m 'species1filename'.mcf 'species2filename'.mcf ... 

	-cas_rdf.py
		>> python cas_rdf.py -f 'xyzfilename'.xyz -m 'species1filename'.mcf 'species2filename'mcf ...

	-cas_angle.py
		>> python cas_angle.py -f 'xyzfilename'.xyz -m 'species1filename'.mcf 'species2filename'mcf ...
	
	cas_dihedral.py
		>> python cas_dihedral.py -f 'xyzfilename'.xyz -m 'species1filename'.mcf 'species2filename'mcf ...

	-cas_plot.py
		>> python cas_plot.py -f 'prpfilename'.prp


(To execute without having to type python, add the path to python in the shebang - first line
of each script.)


NOTE: Make sure the mcf files are in the correct order. Otherwise these scripts will not
work properly. As of now, there is no check for correctness in the order of the mcf files.

For more detail on these scripts and on additional optional flags execute the script with the
help flag (-h) in the command line
