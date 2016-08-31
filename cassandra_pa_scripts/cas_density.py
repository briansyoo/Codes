#!/afs/crc.nd.edu/x86_64_linux/python/2.7.4/bin/python

import numpy as np
import argparse
from argparse import RawTextHelpFormatter
import matplotlib.pyplot as plt

# F2PY Shared Objects
import cas_palib
# Cassandra Classes
from molecule_info import *
import readfiles 


#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
parser = argparse.ArgumentParser("Density Profile Analysis Tool",formatter_class=RawTextHelpFormatter)
parser.add_argument('-f',action='store',help ='Trajectory File (xyz format)')
parser.add_argument('-m',nargs = '+',action = 'store', help='MCF file(s) for each species. Must be in correct order.')
parser.add_argument('-H',action='store',help='Optional: H Matrix File (specify only if file name without extension is different from xyz file name)',default = 'default')
parser.add_argument('-nsl',action='store',help='Optional: Number of slices (default = 100)',default = 100,type =int)
parser.add_argument('-dim',action='store',help='Optional: Axis to analyze (default is Z-axis)',default = 'Z')
parser.add_argument('-b',action='store',help='Optional: Beginning frame (default = 1)', default = 0)
parser.add_argument('-e',action='store',help='Optional: Ending frame (default=last frame)', default=-1)
parser.add_argument('-o',action='store',help='Optional: Output file name',default='density.xvg')
#DEFINE PARSER ARGS
args = parser.parse_args()
nslices = int(args.nsl)
nspecies = len(args.m)
xyzfile = args.f
if args.H == 'default':
	Hfile = xyzfile[:-4]+'.H'
else:
	Hfile = args.H
begin_frame = args.b
end_frame = args.e
outfile = args.o
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#INITIALIZE
mcffiles = []
molecule_info = []
natoms=[]
max_natoms = 0
calc_tot_flag = False
density = np.zeros(nslices,np.float)

#READ MCF FILES
#create an object of molecule_info class for each species containing:
#        self.natoms = int
#        self.atom_name = array
#        self.atom_type = array
#        self.sigma = array
#        self.eps = array
#        self.mass = array
#        self.charge = array
#        self.name = char
#        self.MW = float
#        self.com = float

for i in range(nspecies):
	mcffile = args.m[i]
	mcffiles.append(mcffile)
	molecule_info.append(readfiles.read_mcf_file(mcffile))




# READ H-MATRIX FILE
#Read in the Hmatrix file to obtain the box lengths and the number of molecules
#for each snapshot. 
#Note: Matrices are flattened
H = readfiles.read_H_file(Hfile,nspecies)
Lx_array = H[0]
Ly_array = H[1]
Lz_array = H[2]

#Take the binwidth from the the average length over all frames divided by nslices
binwidth = np.mean(Lz_array)/nslices

#Store number of frames
nframes = int(len(H[3])/nspecies)

#Construct matrix that contains number of each species for each frame
nmolecules_array = np.reshape(H[3],(nframes,nspecies))


# (OPTIONAL) READ INPUT FILE
#out = readfiles.read_inp_file(inpfile)
#molecule = readmcf.read_mcf_file(mcffile)


#OBTAIN THE INDEX OF SPECIES OUTPUT
print """
What would you like to output?
 """
print "(0): System"
for i in range(nspecies):
	print "("+str(i+1)+"): ", molecule_info[i].name 
print "\n"
output_ndx = int(raw_input())-1

#Check if we want to calculate total density profile or not
if output_ndx==-1:
	print "\nObtaining total system density profile\n"
	calc_tot_flag = True 
else:
	print '\n'+"Obtaining density profile for "+molecule_info[output_ndx].name
print "\n"



#Store number of atoms of each species
for i in range(nspecies):
	natoms.append(molecule_info[i].natoms)


#Find max number of atoms in a species
for i in range(nspecies):
	if max_natoms < molecule_info[i].natoms:
		max_natoms = molecule_info[i].natoms



#Now we will fill in the indices with the mass. If no molecule is present,
#the mass should equal 0.

#Initialize mass for each atom in the system
atype_mass = np.zeros((max_natoms,nspecies))

#Loop through species and fill masses of each atom type
for i in range(nspecies):
	for j in range(max_natoms):
		try:
			atype_mass[j][i] = float(molecule_info[i].mass[j]) #convert to kg
		except:
			pass

#If we want to output a species density profile, we will set all other
#species atomtype masses to 0.

if not calc_tot_flag:
	for i in range(nspecies):
		if output_ndx == i:
			continue
		else:
			for j in range(max_natoms):
				atype_mass[j][i] *= np.double(0.0)

#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


#CALCULATE DENSITY PROFILE FROM FORTRAN SUBROUTINE
density = cas_palib.cas_density(xyzfile,
								nslices,
								nspecies,
								nframes,
								Lx_array,Ly_array, Lz_array,
								atype_mass,max_natoms, 
								natoms, nmolecules_array,density)


#PRINT RESULTS
Lzw = np.arange(0,nslices,1)*np.average(Lz_array)/nslices -np.average(Lz_array)/2.0
outfmt = '%12.4f%12.4f\n'
soutfmt = '%12.4f'
outfilew = open(outfile,'w')
print "Lz [Angstrom]", "Density [kg/m^3]\n"
for i in range(nslices):
	outfilew.write(outfmt%(Lzw[i],density[i]))
	print soutfmt%Lzw[i],soutfmt%density[i]
outfilew.close()
print "Average Density [kg/m^3]:\n"
print np.average(density)

#PLOT RESULTS INTO FIGURE
plt.figure(1)
plt.plot(Lzw,density)
plt.savefig(outfile[:-4]+'.png')

