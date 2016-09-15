# Cassandra Classes
from molecule_info import *

def read_atom_info(mcf,natoms):
	i = 0
	molecule = molecule_info()
	for line in mcf:
	
		atom_name = line.split()[1] #atom name
		atom_type = line.split()[2] #atom type
		mass = line.split()[3] #charge
		charge = line.split()[4]
		eps = line.split()[6] #epsilon
		sigma = line.split()[7] #sigma
		
		molecule.add_prop(natoms,atom_name, atom_type, sigma, eps, mass, charge)
		molecule.calc_MW(mass)	
		i+=1
		if i == natoms:
			return molecule

def read_mcf_file(mcffile):
	with open(mcffile,'r') as mcf:
		for line in mcf:
			if "# Atom_Info" in line:
				line = mcf.next()
				natoms = int(line.split()[0])
				molecule = read_atom_info(mcf,natoms)
				mcffile = mcffile.strip('./')
	molecule.add_molecule_name(mcffile.replace('.mcf',''))
	return molecule

def read_H_file(Hfile,nspecies):
	Lx_array = []
	Ly_array = []
	Lz_array = []
	nmolecules_array =[]
	with open(Hfile,'r') as H:
		for line in H:
			line = H.next()
			Lx_array.append(float(line.split()[0]))
			line = H.next()
			Ly_array.append(float(line.split()[1]))
			line = H.next()
			Lz_array.append(float(line.split()[2]))
			line = H.next()
			line = H.next()
			for i in range(nspecies):
				line = H.next()
				nmolecules_array.append(line.split()[1])
		return (Lx_array,Ly_array,Lz_array,nmolecules_array)


def read_inp_file(inpfile):
	nmolecules_array = []
	with open(inpfile,'r') as inp:
		for line in inp:
			if "# Sim_Type" in line:
				line = inp.next()
				sim_type = line.split()[0]
	with open(inpfile,'r') as inp:			
		for line in inp:
			if "# Nbr_Species" in line:
				line = inp.next()
				nspecies = int(line.split()[0])
	with open(inpfile,'r') as inp:
		for line in inp:
			if "# Molecule_Files" in line:
				print nspecies
				for i in range(nspecies):
					line = inp.next()
					nmolecules_array.append(int(line.split()[1]))
	
	return (sim_type,nspecies,nmolecules_array)
				
