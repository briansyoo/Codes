#molecule_info.py


# This library assumes an MCF file has been read and parsed.
# We will store the molecule info into this class.
# Info to be stored:
# 
# -> atom_name
# -> atom_type 
# -> mass
# -> charge
# -> sigma
# -> eps
#			

import numpy as np

class molecule_info:
	def __init__(self):
		self.natoms = 0
		self.atom_name = []
		self.atom_type = []
		self.sigma = []
		self.eps = []
		self.mass = []
		self.charge = []
		self.name = '' 
		self.MW = 0
		self.com = 0

	def add_prop(self,number_atoms, atom_name,atom_type,sigma,eps,mass,charge):
		self.natoms =number_atoms
		self.atom_name.append(atom_name)
		self.atom_type.append(atom_type)
		self.sigma.append(sigma)
		self.eps.append(eps)
		self.mass.append(mass)
		self.charge.append(charge)

	def calc_com(self):
		return
	def calc_MW(self,mass):
		self.MW = self.MW +float(mass)
		#MW = np.sum(self.mass)
		#self.append(MW)
	def add_molecule_name(self,name):
		self.name = name
