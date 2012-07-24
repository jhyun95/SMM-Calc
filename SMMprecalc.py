# -*- coding: utf-8 -*-
"""
For calculation of magnetic properties of paramagnetic single-molecule-magnet-like
materials, explicity paying attention to both intra- and inter-molecule
dipolar interactions, as well as superexchange within the cluster, single-ion-
anisotropy within the cluster, and cluster anisotropy

Daniel M Pajerowski
Daniel@Pajerowski.com
01.06.2012
"""
import SMMcalc_lib
import subprocess
import numpy as np


display_plots = 1
if display_plots == 0:
    print 'not displaying plots'
elif display_plots ==1:
    print 'displaying plot'

#
#read in the user defined parameters from the SMMcalcINPUT.txt file
#
f = open(r"SMMcalcINPUT.txt", "r").read().splitlines()

project_name = f.pop(0)
CIF_file = f.pop(0)
magnetization_file = f.pop(0)
magnetic_ions = f.pop(0).split(' ')
spin_list = [float(s) for s in f.pop(0).split(' ')]

spin_dict = {}
for i in range(np.size(spin_list)):
    spin_dict[magnetic_ions[i]] = spin_list[i]

print "project name: " + project_name

#
#read the CIF file to extract a complete unit cell containing the chosen magnetic ions
#
[atoms, fract_coords] = SMMcalc_lib.CIF_to_cell(file_name = CIF_file, select_atoms = magnetic_ions)

#
#store the atoms and fract_coords files in a temporary file so they can be acessed by other processes
#
SMMcalc_lib.tmp_write(atoms, fract_coords, 'TEMP.txt')

#
#start another process to render the unit cell in real space
#
if display_plots == 1:
    subprocess.Popen("SMMcalc_display.py", shell=True)

#
#choose the atoms to be in the cluster
#
cluster_indices = SMMcalc_lib.cluster_chooser()

#
#make sub-variables that only contain cluster atoms and fract_coords
#
[cluster_atoms, cluster_fract_coords] = \
SMMcalc_lib.cluster_array_maker(cluster_indices = cluster_indices, atoms = atoms, fract_coords = fract_coords)

#
#store the atoms and fract_coords files in a temporary file so they can be acessed by other processes
#
SMMcalc_lib.tmp_write(cluster_atoms, cluster_fract_coords, 'TEMP.txt')


#
#from the cluster_atoms, generate cluster_spins
#
cluster_spins = np.zeros(np.shape(cluster_atoms))
for i in range(np.size(cluster_atoms)):
    cluster_spins[i] = spin_dict[cluster_atoms[i]]


#
#start another process to render the chosen cluster in real space, with new indices
#
if display_plots ==1:
    subprocess.Popen("SMMcalc_display.py", shell=True)

#
#calculate the intra-cluster dipolar factor matrix
#one for each cartesian direction
#
[intra_dipolar_factor, intra_dipolar_indexer, intra_dipolar_revindexer] = SMMcalc_lib.intra_dipole(atoms = cluster_atoms, fract_coords = cluster_fract_coords)


#
#generate the spin operators
#
Sop = SMMcalc_lib.Sop_generator(cluster_spins)

#
#calculate the magnetization
#
#print SMMcalc_lib.MAGCAL(intra_dipolar_factor, intra_dipolar_indexer, cluster_spins, Sop, 2, 10)

#
#calculate the inter-cluster dipolar factor matrix
#one for each cartesian direction
#

#
#calculate the magnetization
#