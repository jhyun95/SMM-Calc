# -*- coding: utf-8 -*-
"""
To display atoms in a separate process from the main code line

I NEED TO SPEED UP THE MAYAVI CODE USING 

obj.scene.disable_render = True
# Do all your scripting that takes ages.
# ...
# Once done, do the following:
obj.scene.disable_render = False


Daniel M Pajerowski
Daniel@Pajerowski.com
01.08.2012
"""
print 'please be patient while Mayavi renders...'

import CifFile  #existing library for reading in CIF format
import os       #operating system
import re       #regular expressions
import numpy as np  #numeric python
from enthought.mayavi import mlab
import SMMcalc_lib

f = open(r"SMMcalcINPUT.txt", "r").read().splitlines()

file_name = f[1]

#convert file_name to acceptable format
file_name = "file://" + file_name
file_name = os.path.normpath(file_name)

#read contents of TEMP file
[atoms, fract_coords] = SMMcalc_lib.tmp_read('TEMP.txt')

#read contents of CIF file
cf = CifFile.ReadCif(file_name)
sample_name = cf.dictionary.keys()[0]

#unit cell in meters (converted from Angstroms)
a =  float(re.sub('\)','',re.sub('\(','',cf[sample_name]['_cell_length_a']))) * 1e-10
b =  float(re.sub('\)','',re.sub('\(','',cf[sample_name]['_cell_length_b']))) * 1e-10
c =  float(re.sub('\)','',re.sub('\(','',cf[sample_name]['_cell_length_c']))) * 1e-10

#angles in degrees
alpha = float(re.sub('\)','',re.sub('\(','',cf[sample_name]['_cell_angle_alpha'])))
beta  = float(re.sub('\)','',re.sub('\(','',cf[sample_name]['_cell_angle_beta'])))
gamma = float(re.sub('\)','',re.sub('\(','',cf[sample_name]['_cell_angle_gamma'])))

#angles in radians
alpha_ = alpha*np.pi/180
beta_ = beta*np.pi/180
gamma_ = gamma*np.pi/180

#need to generate real space vectors...
a1 = np.array([a, 0, 0])
a2 = np.array([b*np.cos(gamma_), b*np.sin(gamma_)*np.sqrt(1-(np.cos(gamma_)*np.cos(beta_)-np.cos(alpha_))**2/(np.sin(gamma_)**2*np.sin(beta_)**2)), b*np.sin(gamma_)*(np.cos(gamma_)*np.cos(beta_)-np.cos(alpha_))/(np.sin(gamma_)*np.sin(beta_))])
a3 = np.array([c*np.cos(beta_), 0, c*np.sin(beta_)])

a1_A = a1*1e10
a2_A = a2*1e10
a3_A = a3*1e10

#start the mayavi pipeline
#    mlab.figure(fgcolor = (0,0,0), bgcolor = (1,1,1))

#define the color dictionary for atoms
color_dict = {'Gd':(0,0,1), 'Cu':(1,0,0)}

#display the bounding box of the unit cell
mlab.plot3d([0,a1_A[0]], [0,a1_A[1]], [0,a1_A[2]], line_width = 4)
mlab.plot3d([0,a2_A[0]], [0,a2_A[1]], [0,a2_A[2]], line_width = 4)
mlab.plot3d([a1_A[0], a1_A[0]+a2_A[0]], [a1_A[1], a1_A[1]+a2_A[1]], [a1_A[2], a1_A[2]+a2_A[2]], line_width = 4)
mlab.plot3d([a2_A[0], a1_A[0]+a2_A[0]], [a2_A[1], a1_A[1]+a2_A[1]], [a2_A[2], a1_A[2]+a2_A[2]], line_width = 4)

mlab.plot3d([a3_A[0]+0,a1_A[0]+a3_A[0]], [a3_A[1]+0,a1_A[1]+a3_A[1]], [a3_A[2]+0,a1_A[2]+a3_A[2]], line_width = 4)
mlab.plot3d([a3_A[0]+0,a2_A[0]+a3_A[0]], [a3_A[1]+0,a2_A[1]+a3_A[1]], [a3_A[2]+0,a2_A[2]+a3_A[2]], line_width = 4)
mlab.plot3d([a3_A[0]+a1_A[0], a1_A[0]+a2_A[0]+a3_A[0]], [a3_A[1]+a1_A[1], a1_A[1]+a2_A[1]+a3_A[1]], [a3_A[2]+a1_A[2], a1_A[2]+a2_A[2]+a3_A[2]], line_width = 4)
mlab.plot3d([a3_A[0]+a2_A[0], a1_A[0]+a2_A[0]+a3_A[0]], [a3_A[1]+a2_A[1], a1_A[1]+a2_A[1]+a3_A[1]], [a3_A[2]+a2_A[2], a1_A[2]+a2_A[2]+a3_A[2]], line_width = 4)

mlab.plot3d([0,a3_A[0]], [0,a3_A[1]], [0,a3_A[2]], line_width = 4)
mlab.plot3d([a1_A[0],a1_A[0]+a3_A[0]], [a1_A[1],a1_A[1]+a3_A[1]], [a1_A[2],a1_A[2]+a3_A[2]], line_width = 4)
mlab.plot3d([a2_A[0],a2_A[0]+a3_A[0]], [a2_A[1],a2_A[1]+a3_A[1]], [a2_A[2],a2_A[2]+a3_A[2]], line_width = 4)
mlab.plot3d([a2_A[0]+a1_A[0],a2_A[0]+a1_A[0]+a3_A[0]], [a2_A[1]+a1_A[1],a2_A[1]+a1_A[1]+a3_A[1]], [a2_A[2]+a1_A[2],a2_A[2]+a1_A[2]+a3_A[2]], line_width = 4)

r_coords = np.zeros(np.shape(fract_coords))
#take the fractional coordinates and generate real coordinates using the basis vectors
for i in range(np.size(atoms)):
    r_coords[i] = a1_A * fract_coords[i][0]  +  a2_A * fract_coords[i][1]  +  a3_A * fract_coords[i][2]
    
    
for i in range(np.size(atoms)):
    mlab.points3d(r_coords[i][0], r_coords[i][1], r_coords[i][2], opacity = 0.15, scale_factor = 1, color = color_dict[atoms[i]], mode = 'cube')
    mlab.text3d(r_coords[i][0], r_coords[i][1], r_coords[i][2], opacity = 0.35, scale = 0.5, text = atoms[i]+', '+str(i), color = color_dict[atoms[i]])

#    placeholder, delete when ready
chosen_atom_indices = []

mlab.show()