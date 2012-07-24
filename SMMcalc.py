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
#import subprocess
import os
import matplotlib
import matplotlib.pyplot
import numpy as np
import scipy.optimize as spop
#import time

def resids(x):
    """minimize the residuals of a fit to the expiermental
    data
    args is the error    
    """
#    t0 = time.time()
    print 'started resids'
    residuals = 0
    for i in range(np.size(data_T)):
        residuals = residuals + 1*(1/data_errt[i])*(SMMcalc_lib.MAGCAL_powder(x[0],x[1],x[2],\
        superexchange_matrix, single_ion_matrix, intra_dipolar_matrix, cluster_spins, Sop, data_T[i], .1)*data_T[i]\
        -data_Mt[i]*data_T[i])**2
    for i in range(np.size(data_H)):
        residuals = residuals + 00*(1/data_errh[i])*(SMMcalc_lib.MAGCAL_powder(x[0],x[1],x[2],\
        superexchange_matrix, single_ion_matrix, intra_dipolar_matrix, cluster_spins, Sop, 2, data_H[i])\
        -data_Mh[i])**2
    print x, residuals
    return residuals

#############def resids_(g, J, D):
#############    """v are the parameters to be varied to make fp(v,x)=y
#############    x is the domain of the function, passed as an argument
#############    y is the experimental data to be fit
#############    fp(v,x) is the model
#############    """
#############    return (fp(g, J, D)-data_I)
#############
#############def fp(params ,two_theta):
#############    """
#############    the model to fit the experimental data
#############    v is the parameter set
#############    x is the domain of the function
#############    """
#############    y0_all = params[0]
#############    u_,v_,w_ = params[1:4]
#############    JzCo = params[4]
#############    JzFe = -abs(params[5]) #ferromagnetism
#############    FM_ = rietveld_lib.FM(atoms, fract_coords, Q_fullprof, JzCo, JzFe, is_orb_Co, is_orb_Fe)
#############    I = np.multiply(rietveld_lib.FNtoI(FM_, Mult_fullprof, TwoTheta_fullprof),  0.00106921*1.1)
#############    y_fit = np.add(np.zeros(np.shape(two_theta)), y0_all)
#############    for j in range(len(Q_fullprof)):
#############        xc = rietveld_lib.bragg_angle( (10.227, 10.227, 10.227), (90,90,90), Q_fullprof[j] )
#############        A = I[j]
#############        w = (  u_*np.tan(xc*np.pi/360)**2  +  v_*np.tan(xc*np.pi/360)  +  w_  )**(0.5)  *  5#   *  width_scalar
#############        y0 = 0
#############        
#############
#############        for i in range(np.size(two_theta)):
#############            y_fit[i] = y_fit[i] + (  rietveld_lib.gaussian(two_theta[i], (y0, A, xc, w))  )
#############            
#############
#############    return y_fit

#
#read in the user defined parameters from the SMMcalcINPUT.txt file
#
f = open(r"SMMcalcINPUT.txt", "r").read().splitlines()

project_name = f.pop(0)
CIF_file = f.pop(0)
magnetization_file = os.path.normpath(f.pop(0)).split(',')
magnetic_ions = f.pop(0).split(' ')
spin_list = [float(s) for s in f.pop(0).split(' ')]
cluster_indices = [int(s) for s in f.pop(0).split(' ')]
D_vector_indices = [int(s) for s in f.pop(0).split(' ')]
superexchange_indices = f.pop(0).split(' ')

#J = []
#for i in superexchange_indices:
#    J[i][0] = 

spin_dict = {}
for i in range(np.size(spin_list)):
    spin_dict[magnetic_ions[i]] = spin_list[i]

print "project name: " + project_name

#
#read the CIF file to extract a complete unit cell containing the chosen magnetic ions
#
[atoms, fract_coords] = SMMcalc_lib.CIF_to_cell(file_name = CIF_file, select_atoms = magnetic_ions)


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


D_vector = SMMcalc_lib.D_vector_calc(D_vector_indices = D_vector_indices, fract_coords = fract_coords)

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
#generate the matrices for hamiltonian
#
intra_dipolar_matrix = SMMcalc_lib.intra_dipolar_matrix(intra_dipolar_factor, intra_dipolar_indexer, intra_dipolar_revindexer, Sop, cluster_spins)
single_ion_matrix = SMMcalc_lib.single_ion_matrix(D_vector, Sop, cluster_spins)
superexchange_matrix = SMMcalc_lib.superexchange_matrix(Sop, cluster_spins)

if magnetization_file == 'simulation':
    print 'just a simulation, folks.  nothing to see here.'
else:
    if np.size(magnetization_file) == 2:
        ft = open(magnetization_file[0], "r").read().splitlines()
        fh = open(magnetization_file[1], "r").read().splitlines()
    
        data_units = ft.pop(0)
        
        data_T = []
        data_Mt = []
        data_errt = []
        for i in range(np.size(ft)):
            data_T.append( float(ft[i].split()[0]) )
            data_Mt.append( float(ft[i].split()[1]) )
            data_errt.append( float(ft[i].split()[2]) )

        data_H = []
        data_Mh = []
        data_errh = []
        for i in range(np.size(fh)):
            data_H.append( float(fh[i].split()[0]) )
            data_Mh.append( float(fh[i].split()[1]) )
            data_errh.append( float(fh[i].split()[2]) )
        
        matplotlib.pyplot.figure(1)
        matplotlib.pyplot.plot(data_T, np.multiply(data_Mt, data_T), 'bo' )
        matplotlib.pyplot.errorbar(data_T, np.multiply(data_Mt, data_T), yerr = np.multiply(data_errt,data_T), ecolor='b' )
        
        matplotlib.pyplot.figure(2)
        matplotlib.pyplot.plot(data_H, data_Mh, 'bo' )
        matplotlib.pyplot.errorbar(data_H, data_Mh, yerr = data_errh, ecolor='b' )

        
        xm = spop.optimize.fmin(resids,x0=(2,8,0,))
    
        print xm
        #
    #    matplotlib.pyplot.plot(data_T, np.multiply(xm[0]/(data_T+xm[1])+xm[2]*0, data_T), 'rx' )
    
    #
    #calculate the magnetization
    #
    
        
        calc_Mt = []
        for i in range(np.size(data_T)):
            calc_Mt.append(SMMcalc_lib.MAGCAL_powder(xm[0],xm[1],xm[2], superexchange_matrix, single_ion_matrix, intra_dipolar_matrix, cluster_spins, Sop, data_T[i], .1))
        matplotlib.pyplot.figure(1)
        matplotlib.pyplot.plot(data_T, np.multiply(calc_Mt, data_T), 'k-' )    
    
        calc_Mh = []
        for i in range(np.size(data_H)):
            calc_Mh.append(SMMcalc_lib.MAGCAL_powder(xm[0],xm[1],xm[2], superexchange_matrix, single_ion_matrix, intra_dipolar_matrix, cluster_spins, Sop, 2, data_H[i]))
        
        matplotlib.pyplot.figure(2)
        matplotlib.pyplot.plot(data_H, calc_Mh, 'k-' )
    
        matplotlib.pyplot.show()
    
        print xm[0],xm[1],xm[2]
#
#calculate the inter-cluster dipolar factor matrix
#one for each cartesian direction
#

#
#calculate the magnetization
#