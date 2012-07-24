# -*- coding: utf-8 -*-
"""
Module containing functions utilized in SMMcalc.

Daniel M Pajerowski
Daniel@Pajerowski.com
01.06.2012
"""

def CIF_to_cell(file_name = "C:/LanthanideSMMs/GdCu SMMs/7.cif", select_atoms = []):
    """read in the contents of a CIF file
    
    have the option of only including a specific sub-set of atoms in the cell
    
    inputs:
        file_name: a string with front slashes for file tree delimiters that directs to the location of the CIF file
        select_atoms: a list of atoms to use, if left blank, all atoms are used
    
    output:
        [atom_label, fractional_coordinates]
        
    """
    import CifFile  #existing library for reading in CIF format
    import os       #operating system
    import re       #regular expressions
    import numpy as np  #numeric python
    import SpaceGroups  #the DANSE spacegroups package
      
    #convert file_name to acceptable format
#    file_name = "file://" + file_name
    file_name = os.path.normpath("file://" + file_name)
    
    #read contents of CIF file
    cf = CifFile.ReadCif(file_name)
    sample_name = cf.dictionary.keys()[0]
    
    #make a complete select_atoms list if a null list was provided
    if select_atoms == []:
        select_atoms = list(set(cf[sample_name]['_atom_site_type_symbol']))
    
    #extract the types of the selected atoms and their fractional coordinates
    i = 0
    atoms = []
    fract_coords = []
    for j in cf[sample_name]['_atom_site_type_symbol']:
        for k in select_atoms:
            if j == k:
    
                x = float(re.sub('\)','',re.sub('\(','',cf[sample_name]['_atom_site_fract_x'][i])))
                y = float(re.sub('\)','',re.sub('\(','',cf[sample_name]['_atom_site_fract_y'][i])))
                z = float(re.sub('\)','',re.sub('\(','',cf[sample_name]['_atom_site_fract_z'][i])))
    
                fract_coords.append(np.array([x,y,z]))  #fraction coordinates
                atoms.append(j)                         #atom types
        i = i + 1
    
    #take the given atoms and populate a complete P1-type unit cell
    
    #first pick out the space group
    space_group_name = re.sub('\)','',re.sub('\(','',cf[sample_name]['_symmetry_space_group_name_H-M']))
    
    print 'be very careful that the correct setting is being used for the space group'
    
    for i in range(np.size(SpaceGroups.SpaceGroupList)):
         if space_group_name == SpaceGroups.SpaceGroupList[i].short_name:
             space_group = SpaceGroups.SpaceGroupList[i]
    
    #then step through each atom, and apply symmetry operators to populate the cell
    for i in range(np.size(atoms)):
        for j in space_group.symop_list:
    #        apply the symmetry operator to 'coords' to make a new vector
            fract_coords_new = np.dot(j.R, fract_coords[i])+j.t
    #        if the new vector is outside the primitive cell, move it back
            for k in [0,1,2]:
                while fract_coords_new[k] >= 1:
                    fract_coords_new[k] = fract_coords_new[k]-1
                while fract_coords_new[k] < 0:
                    fract_coords_new[k] = fract_coords_new[k]+1
    #        add the new position to the list if it isn't already there
    #        check to within a tolerance to avoid numeric errors
            sym_tol = 0.0001
            degtest = 0
            for l in range(np.size(atoms)):
                if (fract_coords_new[0] > fract_coords[l][0]-sym_tol and fract_coords_new[0] < fract_coords[l][0]+sym_tol)\
                and (fract_coords_new[1] > fract_coords[l][1]-sym_tol and fract_coords_new[1] < fract_coords[l][1]+sym_tol)\
                and (fract_coords_new[2] > fract_coords[l][2]-sym_tol and fract_coords_new[2] < fract_coords[l][2]+sym_tol):
                    degtest = 1
                    
            if degtest == 0:
                fract_coords.append(fract_coords_new)
                atoms.append(atoms[i])
    return [atoms, fract_coords]

def atom_chooser(file_name = "C:/LanthanideSMMs/GdCu SMMs/7.cif", atoms = [], fract_coords = []):
    """display a list of atoms using a list of fractional coordinate arrays
    
    unit cell parameters will be extracted from the CIF-file
    
    inputs:
        file_name: a string with front slashes for file tree delimiters that directs to the location of the CIF file
        atoms: string list of atoms index in same manner as fract_coords
        fract_coords: a list of numpy arrays containing the fractional coordinates
    
    output:
        chosen_atom_indices
        
    """
    
    import CifFile  #existing library for reading in CIF format
    import os       #operating system
    import re       #regular expressions
    import numpy as np  #numeric python
    from enthought.mayavi import mlab
      
    #convert file_name to acceptable format
    file_name = "file://" + file_name
    file_name = os.path.normpath(file_name)
    
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
        mlab.points3d(r_coords[i][0], r_coords[i][1], r_coords[i][2], opacity = 0.15, scale_factor = 1, color = (0,0,1), mode = 'cube')
        mlab.text3d(r_coords[i][0], r_coords[i][1], r_coords[i][2], opacity = 0.35, scale = 0.5, text = atoms[i], color = (0,0,1))

#    placeholder, delete when ready
    chosen_atom_indices = []
    
    mlab.show()    
    
    return chosen_atom_indices
    
def tmp_write(atoms = [], fract_coords = [], file_name =[]):
    """write the data to a temp file
    
    inputs: 
        atoms: a list of atoms
        fract_coords: a list of fractional coordinate arrays
        file_name: name of file to write to
    """
    import numpy as np
    
    f = open(file_name, 'w')
    
    for i in range(np.size(atoms)):
        print>>f, atoms[i],fract_coords[i][0],fract_coords[i][1],fract_coords[i][2]
        
    f.close()

def tmp_read(file_name = []):
    """read from a temp data file
    
    inputs:
        file_name: the name of the temp file to read
    
    outputs:
        atoms: a list of atoms
        fract_coords: a list of fractional coordinate arrays
    """
    import numpy as np
    
    f = open(file_name, "r").read().splitlines()

    atoms = [' ']*np.size(f)
    fract_coords = np.zeros((np.size(f),3))
    
    for i in range(np.shape(f)[0]):
        j = f[i].split(' ')
        atoms[i] =  j[0]
        fract_coords[i] = np.array([float(j[1]), float(j[2]), float(j[3])])
        
    return atoms, fract_coords
    
    
def cluster_chooser():
    """query the user for what atoms to be in the magnetic cluster
    
    inputs:
        none
        
    outputs:
        cluster_atoms: list of indices of chosen atoms
    """
    cluster_indices = []
    cluster_inp = []
    while cluster_inp != 'X' and cluster_inp != 'x':    
        cluster_inp = raw_input('define cluster atoms (type X when done):')
        if cluster_inp.isdigit() == True:
            cluster_inp = int(cluster_inp)
            cluster_indices.append(int(cluster_inp))
        elif cluster_inp == 'X' or cluster_inp == 'x':
            print 'OK!'
        else:
            print 'invalid choice, must be a number'
    return cluster_indices

def cluster_array_maker(cluster_indices = [], atoms = [], fract_coords = []):
    """make lists containing cluster information
    
    inputs:
        cluster_indices:
        atoms:
        fract_coords:

    outputs:
        cluster_atoms:
        cluster_fract_coords:
    """
    import numpy as np
    
    cluster_atoms = [' ']*np.size(cluster_indices)
    cluster_fract_coords = np.zeros((np.size(cluster_indices),3))
    j = 0
    for i in cluster_indices:
        cluster_atoms[j] = atoms[i]
        cluster_fract_coords[j] = fract_coords[i]
        j = j+1
        
    return [cluster_atoms, cluster_fract_coords]

def intra_dipole(atoms = [], fract_coords = []):
    """calculate the geometrical part of the dipolar field at atom 'i' due to atom 'j'
    here, the (mu0/4.0/np.pi) prefactor is taken out, as are the scalar aspects
    of the spins
    
    inputs:
        atoms: the list of atoms, only the length is actually utilized
        fract_coords: the fractional coordinates of the atoms in the chosen cluster
                        real-space coordinates will be calculated by reloading
                        CifFile
    
    outputs:
        factor: first index denotes location atom, second index denotes source atom
        indexer: to be used on factor to make it easier to access
                such that 'indexer' takes two arguments, and the first
                refers to the spin in question, and the second refers to the
                spatial part
        revindexer: the opposite of indexer, taking one argument that is
                in reference to the the factor index, and returns
                which spin and spatial part factor is referencing
    """
    import CifFile  #existing library for reading in CIF format
    import os       #operating system
    import re       #regular expressions
    import numpy as np  #numeric python
    
    f = open(r"SMMcalcINPUT.txt", "r").read().splitlines()
    
    file_name = f[1]
    
    #convert file_name to acceptable format
    file_name = "file://" + file_name
    file_name = os.path.normpath(file_name)
    
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
    
    r_coords = np.zeros(np.shape(fract_coords))
    #take the fractional coordinates and generate real coordinates using the basis vectors
    for i in range(np.size(atoms)):
        r_coords[i] = a1 * fract_coords[i][0]  +  a2 * fract_coords[i][1]  +  a3 * fract_coords[i][2]

    indexer = np.zeros((np.size(atoms),3), dtype = np.int8)

    k = 0
    for i in range(np.size(atoms)):   #first index is for which spin
        for j in range(3):          #second index is for spatial part
            indexer[i][j] = int(k)
            k = k+1
    
    revindexer = []
    
    k = 0
    for i in range(np.size(atoms)):   #first index is for which spin
        for j in range(3):          #second index is for spatial part
            revindexer.append(np.array([int(i),int(j)]))
    
    factor = np.zeros((np.size(atoms)*3, np.size(atoms)*3))
    
    unit_vector = []
    
    for i in range(3):
        unit_vector.append(np.zeros(3))
        unit_vector[i][i] = 1
    
    for i in range(np.size(atoms)*3):
        for j in range(np.size(atoms)*3):
            atom1 = revindexer[i][0]
            atom2 = revindexer[j][0]
            if atom1 != atom2:
                r = r_coords[atom1] - r_coords[atom2]
                r_mag = np.sqrt(np.dot(r,r))
                
                m1 = unit_vector[revindexer[i][1]]
                m2 = unit_vector[revindexer[j][1]]
    #            (mu0/4.0/np.pi)  *  
    
                factor[i][j] = np.dot(  (  3.0*r*np.dot(m1,r)/r_mag**5  -  m1/r_mag**3  )  ,  m2)
#                factor[i][j] = (mu0/4.0/np.pi)  * factor[i][j]
    
    return [factor, indexer, revindexer]

def D_vector_calc(D_vector_indices = [], fract_coords = []):
    """calculate the geometrical part of the dipolar field at atom 'i' due to atom 'j'
    here, the (mu0/4.0/np.pi) prefactor is taken out, as are the scalar aspects
    of the spins
    
    inputs:
        D_vector_indices:
        fract_coords: the fractional coordinates of the atoms in the chosen cluster
                        real-space coordinates will be calculated by reloading
                        CifFile
    
    outputs:
        D_vector:
    """
    import CifFile  #existing library for reading in CIF format
    import os       #operating system
    import re       #regular expressions
    import numpy as np  #numeric python
    
    f = open(r"SMMcalcINPUT.txt", "r").read().splitlines()
    
    file_name = f[1]
    
    #convert file_name to acceptable format
    file_name = "file://" + file_name
    file_name = os.path.normpath(file_name)
    
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
    
    #take the fractional coordinates and generate real coordinates using the basis vectors
    
    D_vector = np.add((a1 * fract_coords[D_vector_indices[0]][0]  +  a2 * fract_coords[D_vector_indices[0]][1]  +  a3 * fract_coords[D_vector_indices[0]][2])\
    , -(a1 * fract_coords[D_vector_indices[1]][0]  +  a2 * fract_coords[D_vector_indices[1]][1]  +  a3 * fract_coords[D_vector_indices[1]][2]))

    D_vector = np.multiply(D_vector, np.sqrt(np.dot(D_vector,D_vector))**(-1))

    return D_vector

def Sop_generator(cluster_spins = []):
    import numpy as np
    import copy
    import time
    
    def multi_for(iterables):
         if not iterables:
             yield ()
         else:
             for item in iterables[0]:
                 for rest_tuple in multi_for(iterables[1:]):
                     yield (item,) + rest_tuple
    
    def spin_states(max_spin):
        return np.linspace(-max_spin, max_spin, (2*max_spin+1))
        
    TEMPTIME = time.localtime()
    print 'starting Sop_generator:'+str(TEMPTIME.tm_hour)+':'+str(TEMPTIME.tm_min)+':'+str(TEMPTIME.tm_sec)
    t0 = time.time()   
    
    #in a scalable manner, with z basis, initialize the list to store states
    state = []
    for i in multi_for(map(spin_states, cluster_spins)):
        state.append(list(i))
    
    #initialize the spin-operators
    #i.e. Sp[1][2] gives the Sz operator for the 1st ion in the cluster
    #0 = x, 1 = y, 2 = z, 3 = plus, 4 = minus
    Sop = [  [np.zeros([np.size(state,0),np.size(state,0)]) for j in range(5)]  for k in range(np.size(cluster_spins))]
    
    print 'initialized things:'+str(time.time()-t0)
    
    #for each state, need to calculate matrix elements, hence the double-sum
    for i in range(0,np.size(state,0)):     #state i
        for j in range(0,np.size(state,0)): #state j
            for k in range(np.size(cluster_spins)): #for the kth atom
    
        #        build up the copper spin-operators, the z operator
                if state[i] == state[j]:
                    Sop[k][2][i][j] = state[i][k]  #state[i][k] is the kth Sz of the ith state
                else:
                    Sop[k][2][i][j] = 0
        #        operators acting on state[j][k], the + operator
                p_test = copy.deepcopy(state[j])
                p_test[k] = p_test[k]+1
                if state[i] == p_test:
                    Sop[k][3][i][j] = np.sqrt( (cluster_spins[k]+state[j][k]+1)*(cluster_spins[k]-state[j][k])  )
    #                print 'p_test on '+str(k)
    #                print 'jth state:'+str(state[j])+', ith state:'+str(state[i])+', p_test:'+str(p_test)
    #                print np.sqrt( (cluster_spins[k]+state[j][k]+1)*(cluster_spins[k]-state[j][k])  )
                else:
                    Sop[k][3][i][j] = 0
        #        the - operator
                m_test = copy.deepcopy(state[j])
                m_test[k] = m_test[k]-1
                if state[i] == m_test:
                    Sop[k][4][i][j] = np.sqrt( (cluster_spins[k]-state[j][k]+1)*(cluster_spins[k]+state[j][k])  )
    #                print 'm_test on '+str(k)
    #                print 'jth state:'+str(state[j])+', ith state:'+str(state[i])+', m_test:'+str(m_test)
    #                print np.sqrt( (cluster_spins[k]-state[j][k]+1)*(cluster_spins[k]+state[j][k])   )
                else:
                    Sop[k][4][i][j] = 0
    
    print 'made Sz, S+, and S-:'+str(time.time()-t0)
    
    #from the plus and minus matrices, calculate the x and y matrices
    #note: np.multiply is ok because it is a scalar element-wise multiplication here
    for k in range(np.size(cluster_spins)):
        Sop[k][0] = np.add(Sop[k][3], Sop[k][4])/2 #the x, S-operator
        Sop[k][1] = np.add(Sop[k][4], np.multiply(-1,Sop[k][3]))/ (2j) #the y, S-operator
    
    print 'made Sx, and Sy:'+str(time.time()-t0)
    
    return Sop

def MAGCAL(scale, J, D, superexchange_matrix, single_ion_matrix, intra_dipolar_matrix, cluster_spins, Sop, T, B_field):
    import numpy as np
    import time
    TEMPTIME = time.localtime()

    print 'starting MAGCAL:'+str(TEMPTIME.tm_hour)+':'+str(TEMPTIME.tm_min)+':'+str(TEMPTIME.tm_sec)
    t0 = time.time()
       
    B_hat = np.multiply(B_field, np.sqrt(np.dot(B_field,B_field))**(-1))
    #define the units
    #SI units
    mu0=4*np.pi*10**-7 #dimensionless, permeability of free space
    muB = 9.27400915e-24 #J/T
    kB = 1.3806580*10**-23 #J/K
    NA = 6.02214129*10**23#Avagadro's number
    
    #use the spin matrices to generate the Hamiltonian
    
#    ZEEMAN
    H_zeeman = np.zeros(np.shape(Sop[0][0]))
    for i in range(np.size(cluster_spins)):
        H_zeeman = np.add( H_zeeman , np.multiply(scale*muB*B_field[0],Sop[i][0])  +  \
        np.multiply(scale*muB*B_field[1],Sop[i][1])  +  \
        np.multiply(scale*muB*B_field[2],Sop[i][2]))
        
        """here need to put in the intracluster dipolar energy, remember need to have
        g1*g2*muB^2*mu0/4/np.pi*factor, where factor is for each spin spin interaction
        and is calculated in SMMcalc_lib.intra_dipole
        """
    print 'made the zeeman matrix:'+str(time.time()-t0)


    
    H_exchange = (J*kB)*superexchange_matrix
    
    
    print 'made the superexchange matrix:'+str(time.time()-t0) 
    
    #single ion anisotropy Hamiltonian
    #use J/kB as units for D
    H_D = D*kB*single_ion_matrix

    print 'made the single-ion anisotropy matrix:'+str(time.time()-t0)
    
    #solve for the eigen vectors and eigen energies of the Hamiltonian
    (E_eig, V_eig) = np.linalg.eigh(H_zeeman  +  intra_dipolar_matrix  +  H_exchange  + H_D )
    
    print 'got the eigen-values and eigen-vectors:'+str(time.time()-t0)
    
    #each state has a z-component of the spin for each of the atoms
    #Sz[i][j] is the ith state, and the jth atom
#    Sz_ = np.zeros((np.size(E_eig), np.size(cluster_spins)) , dtype=complex)
#    
#    for i in range(np.size(E_eig)): #the ith energy state
#        for j in range(np.size(cluster_spins)): #the j'th cluster spin
#            Sz_[i][j] = np.dot(np.transpose(V_eig[:,i]).conj(),np.dot(  Sop[j][0]*B_hat[0]+Sop[j][1]*B_hat[1]+Sop[j][2]*B_hat[2]  ,V_eig[:,i]))

    Sz_ = []
    V_eig_t = np.transpose(V_eig).conj()
    for j in range(np.size(cluster_spins)): #the j'th cluster spin
        Sz_.append(  np.diagonal(     np.dot(V_eig_t,np.dot(  Sop[j][0]*B_hat[0] + Sop[j][1]*B_hat[1] + Sop[j][2]*B_hat[2]  ,V_eig))     )  )

    print 'calculated Sz for the states:'+str(time.time()-t0)    
    
    rho = np.exp(np.multiply(1/(kB*T),E_eig))/sum(np.exp(np.multiply(1/(kB*T),E_eig)))
    
    print 'calculated rho:'+str(time.time()-t0)
    
#    Sz_tot = 0
#    Sz_list = np.zeros(np.size(cluster_spins))
#    for i in range(np.size(E_eig)): #for each of the states
#        for j in range(np.size(cluster_spins)): #for each of the Sz operators of the j'th ion
#            Sz_list[j] = Sz_list[j] + Sz_[i][j]*rho[i] # the average Sz value of the j'th ion      
#            Sz_tot = Sz_tot + Sz_[i][j]*rho[i]

    Sz_tot = 0
    Sz_list = np.zeros(np.size(cluster_spins))
    for j in range(np.size(cluster_spins)): #for each of the Sz operators of the j'th ion
        Sz_list[j] = Sz_list[j] + np.dot(Sz_[j],rho) # the average Sz value of the j'th ion      
#        Sz_tot = Sz_tot + Sz_[j][i]*rho[i]

    Sz_tot = np.sum(Sz_list)

#    Sz_tot = 0
#    for i in range(np.size(Sz_list)):
#        Sz_tot = Sz_tot + Sz_list[i] # the sum of the average Sz values
    
    print 'calculated average spin:'+str(time.time()-t0)
    magnetization = Sz_tot*NA*muB*scale
    return magnetization

def MAGCALx(scale, J, D, superexchange_matrix, single_ion_matrix, intra_dipolar_matrix, cluster_spins, Sop, T, B_field):
    import numpy as np
    import time
    TEMPTIME = time.localtime()

    print 'starting MAGCAL:'+str(TEMPTIME.tm_hour)+':'+str(TEMPTIME.tm_min)+':'+str(TEMPTIME.tm_sec)
    t0 = time.time()
       
    #define the units
    #SI units
    mu0=4*np.pi*10**-7 #dimensionless, permeability of free space
    muB = 9.27400915e-24 #J/T
    kB = 1.3806580*10**-23 #J/K
    NA = 6.02214129*10**23#Avagadro's number
    
    #use the spin matrices to generate the Hamiltonian
    
#    ZEEMAN
    H_zeeman = np.zeros(np.shape(Sop[0][0]))
    for i in range(np.size(cluster_spins)):
        H_zeeman = np.add( H_zeeman , np.multiply(scale*muB*B_field,Sop[i][0]) )
        
    print 'made the zeeman matrix:'+str(time.time()-t0)
   
    H_exchange = (J*kB)*superexchange_matrix
   
    print 'made the superexchange matrix:'+str(time.time()-t0) 
    
    #single ion anisotropy Hamiltonian
    #use J/kB as units for D
    H_D = D*kB*single_ion_matrix

    print 'made the single-ion anisotropy matrix:'+str(time.time()-t0)
    
    #solve for the eigen vectors and eigen energies of the Hamiltonian
    (E_eig, V_eig) = np.linalg.eigh(H_zeeman  +  intra_dipolar_matrix  +  H_exchange  + H_D )
    
    print 'got the eigen-values and eigen-vectors:'+str(time.time()-t0)
    
    #each state has a z-component of the spin for each of the atoms
    #Sx[i][j] is the ith state, and the jth atom
    Sx_ = []
    V_eig_t = np.transpose(V_eig).conj()
    for j in range(np.size(cluster_spins)): #the j'th cluster spin
        Sx_.append(  np.diagonal(     np.dot(V_eig_t,np.dot(  Sop[j][0]  ,V_eig))     )  )

    print 'calculated Sx for the states:'+str(time.time()-t0)    
    
    rho = np.exp(np.multiply(1/(kB*T),E_eig))/sum(np.exp(np.multiply(1/(kB*T),E_eig)))
    
    print 'calculated rho:'+str(time.time()-t0)   

    Sx_tot = 0
    Sx_list = np.zeros(np.size(cluster_spins))
    for j in range(np.size(cluster_spins)): #for each of the Sx operators of the j'th ion
        Sx_list[j] = Sx_list[j] + np.dot(Sx_[j],rho) # the average Sx value of the j'th ion      

    Sx_tot = np.sum(Sx_list)
    
    print 'calculated average spin:'+str(time.time()-t0)
    magnetization = Sx_tot*NA*muB*scale
    return magnetization

def MAGCALy(scale, J, D, superexchange_matrix, single_ion_matrix, intra_dipolar_matrix, cluster_spins, Sop, T, B_field):
    import numpy as np
    import time
    TEMPTIME = time.localtime()

    print 'starting MAGCAL:'+str(TEMPTIME.tm_hour)+':'+str(TEMPTIME.tm_min)+':'+str(TEMPTIME.tm_sec)
    t0 = time.time()
       
    #define the units
    #SI units
    mu0=4*np.pi*10**-7 #dimensionless, permeability of free space
    muB = 9.27400915e-24 #J/T
    kB = 1.3806580*10**-23 #J/K
    NA = 6.02214129*10**23#Avagadro's number
    
    #use the spin matrices to generate the Hamiltonian
    
#    ZEEMAN
    H_zeeman = np.zeros(np.shape(Sop[0][0]))
    for i in range(np.size(cluster_spins)):
        H_zeeman = np.add( H_zeeman , np.multiply(scale*muB*B_field,Sop[i][1]) )
        
    print 'made the zeeman matrix:'+str(time.time()-t0)
   
    H_exchange = (J*kB)*superexchange_matrix
   
    print 'made the superexchange matrix:'+str(time.time()-t0) 
    
    #single ion anisotropy Hamiltonian
    #use J/kB as units for D
    H_D = D*kB*single_ion_matrix

    print 'made the single-ion anisotropy matrix:'+str(time.time()-t0)
    
    #solve for the eigen vectors and eigen energies of the Hamiltonian
    (E_eig, V_eig) = np.linalg.eigh(H_zeeman  +  intra_dipolar_matrix  +  H_exchange  + H_D )
    
    print 'got the eigen-values and eigen-vectors:'+str(time.time()-t0)
    
    #each state has a z-component of the spin for each of the atoms
    #Sy[i][j] is the ith state, and the jth atom
    Sy_ = []
    V_eig_t = np.transpose(V_eig).conj()
    for j in range(np.size(cluster_spins)): #the j'th cluster spin
        Sy_.append(  np.diagonal(     np.dot(V_eig_t,np.dot(  Sop[j][1]  ,V_eig))     )  )

    print 'calculated Sy for the states:'+str(time.time()-t0)    
    
    rho = np.exp(np.multiply(1/(kB*T),E_eig))/sum(np.exp(np.multiply(1/(kB*T),E_eig)))
    
    print 'calculated rho:'+str(time.time()-t0)   

    Sy_tot = 0
    Sy_list = np.zeros(np.size(cluster_spins))
    for j in range(np.size(cluster_spins)): #for each of the Sy operators of the j'th ion
        Sy_list[j] = Sy_list[j] + np.dot(Sy_[j],rho) # the average Sy value of the j'th ion      

    Sy_tot = np.sum(Sy_list)
    
    print 'calculated average spin:'+str(time.time()-t0)
    magnetization = Sy_tot*NA*muB*scale
    return magnetization

def MAGCALz(scale, J, D, superexchange_matrix, single_ion_matrix, intra_dipolar_matrix, cluster_spins, Sop, T, B_field):
    import numpy as np
    import time
    TEMPTIME = time.localtime()

    print 'starting MAGCAL:'+str(TEMPTIME.tm_hour)+':'+str(TEMPTIME.tm_min)+':'+str(TEMPTIME.tm_sec)
    t0 = time.time()
       
    #define the units
    #SI units
    mu0=4*np.pi*10**-7 #dimensionless, permeability of free space
    muB = 9.27400915e-24 #J/T
    kB = 1.3806580*10**-23 #J/K
    NA = 6.02214129*10**23#Avagadro's number
    
    #use the spin matrices to generate the Hamiltonian
    
#    ZEEMAN
    H_zeeman = np.zeros(np.shape(Sop[0][0]))
    for i in range(np.size(cluster_spins)):
        H_zeeman = np.add( H_zeeman , np.multiply(scale*muB*B_field,Sop[i][2]) )
        
    print 'made the zeeman matrix:'+str(time.time()-t0)
   
    H_exchange = (J*kB)*superexchange_matrix
   
    print 'made the superexchange matrix:'+str(time.time()-t0) 
    
    #single ion anisotropy Hamiltonian
    #use J/kB as units for D
    H_D = D*kB*single_ion_matrix

    print 'made the single-ion anisotropy matrix:'+str(time.time()-t0)
    
    #solve for the eigen vectors and eigen energies of the Hamiltonian
    (E_eig, V_eig) = np.linalg.eigh(H_zeeman  +  intra_dipolar_matrix  +  H_exchange  + H_D )
    
    print 'got the eigen-values and eigen-vectors:'+str(time.time()-t0)

    rho = np.exp(np.multiply(1/(kB*T),E_eig))/sum(np.exp(np.multiply(1/(kB*T),E_eig)))
    
    print 'calculated rho:'+str(time.time()-t0)  
    
    #each state has a z-component of the spin for each of the atoms
    #Sz[i][j] is the ith state, and the jth atom
    Sz_ = []
    V_eig_t = np.transpose(V_eig).conj()
    for j in range(np.size(cluster_spins)): #the j'th cluster spin
        Sz_.append(  np.diagonal(     np.dot(V_eig_t,np.dot(  Sop[j][2]  ,V_eig))     )  )

    print 'calculated Sz for the states:'+str(time.time()-t0)     

    Sz_tot = 0
    Sz_list = np.zeros(np.size(cluster_spins))
    for j in range(np.size(cluster_spins)): #for each of the Sz operators of the j'th ion
        Sz_list[j] = Sz_list[j] + np.dot(Sz_[j],rho) # the average Sz value of the j'th ion      

    Sz_tot = np.sum(Sz_list)
    
    print 'calculated average spin:'+str(time.time()-t0)
    magnetization = Sz_tot*NA*muB*scale
    return magnetization

def MAGCAL_powder(scale, J, D, superexchange_matrix, single_ion_matrix, intra_dipolar_matrix, cluster_spins, Sop, T, B_field):

    Mx = MAGCALx(scale, J, D, superexchange_matrix, single_ion_matrix, intra_dipolar_matrix, cluster_spins, Sop, T, B_field)
    My = MAGCALy(scale, J, D, superexchange_matrix, single_ion_matrix, intra_dipolar_matrix, cluster_spins, Sop, T, B_field)
    Mz = MAGCALz(scale, J, D, superexchange_matrix, single_ion_matrix, intra_dipolar_matrix, cluster_spins, Sop, T, B_field)
    
    Mave = (Mx+My+Mz)/3
    
    return Mave

def intra_dipolar_matrix(intra_dipolar_factor, intra_dipolar_indexer, intra_dipolar_revindexer, Sop, cluster_spins):
    import numpy as np
    #define the units
    #SI units
    mu0=4*np.pi*10**-7 #dimensionless, permeability of free space
    muB = 9.27400915e-24 #J/T
    
    num_cluster_spins = np.size(cluster_spins)
    
    H_intracluster_dip = np.zeros(np.shape(Sop[0][0]))

    for i in range(num_cluster_spins):
        for j in range(num_cluster_spins):
            for k in range(3): # spatial dimension of i'th
                for l in range(3): # spatial dimension of j'th
                    H_intracluster_dip = H_intracluster_dip +mu0* 2.0*2.0*muB**2*1*(1/np.pi)*(1/4.0)*np.dot(Sop[i][k], Sop[j][l])*intra_dipolar_factor[intra_dipolar_indexer[i][k]][intra_dipolar_indexer[j][l]]
                    
    return H_intracluster_dip


def single_ion_matrix(D_vector, Sop, cluster_spins):
    import numpy as np
    #single ion anisotropy Hamiltonian
    #use J/kB as units for D
    
    H_D = np.multiply(D_vector[0], np.dot(Sop[0][0], Sop[1][0]) )  +  \
    np.multiply(D_vector[1], np.dot(Sop[0][1], Sop[1][1]) )  +  \
    np.multiply(D_vector[2], np.dot(Sop[0][2], Sop[1][2]) )
    
    return H_D

def superexchange_matrix(Sop, cluster_spins):
    import numpy as np
    #single ion anisotropy Hamiltonian
    #use J/kB as units for D
    
    H_exchange = (  np.dot(Sop[0][0], Sop[1][0])  +  np.dot(Sop[0][1], Sop[1][1])  +  np.dot(Sop[0][2], Sop[1][2]))
    
    return H_exchange

