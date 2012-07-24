# -*- coding: utf-8 -*-
"""
Created on Thu Jan 26 10:43:46 2012

@author: -
"""

MvsT_calc = []
MvsT_output = []
T_calc = np.concatenate((np.linspace(2,10,9),np.linspace(15,300,58)))
for T__ in T_calc:
    M__ = SMMcalc_lib.MAGCAL_powder(xm[0],xm[1],xm[2], superexchange_matrix, single_ion_matrix, intra_dipolar_matrix, cluster_spins, Sop, T__, 0.1)
    MvsT_calc.append(M__)
    MvsT_output.append([T__,M__])

np.savetxt('MvsT_calc_ZQ49_022012.txt', MvsT_output)

MvsH_calc = []
MvsH_output = []
B_calc = np.concatenate((np.linspace(.1,1,10),np.linspace(1.5,7,12)))
for B__ in B_calc:
    M__ = SMMcalc_lib.MAGCAL_powder(xm[0],xm[1],xm[2], superexchange_matrix, single_ion_matrix, intra_dipolar_matrix, cluster_spins, Sop, 2, B__)
    MvsH_calc.append(M__)
    MvsH_output.append([B__, M__])

np.savetxt('MvsH_calc_ZQ49_022012.txt', MvsH_output)
