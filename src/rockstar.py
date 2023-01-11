"""
Mass Function Class
"""

import numpy as np


#####TO READ  ROCKSTAR######

def Read_Rockstar(file):
    """
    Read in a rockstar halo catalog and return a dictionnary with all the information stored.
    R is in kpc/h. Need to have the subhalos information. Masses in Msol/h
    """

    Halo_File = []
    with open(file) as f:
        for line in f:
            Halo_File.append(line)
    a = float(Halo_File[1][4:])
    z = 1 / a - 1
    LBox = float(Halo_File[6][10:-7])
    M_part = float(Halo_File[5][16:27])
    Halo_File = Halo_File[16:]  ### Rockstar
    H_Masses, H_Radii, Host_ID, Nbr_part = [], [], [], []
    H_X, H_Y, H_Z = [], [], []
    for i in range(len(Halo_File)):
        line = Halo_File[i].split(' ')
        H_Masses.append(float(line[2]))
        H_X.append(float(line[8]))
        H_Y.append(float(line[9]))
        H_Z.append(float(line[10]))
        H_Radii.append(float(line[5]))
        Host_ID.append(float(line[-1]))
        Nbr_part.append(float(line[7]))

    H_Masses, H_X, H_Y, H_Z, H_Radii = np.array(H_Masses), np.array(H_X), np.array(H_Y), np.array(H_Z), np.array(H_Radii)
    Dict = {'M':H_Masses,'X':H_X,'Y':H_Y,'Z':H_Z, 'R':H_Radii,'z':z,'Lbox':LBox,'Host_ID':Host_ID,'Nbr_part':Nbr_part,'M_part':M_part}

    return Dict


def HMF_Rockstar(Dict,N_min,binn):
    """
    Read in the output of the previous funtion (dictionnary) with all the information stored.
    Output a binned Halo mass function.
    N_min : Minimum number of particle considered per halo
    binn  :  number of mass binn
    """
    masses  = Dict['M']
    Lbox    = Dict['Lbox']
    Host_ID = Dict['Host_ID']
    Mpart   = Dict['M_part']

    Mlist=np.logspace(np.log10(N_min*Mpart),np.log10(np.max(masses)),binn,base=10)
    dlnM=(np.log(Mlist[1])-np.log(Mlist[0]))
    
    Mh = masses[np.where(Host_ID==-1)[0]]
    indexes =np.digitize(Mh,Mlist) ### return the binn number to which each mass belong
    error   = []
    dn_dlnM = []
    for i in range (len(M__list)-1):
        dn_dlnM.append(np.count_nonzero(indexes==i+1)/(dlnM*Lbox**3))
        error.append(np.sqrt(np.count_nonzero(indexes==i+1))/(dlnM*Lbox**3))
    dn_dlnM.append(0)
    error.append(0)
    return Mlist,np.array(dn_dlnM),np.array(error)
    


