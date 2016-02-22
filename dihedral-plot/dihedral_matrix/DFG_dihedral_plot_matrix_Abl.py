import matplotlib
matplotlib.use('Agg')

import mdtraj as md
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import pandas as pd

from msmbuilder import dataset

# load trajectories

Abl_trajectories = dataset.MDTrajDataset("/cbio/jclab/projects/fah/fah-data/munged/no-solvent/10472/run0-clone0.h5")
#Src_trajectories = dataset.MDTrajDataset("/cbio/jclab/projects/fah/fah-data/munged/no-solvent/10471/*.h5")

# load test trajectories

#Abl_trajectories = dataset.MDTrajDataset("../../sim-snippets/dozen_frames_abl.xtc", topology="../../sim-snippets/abl_ref.pdb")
#Src_trajectories = dataset.MDTrajDataset("../../sim-snippets/dozen_frames_src.xtc")

# define DFG dihedral ( this is from Roux umbrella sampling paper and are AlaCbeta, AlaCalpha, AspCalpha, AspCgamma)

Abl_AD = [2257,2255,2265,2270]
Abl_DF = [2267,2265,2277,2282]
Abl_FG = [2279,2277,2297,2300]
Abl_GL = [2300,2297,2304,2309]
Abl_LS = [2306,2304,2323,2330]
Abl_SR = [2325,2323,2334,2339]
Abl_RL = [2336,2334,2358,2363]
Abl_LM = [2360,2358,2377,2382]
Abl_MT = [2379,2377,2394,2400]

Abl_def_list = [('Abl_AD',Abl_AD), ('Abl_DF',Abl_DF), ('Abl_FG',Abl_FG), ('Abl_GL',Abl_GL), ('Abl_LS',Abl_LS), ('Abl_SR',Abl_SR), ('Abl_RL',Abl_RL), ('Abl_LM',Abl_LM), ('Abl_MT',Abl_MT)]

Abl_defs = dict(Abl_def_list)

Abl_DFG = [2257,2255,2265,2270]
Src_DFG = [2190,2188,2198,2203]

def dihedral_making(trajectories,dihedral_indices):

    dihedral = []

    for traj in trajectories:

        dihedral.append(md.compute_dihedrals(traj,[dihedral_indices]))

    flattened_dihedral = np.asarray([val for sublist in dihedral for val in sublist])

    #import math
    #flattened_dihedral_rotate  =  [A-(2*math.pi) if A >= 1.9 else A for A in flattened_dihedral]

    flattened_dihedral = [A for A in flattened_dihedral]

    #return flattened_dihedral_rotate
    return flattened_dihedral

keys = [item[0] for item in Abl_def_list]

dihedral_results = [dihedral_making(Abl_trajectories,Abl_defs[key]) for key in keys]

dataframe = pd.DataFrame(data=dihedral_results,index=keys)
dataframe = dataframe.T.astype(float)
print dataframe

dataframe.to_csv('Abl_dihedrals.csv')
#np.save('Src_dihedrals.npy',Src_dihedrals)

#fig, axes = plt.subplots(nrows=3, ncols=3, figsize=(12, 12))
#dataframe.plot(kind='hist', alpha=0.5)
plt.figure()
dataframe.diff().hist(alpha=0.5,bins=50)

#sns.distplot(Abl_rotate, color="r",label="Abl")
###sns.distplot(Src_rotate, color="b",label="Src")
#plt.xlabel('Dihedral (radians)')
#plt.ylabel('Occupancy')
#plt.legend()

plt.savefig('abl_dihedral_matrix.png',dpi=1000)

