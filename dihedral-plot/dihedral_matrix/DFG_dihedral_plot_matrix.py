import matplotlib
matplotlib.use('Agg')

import mdtraj as md
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import pandas as pd

from msmbuilder import dataset

# load trajectories

#Abl_trajectories = dataset.MDTrajDataset("/cbio/jclab/projects/fah/fah-data/munged/no-solvent/10472/run0-clone0.h5")
Src_trajectories = dataset.MDTrajDataset("/cbio/jclab/projects/fah/fah-data/munged/no-solvent/10471/run0-clone0.h5")

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

Src_AD = [2190,2188,2198,2203]
Src_DF = [2200,2198,2210,2215]
Src_FG = [2212,2210,2230,2233]
Src_GL = [2233,2230,2237,2242]
Src_LA = [2239,2237,2256,2262]
Src_AR = [2258,2256,2266,2271]
Src_RL = [2268,2266,2290,2295]
Src_LI = [2292,2290,2309,2313]
Src_IG = [2311,2309,2328,2333]

Src_def_list = [('Src_AD',Src_AD),('Src_DF',Src_DF),('Src_FG',Src_FG),('Src_GL',Src_GL),('Src_LA',Src_LA),('Src_AR',Src_AR),('Src_RL',Src_RL),('Src_LI',Src_LI),('Src_IG',Src_IG)]

Src_defs = dict(Src_def_list)

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

#Making Abl Dataframe
#Abl_keys = [item[0] for item in Abl_def_list]
#Abl_dihedral_results = [dihedral_making(Abl_trajectories,Abl_defs[key]) for key in Abl_keys]
#Abl_dataframe = pd.DataFrame(data=Abl_dihedral_results,index=Abl_keys)
#Abl_dataframe = Abl_dataframe.T.astype(float)
#print Abl_dataframe

#Making Src Dataframe
Src_keys = [item[0] for item in Src_def_list]
Src_dihedral_results = [dihedral_making(Src_trajectories,Src_defs[key]) for key in Src_keys]
Src_dataframe = pd.DataFrame(data=Src_dihedral_results,index=Src_keys)
Src_dataframe = Src_dataframe.T.astype(float)
print Src_dataframe

#Abl_dataframe.to_csv('Abl_dihedrals.csv')
Src_dataframe.to_csv('Src_dihedrals.csv')

#fig, axes = plt.subplots(nrows=3, ncols=3, figsize=(12, 12))
#dataframe.plot(kind='hist', alpha=0.5)
plt.figure()
Src_dataframe.diff().hist(alpha=0.5,bins=50,color='r')

#sns.distplot(Abl_rotate, color="r",label="Abl")
###sns.distplot(Src_rotate, color="b",label="Src")
#plt.xlabel('Dihedral (radians)')
#plt.ylabel('Occupancy')
#plt.legend()

plt.savefig('src_dihedral_matrix.png',dpi=1000)

