import matplotlib
matplotlib.use('Agg')

import mdtraj as md
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

from msmbuilder import dataset

# load trajectories

Abl_in = dataset.MDTrajDataset("Abl_DFG_in/ABL1*.pdb")
Abl_out = dataset.MDTrajDataset("Abl_DFG_out/ABL1*.pdb")

Src_in = dataset.MDTrajDataset("Src_DFG_in/SRC*.pdb")
Src_out = dataset.MDTrajDataset("Src_DFG_out/SRC*.pdb")

# define DFG dihedral ( this is from Roux umbrella sampling paper and are AlaCbeta, AlaCalpha, AspCalpha, AspCgamma)
#These are with hydrogens
#Abl_DFG = [2257,2255,2265,2270]
#Src_DFG = [2190,2188,2198,2203]
# Below are the dihedral coordinates in PDBs without hydrogens
Abl_DFG = [1117,1116,1121,1123]
Src_DFG = [1074,1073,1078,1080]

def DFG_dihedral(trajectories,def_DFG):

    dihedral = []

    for traj in trajectories:

        dihedral.append(md.compute_dihedrals(traj,[def_DFG]))

    flattened_dihedral = np.asarray([val for sublist in dihedral for val in sublist])

    return [flattened_dihedral]

[Abl_in_dihedrals] = DFG_dihedral(Abl_in, Abl_DFG)
[Abl_out_dihedrals] = DFG_dihedral(Abl_out, Abl_DFG)
[Src_in_dihedrals] = DFG_dihedral(Src_in, Src_DFG)
[Src_out_dihedrals] = DFG_dihedral(Src_out, Src_DFG)

print len(Abl_in)

import math

#def dih_rotate(dihedral):
#     dihedral_rotate = [A-(2*math.pi) if A >= 1.9 else A for A in dihedral]
#     return dihedral_rotate

Abl_in_rotate =  [A-(2*math.pi) if A >= 1.9 else A for A in Abl_in_dihedrals]
Abl_out_rotate =  [A-(2*math.pi) if A >= 1.9 else A for A in Abl_out_dihedrals]
Src_in_rotate =  [S-(2*math.pi) if S >= 1.9 else S for S in Src_in_dihedrals]
Src_out_rotate =  [S-(2*math.pi) if S >= 1.9 else S for S in Src_out_dihedrals]

#sns.distplot(Abl_rotate, color="r",label="Abl")
#sns.distplot(Src_rotate, color="b",label="Src")

for i in range(len(Abl_in_rotate)):
	plt.axvline(Abl_in_rotate[i], color="r")
for i in range(len(Abl_out_rotate)):
        plt.axvline(Abl_out_rotate[i], color="m")
for i in range(len(Src_in_rotate)):
	plt.axvline(Src_in_rotate[i], color="b")
for i in range(len(Src_out_rotate)):
        plt.axvline(Src_out_rotate[i], color="g")

plt.xlabel('Dihedral (radians)')
plt.ylabel('Occupancy')
#plt.legend()

plt.savefig('DFG_dihedral_plot_autoparameter_lines.png')

