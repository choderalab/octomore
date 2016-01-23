import matplotlib
matplotlib.use('Agg')

import mdtraj as md
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

from msmbuilder import dataset

# load trajectories

Abl_trajectories = dataset.MDTrajDataset("/cbio/jclab/projects/fah/fah-data/munged/no-solvent/10472/*.h5")
Src_trajectories = dataset.MDTrajDataset("/cbio/jclab/projects/fah/fah-data/munged/no-solvent/10471/*.h5")

Abl_trajectories_0 = dataset.MDTrajDataset("/cbio/jclab/projects/fah/fah-data/munged/no-solvent/10472/*clone0.h5")
Src_trajectories_0 = dataset.MDTrajDataset("/cbio/jclab/projects/fah/fah-data/munged/no-solvent/10471/*clone0.h5")

# define DFG dihedral ( this is from Roux umbrella sampling paper and are AlaCbeta, AlaCalpha, AspCalpha, AspCgamma)

Abl_DFG = [2257,2255,2265,2270]
Src_DFG = [2190,2188,2198,2203]

def DFG_dihedral(trajectories,def_DFG):

    dihedral = []

    for traj in trajectories:

        dihedral.append(md.compute_dihedrals(traj,[def_DFG]))

    flattened_dihedral = np.asarray([val for sublist in dihedral for val in sublist])

    return [flattened_dihedral]

Abl_firstframes = [traj[0] for traj in Abl_trajectories_0]
Src_firstframes = [traj[0] for traj in Src_trajectories_0]

#[Abl_dihedral] = DFG_dihedral(Abl_trajectories, Abl_DFG)
#[Src_dihedral] = DFG_dihedral(Src_trajectories, Src_DFG)

[Abl_lines] = DFG_dihedral(Abl_firstframes, Abl_DFG)
[Src_lines] = DFG_dihedral(Src_firstframes, Src_DFG)

#np.save('Abl_dihedral_newsims.npy',Abl_dihedral)
#np.save('Src_dihedral_newsims.npy',Src_dihedral)

Abl_dihedral = np.load('../repeat/Abl_dihedral.npy')
Src_dihedral = np.load('../repeat/Src_dihedral.npy')

import math

Abl_rotate =  [A-(2*math.pi) if A >= 1.9 else A for A in Abl_dihedral]
Src_rotate =  [S-(2*math.pi) if S >= 1.9 else S for S in Src_dihedral]

Abl_line_rotate =  [A-(2*math.pi) if A >= 1.9 else A for A in Abl_lines]
Src_line_rotate =  [S-(2*math.pi) if S >= 1.9 else S for S in Src_lines]

sns.distplot(Abl_rotate, color="r",label="Abl")
sns.distplot(Src_rotate, color="b",label="Src")

Abl_line_rotate = np.asarray(Abl_line_rotate)
Src_line_rotate = np.asarray(Src_line_rotate)

print Abl_line_rotate

#plt.axvline(Abl_line_rotate.all(), color="r")
#plt.axvline(Src_line_rotate.all(), color="b")

for i in range(len(Abl_line_rotate)):
	plt.axvline(Abl_line_rotate[i], color="r")
for i in range(len(Src_line_rotate)):
	plt.axvline(Src_line_rotate[i], color="b")

plt.xlabel('Dihedral (radians)')
plt.ylabel('Occupancy')
plt.legend()

plt.savefig('abl_src_DFG_dihedral_hist_oldsims_lines.png',dpi=1000)

