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

# define DFG dihedral ( this is shifted one over from Roux umbrella sampling paper and are AspCbeta, AspCalpha, PheCalpha, PheCgamma instead of AlaCbeta, AlaCalpha, AspCalpha, AspCgamma)

# You can select these residues either by looking at the PDB or
#confirm that you've picked the right atoms using this mdtraj syntax:
# traj = md.load('/cbio/jclab/projects/fah/fah-data/munged/no-solvent/10471/run0-clone0.h5')
# topology = traj.topology
#atom = topology.atom(2198)
#print('''Hi! I am the %sth atom, and my name is %s. 
#I am a %s atom with %s bonds. 
#I am part of an %s residue.''' % ( atom.index, atom.name, atom.element.name, atom.n_bonds, atom.residue.name))

#These are the AD atom numbers
#Abl_DFG = [2257,2255,2265,2270]
#Src_DFG = [2190,2188,2198,2203]

#These are the DF atom numbers
Abl_DFG = [2267,2265,2277,2282]
Src_DFG = [2200,2198,2210,2215]

def DFG_dihedral(trajectories,def_DFG):

    dihedral = []

    for traj in trajectories:

        dihedral.append(md.compute_dihedrals(traj,[def_DFG]))

    flattened_dihedral = np.asarray([val for sublist in dihedral for val in sublist])

    return [flattened_dihedral]

[Abl_dihedral] = DFG_dihedral(Abl_trajectories, Abl_DFG)
[Src_dihedral] = DFG_dihedral(Src_trajectories, Src_DFG)

np.save('Abl_dihedral_nowDF.npy',Abl_dihedral)
np.save('Src_dihedral_nowDF.npy',Src_dihedral)

import math

Abl_rotate =  [A-(2*math.pi) if A >= 1.9 else A for A in Abl_dihedral]
Src_rotate =  [S-(2*math.pi) if S >= 1.9 else S for S in Src_dihedral]

sns.distplot(Abl_rotate, color="r",label="Abl")
sns.distplot(Src_rotate, color="b",label="Src")
plt.xlabel('Dihedral (radians)')
plt.ylabel('Occupancy')
plt.legend()

plt.savefig('abl_src_DFG_dihedral_hist_nowDF.png',dpi=1000)

