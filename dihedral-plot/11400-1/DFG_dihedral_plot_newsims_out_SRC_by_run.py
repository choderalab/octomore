import matplotlib
matplotlib.use('Agg')

import mdtraj as md
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

from msmbuilder import dataset

# load clone0 of each trajectory

Src_trajectories_0 = dataset.MDTrajDataset("/cbio/jclab/projects/fah/fah-data/munged2/no-solvent/11401/*clone0.h5")

# define DFG dihedral (this is from Roux umbrella sampling paper and are AlaCbeta, AlaCalpha, AspCalpha, AspCgamma)

Src_DFG = [2190,2188,2198,2203]

Src_runs = len(Src_trajectories_0)

Src_project = 11401

Src_dihedral = np.load('Src_dihedral_newsims_separate.npy')

#funciton
def DFG_dihedral(trajectories,def_DFG):

    dihedral = []

    for traj in trajectories:

        dihedral.append(md.compute_dihedrals(traj,[def_DFG]))

    flattened_dihedral = np.asarray([val for sublist in dihedral for val in sublist])

    return [flattened_dihedral]

# Calculate dihedral at the first frame of clone0 of each run

Src_firstframes = [traj[0] for traj in Src_trajectories_0]

[Src_lines] = DFG_dihedral(Src_firstframes, Src_DFG)

# Rotate dihedral so histogram doesn't get cut in figure
import math

Src_rotate = [ np.array( [S-(2*math.pi) if S >= 1.9 else S for S in run ] ) for run in Src_dihedral[0] ]

Src_line_rotate =  [S-(2*math.pi) if S >= 1.9 else S for S in Src_lines]

# Define which sims start in DFG-in vs DFG-out conformation
Src_line_rotate = np.asarray(Src_line_rotate)

Src_DFG_in = np.where(Src_line_rotate > -0.5)
Src_DFG_out = np.where(Src_line_rotate < -0.5)

print Src_DFG_out

# Define 8 src_out runs.
#src_out = np.vstack([ Src_rotate[index] for index in Src_DFG_out[0]])
src_out_0 = np.vstack([ Src_rotate[Src_DFG_out[0][0]] ])
src_out_1 = np.vstack([ Src_rotate[Src_DFG_out[0][1]] ])
src_out_2 = np.vstack([ Src_rotate[Src_DFG_out[0][2]] ])
src_out_3 = np.vstack([ Src_rotate[Src_DFG_out[0][3]] ])
src_out_4 = np.vstack([ Src_rotate[Src_DFG_out[0][4]] ])
src_out_5 = np.vstack([ Src_rotate[Src_DFG_out[0][5]] ])
src_out_6 = np.vstack([ Src_rotate[Src_DFG_out[0][6]] ])
src_out_7 = np.vstack([ Src_rotate[Src_DFG_out[0][7]] ])

#Plot just vertical lines for Src DFG out

Src_line_rotate_out = np.vstack([ Src_line_rotate[index] for index in Src_DFG_out[0]])

plt.axvline(Src_line_rotate_out[0], color="g", ymin=0.95)
plt.axvline(Src_line_rotate_out[1], color="m", ymin=0.95)
plt.axvline(Src_line_rotate_out[2], color="b", ymin=0.95)
plt.axvline(Src_line_rotate_out[3], color="r", ymin=0.95)
plt.axvline(Src_line_rotate_out[4], color="c", ymin=0.95)
plt.axvline(Src_line_rotate_out[5], color="y", ymin=0.95)
plt.axvline(Src_line_rotate_out[6], color="k", ymin=0.95)
plt.axvline(Src_line_rotate_out[7], color="w", ymin=0.95)

# Plot histogram with special seaborn sauce
sns.distplot(src_out_0, color="g",label="Src - out (%s)" %Src_DFG_out[0][0] )
sns.distplot(src_out_1, color="m",label="Src - out (%s)" %Src_DFG_out[0][1] )
sns.distplot(src_out_2, color="b",label="Src - out (%s)" %Src_DFG_out[0][2] )
sns.distplot(src_out_3, color="r",label="Src - out (%s)" %Src_DFG_out[0][3] )
sns.distplot(src_out_4, color="c",label="Src - out (%s)" %Src_DFG_out[0][4] )
sns.distplot(src_out_5, color="y",label="Src - out (%s)" %Src_DFG_out[0][5] )
sns.distplot(src_out_6, color="k",label="Src - out (%s)" %Src_DFG_out[0][6] )
sns.distplot(src_out_7, color="w",label="Src - out (%s)" %Src_DFG_out[0][7] )

plt.xlabel('Dihedral (radians)')
plt.ylabel('Occupancy')
plt.legend()

plt.savefig('src_DFG_dihedral_hist_newsims_out_byrun.png',dpi=1000)

