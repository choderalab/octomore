import matplotlib
matplotlib.use('Agg')

import mdtraj as md
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

from msmbuilder import dataset

# load clone0 of each trajectory

Abl_trajectories_0 = dataset.MDTrajDataset("/cbio/jclab/projects/fah/fah-data/munged2/no-solvent/11400/*clone0.h5")
Src_trajectories_0 = dataset.MDTrajDataset("/cbio/jclab/projects/fah/fah-data/munged2/no-solvent/11401/*clone0.h5")

# define DFG dihedral (this is from Roux umbrella sampling paper and are AlaCbeta, AlaCalpha, AspCalpha, AspCgamma)

Abl_DFG = [2257,2255,2265,2270]
Src_DFG = [2190,2188,2198,2203]

Abl_runs = len(Abl_trajectories_0)
Src_runs = len(Src_trajectories_0)

Abl_project = 11400
Src_project = 11401

def DFG_dihedral(trajectories,def_DFG):

    dihedral = []

    for traj in trajectories:

        dihedral.append(md.compute_dihedrals(traj,[def_DFG]))

    flattened_dihedral = np.asarray([val for sublist in dihedral for val in sublist])

    return [flattened_dihedral]

def DFG_dihedral_byrun(project,runs,def_DFG):
    
    dihedral = []
    dihedral_combinetrajs = []
    print "Working on project %s." % project

    for run in range(runs):
     
        trajectories = dataset.MDTrajDataset("/cbio/jclab/projects/fah/fah-data/munged2/no-solvent/%d/run%d-clone*.h5" % (project,run))
        print "Run %s has %s trajectories." % (run,len(trajectories))        

        for traj in trajectories:

            dihedral_combinetrajs.append(md.compute_dihedrals(traj,[def_DFG]))
        # flatten
        dihedral_combinetrajs = [val for sublist in dihedral_combinetrajs for val in sublist]

        dihedral.append(dihedral_combinetrajs) 
        dihedral_combinetrajs = []

    dihedral = np.asarray([dihedral])

    return [dihedral]

#[Abl_dihedral] = DFG_dihedral_byrun(Abl_project, Abl_runs, Abl_DFG)
#[Src_dihedral] = DFG_dihedral_byrun(Src_project, Src_runs, Src_DFG)

#np.save('Abl_dihedral_newsims_separate.npy',Abl_dihedral)
#np.save('Src_dihedral_newsims_separate.npy',Src_dihedral)

Abl_dihedral = np.load('Abl_dihedral_newsims_separate.npy')
Src_dihedral = np.load('Src_dihedral_newsims_separate.npy')

# Calculate dihedral at the first frame of clone0 of each run

Abl_firstframes = [traj[0] for traj in Abl_trajectories_0]
Src_firstframes = [traj[0] for traj in Src_trajectories_0]

[Abl_lines] = DFG_dihedral(Abl_firstframes, Abl_DFG)
[Src_lines] = DFG_dihedral(Src_firstframes, Src_DFG)

# Rotate dihedral so histogram doesn't get cut in figure
import math

Abl_rotate = [ np.array( [A-(2*math.pi) if A >= 1.9 else A for A in run ] ) for run in Abl_dihedral[0] ]
Src_rotate = [ np.array( [S-(2*math.pi) if S >= 1.9 else S for S in run ] ) for run in Src_dihedral[0] ]

Abl_line_rotate =  [A-(2*math.pi) if A >= 1.9 else A for A in Abl_lines]
Src_line_rotate =  [S-(2*math.pi) if S >= 1.9 else S for S in Src_lines]

# Define which sims start in DFG-in vs DFG-out conformation
Abl_line_rotate = np.asarray(Abl_line_rotate)
Src_line_rotate = np.asarray(Src_line_rotate)

Abl_DFG_in = np.where(Abl_line_rotate > -0.5)
Abl_DFG_out = np.where(Abl_line_rotate < -0.5)

Src_DFG_in = np.where(Src_line_rotate > -0.5)
Src_DFG_out = np.where(Src_line_rotate < -0.5)

# Accumulate these in a loop.
abl_in = np.vstack([ Abl_rotate[index] for index in Abl_DFG_in[0]])
src_in = np.vstack([ Src_rotate[index] for index in Src_DFG_in[0]])
abl_out = np.vstack([ Abl_rotate[index] for index in Abl_DFG_out[0]])
src_out = np.vstack([ Src_rotate[index] for index in Src_DFG_out[0]])

# Plot vertical lines at dihedral for first frames
#for i in range(len(Abl_line_rotate)):
#        plt.axvline(Abl_line_rotate[i], color="r")
#for i in range(len(Src_line_rotate)):
#        plt.axvline(Src_line_rotate[i], color="b")
#plt.vline can set y range.

#Plot just vertical lines for Src DFG out

Src_line_rotate_out = np.vstack([ Src_line_rotate[index] for index in Src_DFG_out[0]])
Src_line_rotate_in = np.vstack([ Src_line_rotate[index] for index in Src_DFG_in[0]])

for i in range(len(Src_line_rotate_in)):
        plt.axvline(Src_line_rotate_in[i], color="b", ymin=0.95)
for i in range(len(Src_line_rotate_out)):
        plt.axvline(Src_line_rotate_out[i], color="g", ymin=0.95)

# Plot histogram with special seaborn sauce
#sns.distplot(abl_in, color="r",label="Abl - in (%s) " %len(Abl_DFG_in[0]) )
#sns.distplot(abl_out, color="m",label="Abl - out (%s) " %len(Abl_DFG_out[0]) )
sns.distplot(src_in, color="b",label="Src - in (%s) " %len(Src_DFG_in[0]) )
sns.distplot(src_out, color="g",label="Src - out (%s)" %len(Src_DFG_out[0]) )

plt.xlabel('Dihedral (radians)')
plt.ylabel('Occupancy')
plt.legend()

plt.savefig('src_DFG_dihedral_hist_newsims_NYAS.png',dpi=1000)

