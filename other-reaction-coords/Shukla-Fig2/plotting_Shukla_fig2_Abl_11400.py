# import libraries
import matplotlib
matplotlib.use('Agg')

import mdtraj as md
import matplotlib.pyplot as plt
import numpy as np

from msmbuilder import dataset

import seaborn as sns
sns.set_style("whitegrid")
sns.set_context("poster")

# load trajectories
Abl_trajectories = dataset.MDTrajDataset("/cbio/jclab/projects/fah/fah-data/munged/no-solvent/11400/*.h5")

# import 2SRC structure to compare to
SRC2 = md.load("ABL_2SRC_A.pdb")

# Define hydrogen bond coordinates (0-indexed)
KER_abl = [[29,44],[44,144]]

# Define Activation loop (resid)
Aloop_abl = [140,160]

def shukla_coords(trajectories,KER,Aloop,SRC2):

    difference = []
    rmsd = []

    for traj in trajectories:

        # append difference
        k295e310 = md.compute_contacts(traj, [KER[0]])
        e310r409 = md.compute_contacts(traj, [KER[1]])
        difference.append(10*(e310r409[0] - k295e310[0])) # 10x because mdtraj is naturally in nm

        # append rmsd
        Activation_Loop_SRC2 = SRC2.top.select("backbone and (resid %s to %s)" %(140,160))
        Activation_Loop_kinase = traj.top.select("backbone and (resid %s to %s)" %(Aloop[0],Aloop[1]))

        SRC2_cut = SRC2.atom_slice(Activation_Loop_SRC2)
        traj_cut = traj.atom_slice(Activation_Loop_kinase)

        rmsd.append(10*(md.rmsd(traj_cut,SRC2_cut,frame=0))) # 10x because mdtraj is naturaly in nm

    # flatten list of arrays
    flattened_difference = np.asarray([val for sublist in difference for val in sublist])
    flattened_rmsd = np.asarray([val for sublist in rmsd for val in sublist])

    return [flattened_rmsd, flattened_difference]

# generate data
[rmsd_abl,difference_abl] = shukla_coords(Abl_trajectories,KER_abl,Aloop_abl,SRC2)

#save rmsd and difference data
np.save('rmsd_abl.npy',rmsd_abl)
np.save('difference_abl.npy',difference_abl)

#plot
sns.kdeplot(rmsd_abl,difference_abl[:,0],shade=True,log=True,cmap="Blues",shade_lowest=False,cut=10)

plt.xlabel('RMSD Activation Loop ($\AA$)')
plt.ylabel('d(E286-R386) - d(K271-E286) ($\AA$)')
plt.ylim(-20,20)
plt.xlim(0,10)
plt.title('Abl sims - 11400')

plt.savefig('plotting_Shukla_fig2_Abl_11400.png',dpi=1000)
