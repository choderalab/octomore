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
#Shukla_trajectories = dataset.MDTrajDataset("/cbio/jclab/projects/shukla-trajectories/*part.pdb")

# import 2SRC structure to compare to
#SRC2 = md.load("SRC_2SRC_A.pdb")

# Define hydrogen bond coordinates (0-indexed)
#KER_src = [[35,50],[50,149]]

# Define Activation loop (resid)
#Aloop_src = [145,165]

#def shukla_coords(trajectories,KER,Aloop,SRC2):

#    difference = []
#    rmsd = []

#    for traj in trajectories:

        # append difference
#        k295e310 = md.compute_contacts(traj, [KER[0]])
#        e310r409 = md.compute_contacts(traj, [KER[1]])
#        difference.append(10*(e310r409[0] - k295e310[0])) # 10x because mdtraj is naturally in nm

        # append rmsd
#        Activation_Loop_SRC2 = SRC2.top.select("backbone and (resid %s to %s)" %(138,158))
#        Activation_Loop_kinase = traj.top.select("backbone and (resid %s to %s)" %(Aloop[0],Aloop[1]))

#        SRC2_cut = SRC2.atom_slice(Activation_Loop_SRC2)
#        traj_cut = traj.atom_slice(Activation_Loop_kinase)

#        rmsd.append(10*(md.rmsd(traj_cut,SRC2_cut,frame=0))) # 10x because mdtraj is naturaly in nm

    # flatten list of arrays
#    flattened_difference = np.asarray([val for sublist in difference for val in sublist])
#    flattened_rmsd = np.asarray([val for sublist in rmsd for val in sublist])

#    return [flattened_rmsd, flattened_difference]

# generate data
#[rmsd_src,difference_src] = shukla_coords(Shukla_trajectories,KER_src,Aloop_src,SRC2)

#save rmsd and difference data
#np.save('rmsd_shukla.npy',rmsd_src)
#np.save('difference_shukla.npy',difference_src)

rmsd_src = np.load('rmsd_shukla.npy')
difference_src = np.load('difference_shukla.npy')

#plot
#sns.kdeplot(rmsd_src,difference_src[:,0],shade=True,log=True,cmap="Greens",shade_lowest=False)
g = (sns.jointplot(rmsd_src,difference_src[:,0], kind='scatter',color="g",stat_func=None)
   .set_axis_labels("RMSD Activation Loop ($\AA$)", "d(E310-R409) - d(K295-E310) ($\AA$)"))

#plt.xlabel('RMSD Activation Loop ($\AA$)')
#plt.ylabel('d(E310-R409) - d(K295-E310) ($\AA$)')
plt.ylim(-20,20)
plt.xlim(0,10)
#plt.title('Shukla sim')
plt.tight_layout()

plt.savefig('plotting_Shukla_fig2_Shukla_sim_scatter.png',dpi=1000)
