import matplotlib
matplotlib.use('Agg')

import matplotlib.pyplot as plt
import numpy as np
import seaborn
import mdtraj as md
import re
from glob import glob

projects = dict()

projects['ABL'] = 11400
projects['SRC'] = 11401

DFG = dict()

DFG['ABL'] = [2257,2255,2265,2270]
DFG['SRC'] = [2190,2188,2198,2203]

def DFG_dihedral_bytraj(files,def_DFG):
    
    dihedrals_with_labels = []
    
    for filename in files:
        dihedrals = []
        run = re.search('run([^-]+)',filename).group(1)
        clone = re.search('clone([^.]+)',filename).group(1)
        #eventually I will add these three pieces of information maybe I need to make a class to do this?
    
        traj = md.load(filename)
        
        dihedrals.append(md.compute_dihedrals(traj,[def_DFG]))
       
        short_filename = filename.split('/')[-1]
 
        dihedrals_with_labels.append([[short_filename,run,clone],dihedrals])
       
    return dihedrals_with_labels

for i in range(100):

    print 'working on run%s ...'%i

    files = glob('/cbio/jclab/projects/fah/fah-data/munged3/no-solvent/11400/run%s-*.h5'%i)

    dihedrals = DFG_dihedral_bytraj(files,DFG['ABL'])
    np.save('dihedrals_ABL_bytraj_run%s.npy'%i,dihedrals)

    colors = matplotlib.cm.summer(np.linspace(0, 1, len(dihedrals)))

    plt.figure(figsize=(18,14))

    frame_add = [0]
    for j in range(len(dihedrals)):
        frame_add = [frame_add[-1] + a for a in range(len(dihedrals[j][1][0]))]
        plt.plot(frame_add,dihedrals[j][1][0],'o',color=colors[j],label='%s'%dihedrals[j][0][0])
    plt.axhline(y=-0.5,color='red',label='cut-off')
    plt.legend(loc=0)
    plt.savefig('ABL_11400_DFG_by_traj_run%s.png'%i)


