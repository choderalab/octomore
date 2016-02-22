# Sonya Hanson - February 17, 2016

# This makes a pymol alignment file for Src and Abl structures so it's
# easier to manually decide if they're DFG_in or DFG_out.

# Usage: python align.py SRC* 

import sys

# Write pymol alignment file.
# Usage  'open -a MacPyMOL *.pdb' and then @SRC_align.pml in pymol.
structures = sys.argv[1:]

names = []

for struc in structures:
   names.append(struc.split('.')[0])

input = names[0].split('_')[0]

# Defiine what structure we're aligning to.

ref = names[0]

# Make and write file.

f=open('%s_align.pml' % input,'w')

for name in names[1:]:
    f.write ('align %s, %s ;\n' % (name,ref))

f.close()
