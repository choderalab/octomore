#!/bin/bash

# Sonya Hanson: February 12, 2016
# This script copies a lot of model.pdb.gz's and changes their name according to the folder they were in.
# Usage: bash autoparams_pdbs.sh /cbio/jclab/projects/parton/kinome-ensembler/models/ABL1_HUMAN_D0/ABL1_*

for f in $@
do
  NAME=$( echo "$f"|rev|cut -d / -f1|rev )
  cp $f/model.pdb.gz $NAME.pdb.gz
done
