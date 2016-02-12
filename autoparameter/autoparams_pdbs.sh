#!/bin/bash
for f in $@
do
  NAME=$( echo "$f"|rev|cut -d / -f1|rev )
  cp $f/model.pdb.gz $NAME.pdb.gz
done
