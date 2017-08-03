#!/bin/bash
# Download and organize latest interactome3d representative dataset

# Because we don't have access to the directory, we request interactions
# and proteins from _01 to _31. There may be more in the future
wget http://interactome3d.irbbarcelona.org/user_data/human/download/representative/{interactions,proteins}_{01..31}.tgz

mkdir pdb_structures

for f in ./*.tgz;
  do tar xzvf $f -C pdb_structures;
done

# cleanup
rm ./*.tgz
