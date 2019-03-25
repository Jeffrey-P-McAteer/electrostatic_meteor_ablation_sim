#!/bin/tcsh
rsh n01 /bin/rm -rf /scratch/meers/\* &
rsh n02 /bin/rm -rf /scratch/meers/\* &
rsh n03 /bin/rm -rf /scratch/meers/\* &
rsh n04 /bin/rm -rf /scratch/meers/\* &
rsh n05 /bin/rm -rf /scratch/meers/\* &
rsh n06 /bin/rm -rf /scratch/meers/\* &
rsh n07 /bin/rm -rf /scratch/meers/\* &
rsh n08 /bin/rm -rf /scratch/meers/\* &
rsh n09 /bin/rm -rf /scratch/meers/\* &
rsh n11 /bin/rm -rf /scratch/meers/\* &
rsh n12 /bin/rm -rf /scratch/meers/\* &
rsh n13 /bin/rm -rf /scratch/meers/\* &
rsh n14 /bin/rm -rf /scratch/meers/\* &
rsh n15 /bin/rm -rf /scratch/meers/\* &
rsh n16 /bin/rm -rf /scratch/meers/\* &
rsh n17 /bin/rm -rf /scratch/meers/\* &
rsh n18 /bin/rm -rf /scratch/meers/\* &

