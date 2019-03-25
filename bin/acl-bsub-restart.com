#!/bin/csh
#BSUB -N
#BSUB -n 32
#BSUB -F 1000000
#BSUB -R "span[ptile=32]" 
#BSUB -q small
#BSUB -W 11:59
#	Submit a job to the ACL SGI machines

set verbose = 1
setenv HOSTTYPE iris4d
setenv SCRATCH /scratch/meers/eppic/

# Print out environment variables:
setenv

if !(-e $SCRATCH) mkdir $SCRATCH

cd $SCRATCH
/bin/rm -r *

if !(-e $SCRATCH/data) mkdir $SCRATCH/data

#Copy restart and output data from n02
rcp -r n02:$SCRATCH $SCRATCH

cd ~/ebeam/eppic
cp -rp bin src doc lib include $SCRATCH/
cp -p data/eppic.i $SCRATCH/data/

cd $SCRATCH

echo compiling
cd src
make ../bin/eppic.x
cd ..

onintr sigerr
date
mpirun bin/eppic.x data/eppic.i >& data/eppic.log
date
#ssrun -fpcsamp bin/eppic.x data/eppic.i > data/eppic.log
echo Job Done
tail -100 data/eppic.log

exit

sigerr:
    echo "Error signal caught - no storage file made"
exit(-1)
