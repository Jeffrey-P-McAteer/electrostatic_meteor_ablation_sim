#!/bin/bash
#SBATCH -J eppic           # job name
#SBATCH -o eppic.o%j       # output and error file name (%j expands to jobID)
#SBATCH -N 16            # total number of nodes requested
#SBATCH -n 1024            # total number of mpi tasks requested
#SBATCH -p normal         # queue (partition) -- normal, development, etc.
#SBATCH -t 30:00:00         # run time (hh:mm:ss) 
#SBATCH --mail-user=gsugar@stanford.edu
#SBATCH --mail-type=begin  # email me when the job starts
#SBATCH --mail-type=end    # email me when the job finishes

ATP_ENABLED=1
export ATP_ENABLED
export APRUN_XFER_LIMITS=1
ulimit -c unlimited

INPUTFILE=eppic.i
INPUTDIR=20181115_meteor_nob_nn_1ev
SCRATCHDIR=$SCRATCH
RUNDIR=$SCRATCHDIR
EPPIC=$HOME/eppic/src/eppic.x
JOBFILE=$RUNDIR/$INPUTDIR/$INPUTFILE

# Tar the source and the input
echo "Compressing source code and eppic.i into eppic.tgz"
cd $HOME/eppic
tar chf - src/*.h* src/*.c* src/makefile src/make.inc.knl src/make.options src/classes */*.com $JOBFILE bin |gzip - > $RUNDIR/$INPUTDIR/eppic.tgz

cd $RUNDIR/$INPUTDIR
ibrun $EPPIC $JOBFILE
echo "eppic finished."
#echo "now calling idl."
#idl analyze_meteor.pro
#echo "idl post-analysis complete"
