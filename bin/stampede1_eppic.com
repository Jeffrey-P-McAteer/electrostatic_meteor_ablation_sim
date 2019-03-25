#!/bin/bash
#SBATCH -J eppic           # job name
#SBATCH -o eppic.o%j       # output and error file name (%j expands to jobID)
#SBATCH -n 64            # total number of mpi tasks requested
#SBATCH -p normal         # queue (partition) -- normal, development, etc.
#SBATCH -t 2:00:00         # run time (hh:mm:ss) 
#SBATCH --mail-user=gsugar@stanford.edu
#SBATCH --mail-type=begin  # email me when the job starts
#SBATCH --mail-type=end    # email me when the job finishes

ATP_ENABLED=1
export ATP_ENABLED
export APRUN_XFER_LIMITS=1
ulimit -c unlimited

INPUTFILE=eppic.i
INPUTDIR=20170607_MeteorAblate_onecoll_colltype5
SCRATCHDIR=$SCRATCH
RUNDIR=$SCRATCHDIR
EPPIC=$HOME/eppic/src/eppic.x
JOBFILE=$RUNDIR/$INPUTDIR/$INPUTFILE

cd $RUNDIR/$INPUTDIR
ibrun $EPPIC $JOBFILE
