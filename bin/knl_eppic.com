#!/bin/bash
#SBATCH -J eppic           # job name
#SBATCH -o eppic.o%j       # output and error file name (%j expands to jobID)
#SBATCH -n 4096            # total number of mpi tasks requested
#SBATCH -N 64
#SBATCH -p normal          # queue (partition) -- normal, development, etc.
#SBATCH -t 01:00:00        # run time (hh:mm:ss) 
#SBATCH --mail-user=meerso@bu.edu
#SBATCH --mail-type=begin  # email me when the job starts
#SBATCH --mail-type=end    # email me when the job finishes

ATP_ENABLED=1
export ATP_ENABLED
export APRUN_XFER_LIMITS=1
ulimit -c unlimited
export GMON_OUT_PREFIX gout.*

INPUTFILE=eppic.i
INPUTDIR=Meteor3Dp160E2mV170130
SCRATCHDIR=$SCRATCH
RUNDIR=$SCRATCHDIR
EPPIC=$HOME/eppic/src/eppic.x
JOBFILE=$RUNDIR/$INPUTDIR/$INPUTFILE

cd $RUNDIR/$INPUTDIR

cp -p ~/eppic/data/eppic.tgz . $RUNDIR/$INPUTDIR/ 
ibrun $EPPIC $JOBFILE

