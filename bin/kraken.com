#!/bin/bash
#PBS -A TG-ATM110034
#PBS -l size=4104,walltime=14:00:00
#PBS -j oe
#PBS -m e
#PBS -M meerso@bu.edu

ATP_ENABLED=1
export ATP_ENABLED
export APRUN_XFER_LIMITS=1
ulimit -c unlimited

INPUTFILE=FBI3D.in
INPUTDIR=FBI3D2048E105nu.5nuiSqrtV120704

SCRATCHDIR=/lustre/scratch/moppenhe
RUNDIR=$SCRATCHDIR/eppic_runs
EPPIC=$SCRATCHDIR/eppic/src/eppic.x
JOBFILE=$RUNDIR/$INPUTDIR/$INPUTFILE


cd $PBS_O_WORKDIR
aprun -n 4096 $EPPIC $JOBFILE
