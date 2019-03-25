#!/bin/csh
#BSUB -P ephsphol
#	Submit a job to the SCV machines

set verbose = 1

showproj

limit cputime 299m
limit coredumpsize 0

setenv SCRATCH /project2/eregion/eppic/

# Print out environment variables:
#set
setenv

#Determine the nuber of processors from the que name
set NP=`echo $LSB_QUEUE | perl -e 'print $1 if <>=~/-mp(\d+)/'`
if ( "X$NP" == "X" ) set NP = 1
echo Running on $NP processors

cd $SCRATCH

echo compiling
cd src
gmake ../bin/eppic.x ../data/eppic.tgz
cd ..

set logfile='data/eppic.log'

# Determine if this is a restart run
set ppicdir = $SCRATCH
set infile = $ppicdir/data/eppic.i
grep "iread = 1" $infile
if !($status) then
    set restart = 1
else
    set restart = 0
endif
echo restart status $restart

if !($restart) /bin/rm $logfile
echo Running on $HOST in `pwd`                             >>! $logfile

echo --------------------------------------------------------- >> $logfile
echo start : `date`                                            >> $logfile
echo --------------------------------------------------------- >> $logfile
    
date

onintr -

if ($HOSTTYPE == "IRIX64") mpirun -np $NP bin/eppic.x data/eppic.i >>& $logfile

if ($HOSTTYPE == "AIX") poe bin/eppic.x data/eppic.i >>& $logfile -procs $NP

#set status_run = $status
#echo $stutus_run

#onintr

#Check run exit status
date
#ssrun -fpcsamp bin/eppic.x data/eppic.i > data/eppic

echo --------------------------------------------------------- >> $logfile
echo stop : `date`                                             >> $logfile
echo --------------------------------------------------------- >> $logfile

#Look in the log file to determine if eppic completed successfully
#grep "terminating all" $logfile
#if !($status) then
#    set ended = 1
#else
#    set ended = 0
#endif
##make this into a restart run, if not already one:
#if !($restart) then
#    sed -e 's/iread = 0/iread = 1/' $infile >! eppic-tmp.i
#    /bin/mv -f eppic-tmp.i $infile
#endif
#
#bsub -q $LSB_QUEUE < $SCRATCH/bin/scv-bsub.com

#if we have not completed our run, restart:
cp -fp $SCRATCH/data/eppic.i $SCRATCH/data/eppic_prev.i
mv -f $SCRATCH/data/eppic2.i $SCRATCH/data/eppic.i
if ($status) then
    exit
endif
mv -f $SCRATCH/data/eppic3.i $SCRATCH/data/eppic2.i
mv -f $SCRATCH/data/eppic4.i $SCRATCH/data/eppic3.i
mv -f $SCRATCH/data/eppic5.i $SCRATCH/data/eppic4.i
mv -f $SCRATCH/data/eppic6.i $SCRATCH/data/eppic5.i
bsub -q p4-mp16 < $SCRATCH/bin/scv-bsub.com

tail -100 $SCRATCH/data/eppic.log

exit

sigerr:
    echo "Error signal caught"
exit(-1)
