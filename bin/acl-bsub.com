#!/bin/csh
#BSUB -N
#BSUB -n 16
#BSUB -R "span[ptile=64]" 
#BSUB -F 1000000
#BSUB -q small
#BSUB -W 720
#	Submit a job to the ACL SGI machines


set verbose = 1

limit cputime    715m
limit coredumpsize 0

setenv HOSTTYPE iris4d
setenv SCRATCH /scratch/meers/eppic

# Print out environment variables:
#setenv

#Make essential directories
if !(-e $SCRATCH) mkdir $SCRATCH
if !(-e $SCRATCH/data) mkdir $SCRATCH/data

#delete unneeded material
cd $SCRATCH
/bin/rm -r src bin doc lib include

# Determine if this is a restart run
set ppicdir = $HOME/eppic
set infile = $ppicdir/data/eppic.i
grep "iread = 1" $infile
if !($status) then
    set restart = 1
else
    set restart = 0
endif
echo restart status $restart

#Copy restart and output data from last machine run on or n02, if restart
if ($restart) then
    if (-e $ppicdir/bin/eppic_machine) then
	set last_machine = `cat $ppicdir/bin/eppic_machine`
	echo last_machine was $last_machine 
    else
        set last_machine = "n02"
    endif
    if ($last_machine != $HOSTNAME) cp -rp /n/$last_machine/$SCRATCH/* $SCRATCH/
endif
#Set the current machine
echo $HOSTNAME > $ppicdir/bin/eppic_machine

#Copy over essentials from home directory:
cd $ppicdir
cp -rp bin src doc lib include $SCRATCH/
cp -p data/eppic.i $SCRATCH/data/
cp -p data/eppic.tgz $SCRATCH/data/

cd $SCRATCH

echo compiling
cd src
make ../bin/eppic.x
cd ..

#onintr sigerr

set logfile='data/eppic.log'

if !($restart) /bin/rm $logfile
echo Running on $HOSTNAME in `pwd`                             >>! $logfile

echo --------------------------------------------------------- >> $logfile
echo start : `date`                                            >> $logfile
echo --------------------------------------------------------- >> $logfile
    
date
mpirun bin/eppic.x data/eppic.i >>& data/eppic.log

#Check run exit status
date
#ssrun -fpcsamp bin/eppic.x data/eppic.i > data/eppic.log

echo --------------------------------------------------------- >> $logfile
echo stop : `date`                                             >> $logfile
echo --------------------------------------------------------- >> $logfile

#Look in the log file to determine if eppic completed successfully
grep "terminating all" $logfile
if !($status) then
    set ended = 1
else
    set ended = 0
endif

tail -100 data/eppic.log

#copy results to n02 (after clearing the approriate space)
#if !($HOSTNAME == "n02") then
#    rsh n02 /bin/mv $SCR/eppic $SCR/eppic-$$
#    /bin/rm -f core*
#    rcp -r $SCR/eppic n02:$SCR/
#endif

if ($ended) /bin/rm -f $ppicdir/bin/eppic_machine
if ($ended) exit

#make this into a restart run, if not already one:
if !($restart) then
    sed -e 's/iread = 0/iread = 1/' $infile >! eppic-tmp.i
    /bin/mv -f eppic-tmp.i $infile
    bsub < $ppicdir/bin/acl-bsub.com
endif



exit

sigerr:
    echo "Error signal caught"
exit(-1)
