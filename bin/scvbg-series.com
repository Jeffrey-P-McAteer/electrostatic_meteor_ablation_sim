#!/bin/csh
# This script checks to see that the queue is empty and if it is and 
# if their is an eppic.i.1 file, move that to eppic.i and run it.
set verbose = 1

setenv SCRATCH /project2/eregion/eppic/

llq |grep $USER
while (!($status))
    sleep 10s
    llq |grep $USER
end

#No job is in the queue
echo No job found in queue 

cp -fp $SCRATCH/data/eppic.i $SCRATCH/data/eppic_prev.i
mv -f $SCRATCH/data/eppic.i.1 $SCRATCH/data/eppic.i
if ($status) then
    exit
endif

#mv all the files up the queue
set i = 2
while (!($status))
    set j = i-1
    mv -f $SCRATCH/data/eppic.i.$i $SCRATCH/data/eppic.i.$j
end

llsubmit $SCRATCH/bin/scv-jcf.com

sleep 10m

#if the job is still running then reissue the command

$SCRATCH/bin/scvbg-series.com &

exit


