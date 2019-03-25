#!/bin/bash

# usage
usage() 
{
    cat <<EOF 
    bglsubmit.sh [OPTIONS] inputfile
    
    OPTIONS:
    -h                    print this help messag
    -r nrestart           Number of restart jobs to do
    -m mode               Bluegene mode to run, VN or CO 
    -n nproc              Number of processors to use
    -w wall_clock_time    Queue time limit of each job, in  minutes
    -t time_limit         Simulation Time limit, defaults to wall_clock_time

EOF
}
# options
# -- defaults
mode="VN"
nproc="64"
exec="../src/eppic"
input="eppic.i"
wallclock_limit="300"
working_dir="$PWD"
nrestart=0
time_limit=-1
# test if there is atleast the required inputfile
if test $# -gt 0; then
    # now that that is settled, check for options
    optcount=0
    while getopts "hr:m:n:w:t:" OPTION; do 
	optcount=$(( $optcount + 1 ))
	case $OPTION in
	    h)
		usage
		exit 0
		;;
	    r) 
		nrestart=$OPTARG
		optcount=$(( $optcount + 1 ))
		;;
	    m)  
		mode=$OPTARG
		optcount=$(( $optcount + 1 ))
		;;
	    n)  
		nproc=$OPTARG
		optcount=$(( $optcount + 1 ))
		;;
	    w)  
		wallclock_limit="$OPTARG:00"
		optcount=$(( $optcount + 1 ))
		;;
	    t)  
		time_limit=$OPTARG
		optcount=$(( $optcount + 1 ))
		;;

	    ?)
		usage
		exit 1
	esac
    done
    # with options checked, check again for required input file
    optcount=$(( $optcount + 1 ))
    if test $optcount -gt $#; then
        echo "Missing input file!"
	    usage
		exit 1
	else
	    eval input=\$$optcount
	fi
else
    echo "Missing input file!"
    usage
    exit 1
fi
if test $time_limit -eq -1; then
    time_limit=`echo "0.9*$wallclock_limit" | bc`
fi

# setup restart script

isrestart=`grep "iread *=" $input | awk '{ print $3 }'`
if test $isrestart -gt 0; then
  restart_input="$input"
else
  restart_input="$input.restart"
  cp $input $restart_input
  sed -i -e "s/iread *= *0/iread = 1/" $restart_input
  grep 'time_limit' $input > /dev/null
  found=$?
  if test $found -eq 0; then
      sed -i -e "s/time_limit *=.*/time_limit = $time_limit/" $input
      sed -i -e "s/time_limit *=.*/time_limit = $time_limit/" $restart_input
  else
      echo "time_limit = $time_limit" >> $input
      echo "time_limit = $time_limit" >> $restart_input
  fi
fi


# setup submission script
# -- first run
cat <<EOF
### variables custom to eppic
# @ inputfile = $input
# @ restartfile = $restart_input
# @ working_dir = $working_dir
# @ nproc = $nproc
# @ mode = $mode
# @ exec = $exec
# 
### variables requrired by bluegene for all runs
# @ job_type = parallel
# @ executable = /bgl/BlueLight/ppcfloor/bglsys/bin/mpirun
# @ initialdir = \$(working_dir)
# @ input = /dev/null
# @ output = \$(inputfile).\$(jobid).\$(step_name).out
# @ error = \$(inputfile).\$(jobid).\$(step_name).err
# @ notification = complete
#
### first job
# @ step_name = job0
# @ arguments = -verbose 1 -np \$(nproc) -mode \$(mode) -cwd \$(working_dir) -exe \$(exec) -args \$(working_dir)/\$(inputfile)
# @ wall_clock_limit = $wallclock_limit:00
# @ queue
#
EOF

# -- each restart run
if test $nrestart -gt 0; then
    echo "### restart jobs"
    restart_count=0;
    while test $restart_count -lt $nrestart; do
    
	restart_count=$(( $restart_count + 1 ))
	cat <<EOF
### restart job $restart_count
# @ step_name = job$restart_count
# @ dependency = (job$(( $restart_count - 1)) >= 0)
# @ arguments = -verbose 1 -np \$(nproc) -mode \$(mode) -cwd \$(working_dir) -exe \$(exec) -args \$(working_dir)/\$(restartfile)
# @ wall_clock_limit = $wallclock_limit:00
# @ queue
#
EOF
    done
fi
