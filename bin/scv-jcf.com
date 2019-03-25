
# Sample LoadLeveler Job Control File
# (10/17/05
#
### required, do not change
# @ job_type = parallel
#
### required, do not change
# @ executable = /bgl/BlueLight/ppcfloor/bglsys/bin/mpirun
#
### required, the arguments to mpirun, e.g.
# @ arguments = -verbose 1 -np 1024 -mode VN -cwd /project/eregion/meerso/eppic_runs/chromo130801/ -exe /project/eregion/meerso/eppic_dev/src/eppic.x -args /project/eregion/meerso/eppic_runs/chromo130801/eppic.i
#
### required, you must specify a runtime limit for your job
### format: hours:minutes:seconds or minutes:seconds or seconds
# @ wall_clock_limit = 600:00
#
### a directory for mpirun io files (defaults to cwd at time of submission)
### @ initialdir = 
#
### the next 3 are pathnames for mpirun's stdin, stdout, and stderr
### they can be fullpaths or relative to initialdir
### $(jobid) will be assigned by ll
# @ input = /dev/null
# @ output = eppic$(jobid).out
# @ error = eppic$(jobid).err
#
### send email when job completes
# @ notification = complete
#
### required, do not change
# @ queue
