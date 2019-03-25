!#/bun/csh
set datestamp = `date +%F-%I-%M`
mpirun -noallocate -partition R001 -verbose 1 -np 1024 -mode VN -cwd /project/eregion/meerso/eppic -exe /project/eregion/meerso/eppic/bin/eppic.x -args data/eppic.i > & data/eppic$datestamp.log 