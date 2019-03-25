

@eppic.i

if (n_elements(ny) eq 0) then nz = 1
if (n_elements(nz) eq 0) then nz = 1

if (n_elements(nout_avg) eq 0) then nout_avg =1
print,nx,ny,nz,nout_avg,nx/nout_avg
nx_out=LONG64(nx/nout_avg)
if (nx_out lt 1) then nx_out = 1
ny_out=LONG64(ny/nout_avg)
if (ny_out lt 1) then ny_out = 1
nz_out=LONG64(nz/nout_avg)
if (nz_out lt 1) then nz_out = 1
ncells = nx_out*ny_out*nz_out
print,ncells,nx_out,ny_out,nz_out

;.compile timesteps
totaltime = timesteps("den0.bin",ncells,nsubdomains,double=True)
