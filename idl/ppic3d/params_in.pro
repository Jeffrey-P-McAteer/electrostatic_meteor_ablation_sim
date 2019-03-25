; This script sets and reads parameters from the simulation for postprocessing

print, 'inside params_in'
; Default parameters, not necessarily read in later:
ny=1
dy=1
nz=1
dz=1

if (file_info('eppic.i')).exists then $
  if (file_info('eppic.i')).mtime gt (file_info('ppic3d.i')).mtime or $
  (file_info('ppic3d.i')).size eq 0   then $
  file_copy, 'eppic.i', 'ppic3d.i', /overwrite

@ppic3d.i

if n_elements(ndim_space) eq 0 then ndim_space=1+(ny gt 1)+(nz gt 1)
if ndim_space lt 3 then nz=1
if ndim_space lt 3 then dz=0.
if ndim_space lt 2 then ny=1
if ndim_space lt 2 then dy=0.

if n_elements(nsubdomains) eq 0 then nsubdomains=1 

nx2=long(nx*nsubdomains/nout_avg) > 1
ny2=long(ny/nout_avg) > 1
nz2=long(nz/nout_avg) > 1

if n_elements(iskip) eq 0 then iskip = 1

if boundary_type_y ne 0 then ny2 = ny2+1
if boundary_type_z ne 0 then nz2 = nz2+1
