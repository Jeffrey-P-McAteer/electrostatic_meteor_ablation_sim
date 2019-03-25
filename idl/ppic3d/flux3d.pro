; Read in fluxsities
@params_in.pro

if n_elements(reread) eq 0 then reread=1

if (reread or n_elements(fluxx0) lt 2) then $
 fluxx0=readarray('fluxx0.bin',[nx2,ny2,nz2],order=[2,1,0],/binary,skip=iskip)
if (reread or n_elements(fluxy0) lt 2) then $
 fluxy0=readarray('fluxy0.bin',[nx2,ny2,nz2],order=[2,1,0],/binary,skip=iskip)
if (reread or n_elements(fluxz0) lt 2) then $
 fluxz0=readarray('fluxz0.bin',[nx2,ny2,nz2],order=[2,1,0],/binary,skip=iskip)


if ((reread or n_elements(fluxx1) lt 2) and ndist ge 1) then $
 fluxx1=readarray('fluxx1.bin',[nx2,ny2,nz2],order=[2,1,0],/binary,skip=iskip)
if ((reread or n_elements(fluxy1) lt 2) and ndist ge 1) then $
 fluxy1=readarray('fluxy1.bin',[nx2,ny2,nz2],order=[2,1,0],/binary,skip=iskip)
if ((reread or n_elements(fluxz1) lt 2) and ndist ge 1) then $
 fluxz1=readarray('fluxz1.bin',[nx2,ny2,nz2],order=[2,1,0],/binary,skip=iskip)

if ((reread or n_elements(fluxz2) lt 2) and ndist ge 2) then $
 fluxx2=readarray('fluxx2.bin',[nx2,ny2,nz2],order=[2,1,0],/binary,skip=iskip)
if ((reread or n_elements(fluxy2) lt 2) and ndist ge 2) then $
 fluxy2=readarray('fluxy2.bin',[nx2,ny2,nz2],order=[2,1,0],/binary,skip=iskip)
if ((reread or n_elements(fluxz2) lt 2) and ndist ge 2) then $
 fluxz2=readarray('fluxz2.bin',[nx2,ny2,nz2],order=[2,1,0],/binary,skip=iskip)


