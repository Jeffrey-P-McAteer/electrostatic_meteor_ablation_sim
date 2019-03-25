; Read in nvsqrsities
@params_in.pro

if (n_elements(nvsqrx0) lt 2) then $
 nvsqrx0=readarray('nvsqrx0.bin',[nx2,ny2,nz2],order=[2,1,0],/binary,skip=iskip)
if (n_elements(nvsqry0) lt 2) then $
 nvsqry0=readarray('nvsqry0.bin',[nx2,ny2,nz2],order=[2,1,0],/binary,skip=iskip)
if (n_elements(nvsqrz0) lt 2) then $
 nvsqrz0=readarray('nvsqrz0.bin',[nx2,ny2,nz2],order=[2,1,0],/binary,skip=iskip)


if (ndist ge 1 and n_elements(nvsqrx1) lt 2) then $
 nvsqrx1=readarray('nvsqrx1.bin',[nx2,ny2,nz2],order=[2,1,0],/binary,skip=iskip)
if (ndist ge 1 and n_elements(nvsqry1) lt 2) then $
 nvsqry1=readarray('nvsqry1.bin',[nx2,ny2,nz2],order=[2,1,0],/binary,skip=iskip)
if (ndist ge 1 and n_elements(nvsqrz1) lt 2) then $
 nvsqrz1=readarray('nvsqrz1.bin',[nx2,ny2,nz2],order=[2,1,0],/binary,skip=iskip)

if (ndist ge 2 and n_elements(nvsqrx2) lt 2) then $
 nvsqrx2=readarray('nvsqrx2.bin',[nx2,ny2,nz2],order=[2,1,0],/binary,skip=iskip)
if (ndist ge 2 and n_elements(nvsqry2) lt 2) then $
 nvsqry2=readarray('nvsqry2.bin',[nx2,ny2,nz2],order=[2,1,0],/binary,skip=iskip)
if (ndist ge 2 and n_elements(nvsqrz2) lt 2) then $
 nvsqrz2=readarray('nvsqrz2.bin',[nx2,ny2,nz2],order=[2,1,0],/binary,skip=iskip)


