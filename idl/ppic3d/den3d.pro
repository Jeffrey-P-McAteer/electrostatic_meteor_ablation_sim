; Read in densities

@params_in.pro

if (n_elements(den0) lt 2) then $
 den0=readarray('den0.bin',[nx2,ny2,nz2],order=[2,1,0],/binary,skip=iskip)

if (ndist ge 2 and n_elements(den1) lt 2) then $
 den1=readarray('den1.bin',[nx2,ny2,nz2],order=[2,1,0],/binary,skip=iskip)

if (ndist ge 3 and n_elements(den2) lt 2) then $
 den2=readarray('den2.bin',[nx2,ny2,nz2],order=[2,1,0],/binary,skip=iskip)

if (ndist ge 4 and n_elements(den3) lt 2) then $
 den3=readarray('den3.bin',[nx2,ny2,nz2],order=[2,1,0],/binary,skip=iskip)

densize=size(den0)


