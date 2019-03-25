; Read in densities

@eppic.i
@params_in.pro
nx3=nx2/nsubdomains

;if n_elements(ndim_space) eq 0 then ndim_space=1+(ny gt 1)+(nz gt 1)
;if ndim_space lt 3 then nz=1
;if ndim_space lt 3 then dz=0.
;if ndim_space lt 2 then ny=1
;if ndim_space lt 2 then dy=0.

;if n_elements(nsubdomains) eq 0 then nsubdomains=1 

;nx2=long(nx*nsubdomains/nout_avg) > 1
;ny2=long(ny/nout_avg) > 1
;nz2=long(nz/nout_avg) > 1

if n_elements(iskip) eq 0 then iskip = 1

nx3=nx2/nsubdomains

 den0=read_domains('den0.bin',[nx3,ny2,nz2], ndomains=nsubdomains, $
                   order=[2,1,0],/binary,skip=iskip)

if (ndist gt 1) then $
 den1=read_domains('den1.bin',[nx3,ny2,nz2], ndomains=nsubdomains, $
                   order=[2,1,0],/binary,skip=iskip)

if (ndist gt 2) then $
 den2=read_domains('den2.bin',[nx3,ny2,nz2], ndomains=nsubdomains, $
                   order=[2,1,0],/binary,skip=iskip)

if (ndist gt 3) then $
 den3=read_domains('den3.bin',[nx3,ny2,nz2], ndomains=nsubdomains, $
                   order=[2,1,0],/binary,skip=iskip)

densize=size(den0)

xv=findgen(densize(1))*dx*nout_avg
yv=findgen(densize(2))*dy*nout_avg
zv=findgen(densize(3))*dz*nout_avg
;Time coordinate
tv=findgen(densize(4))*dt*nout*iskip

; Adjust nx to contain all subdomains.
;nx=nx*nsubdomains
