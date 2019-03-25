; Read in densities
.compile read_domains_new.pro
@eppic.i
@params_in.pro
nx3=nx2/nsubdomains

if (n_elements(istart) eq 0) then istart=0
if (n_elements(iend) eq 0) then iend=-1

sizepertime=long64((nx3*1.0*ny2))
if(ndim_space eq 3) then sizepertime = long64(sizepertime*nz2)

;print,"start should be",long64(istart*sizepertime)
;print,"end should be  ",long64(iend*sizepertime)



den0=read_domains('den0.bin',[nx3,ny2,nz2], ndomains=nsubdomains, $
                  order=[2,1,0],/binary,skip=iskip,$
                  first=long64(istart*sizepertime),$
                  last=long64(iend*sizepertime))


if (ndist gt 1) then $
 den1=read_domains('den1.bin',[nx3,ny2,nz2], ndomains=nsubdomains, $
                   order=[2,1,0],/binary,skip=iskip,$
                   first=long64(istart*sizepertime),last=long64(iend*sizepertime))

if (ndist gt 2) then $
 den2=read_domains('den2.bin',[nx3,ny2,nz2], ndomains=nsubdomains, $
                   order=[2,1,0],/binary,skip=iskip,$
                   first=long64(istart*sizepertime),last=long64(iend*sizepertime))

if (ndist gt 3) then $
 den3=read_domains('den3.bin',[nx3,ny2,nz2], ndomains=nsubdomains, $
                   order=[2,1,0],/binary,skip=iskip,$
                   first=long64(istart*sizepertime),last=long64(iend*sizepertime))

densize=size(den0)

xv=findgen(densize(1))*dx*nout_avg
yv=findgen(densize(2))*dy*nout_avg
zv=findgen(densize(3))*dz*nout_avg
;Time coordinate
tv=findgen(densize(4))*dt*nout*iskip

; Adjust nx to contain all subdomains.
nx=nx*nsubdomains
