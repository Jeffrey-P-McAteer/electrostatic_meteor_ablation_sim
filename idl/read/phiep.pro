; input the potential in 1 - 3 -D .  Also set the x,y and z scales
@eppic.i
@params_in.pro

array_size=nx2/nsubdomains*ny2*nz2 ; Size of a single array (on 1 domain)
if n_elements(istart) eq 0 then istart=0
if n_elements(istop) eq 0 then istop=-1

phi=read_domains('phi.bin',[nx2/nsubdomains,ny2,nz2], ndomains=nsubdomains, $
                 order=[2,1,0],/binary,skip=iskip,first=istart*array_size, $
                 last=istop*array_size)
phisize=size(phi)

; Adjust nx to contain all subdomains.
nx=nx*nsubdomains

if phisize(0) gt 1 then nphi=phisize(2) else nphi=1

xv=findgen(phisize(1))*dx*nout_avg
yv=findgen(phisize(2))*dy*nout_avg
zv=findgen(phisize(3))*dz*nout_avg
;Time coordinate
tv=findgen(phisize(4))*dt*nout*iskip

