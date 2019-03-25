; input the potential in 1 - 3 -D .  Also set the x,y and z scales
@eppic.i
@params_in.pro

if n_elements(boundary_type) eq 0 then boundary_type=0
if (boundary_type ne 0) then padleft=1
if (boundary_type ne 0) then padright=2

if (nout_avg gt 1) then padleft=0
if (nout_avg gt 1) then padright=0

if n_elements(istart) eq 0 then istart=0
if (n_elements(iend) eq 0) then iend=-1

if (n_elements(padleft) eq 0) then padleft=0
if (n_elements(padright) eq 0) then padright = 0

array_size=(nx2/nsubdomains+padleft+padright)*ny2*nz2 ; Size of a single array (on 1 domain)

phi=read_domains('phi.bin',[nx2/nsubdomains,ny2,nz2], ndomains=nsubdomains, $
                 order=[2,1,0],/binary,skip=iskip,first=long64(istart)*long64(array_size), $
                 last=long64(iend)*long64(array_size),non_periodic_x=[padleft,padright]);,/double)
phisize=size(phi)

; Adjust nx to contain all subdomains.
nx=nx*nsubdomains

if phisize(0) gt 1 then nphi=phisize(2) else nphi=1

xv=findgen(phisize(1))*dx*nout_avg
yv=findgen(phisize(2))*dy*nout_avg
zv=findgen(phisize(3))*dz*nout_avg
;Time coordinate
tv=findgen(phisize(4))*dt*nout*iskip

