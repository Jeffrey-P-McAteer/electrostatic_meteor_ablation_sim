; input the potential in 1 - 3 -D .  Also set the x,y and z scales
@eppic.i
@params_in.pro

array_size=nx2/nsubdomains*ny2*nz2 ; Size of a single array (on 1 domain)
if n_elements(istart) eq 0 then istart=0
if n_elements(istop) eq 0 then istop=-1

phistore=read_domains('phistore1.bin',[nx2/nsubdomains,ny2,nz2], ndomains=nsubdomains, $
                 order=[2,1,0],/binary,skip=iskip,first=istart*array_size, $
                 last=istop*array_size)
phistore_size=size(phistore)

; Adjust nx to contain all subdomains.
nx=nx*nsubdomains

if phistore_size(0) gt 1 then nphistore=phistore_size(2) else nphistore=1

xv=findgen(phistore_size(1))*dx*nout_avg
yv=findgen(phistore_size(2))*dy*nout_avg
zv=findgen(phistore_size(3))*dz*nout_avg
;Time coordinate
tv=findgen(phistore_size(4))*dt*nout*iskip

