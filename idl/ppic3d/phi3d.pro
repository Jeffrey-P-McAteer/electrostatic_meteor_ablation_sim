; input the potential in 1 - 3 -D .  Also set the x,y and z scales
@params_in.pro
if (n_elements(phi) lt 2) then $
  phi=readarray('phi.bin',[nx2,ny2,nz2],order=[2,1,0],/binary,skip=iskip)

phisize=size(phi)
if phisize(0) gt 1 then nphi=phisize(2) else nphi=1

xv=findgen(phisize(1))*dx*nout_avg
yv=findgen(phisize(2))*dy*nout_avg
zv=findgen(phisize(3))*dz*nout_avg
;Time coordinate
tv=findgen(phisize(4))*dt*nout*iskip

