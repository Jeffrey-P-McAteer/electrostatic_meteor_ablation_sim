; input the potential in 1 - 3 -D .  Also set the x,y and z scales
@params_in.pro
if (n_elements(charge) eq 0) then $
  charge=readarray('charge.bin',[nx2,ny2,nz2],order=[2,1,0],/binary,skip=iskip)

chargesize=size(charge)

xv=findgen(chargesize(1))*dx*nout_avg
yv=findgen(chargesize(2))*dy*nout_avg
zv=findgen(chargesize(3))*dz*nout_avg
;Time coordinate
tv=findgen(chargesize(4))*dt*nout*iskip

