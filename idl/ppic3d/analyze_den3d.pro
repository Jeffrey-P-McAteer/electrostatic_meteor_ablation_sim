;Generate density images

iskip=1

;prevent the system from reading den0 or den1
den0=[0.,0.]
den1=den0
@den3d
delvar,den0,den1

; Only analyze den2 - presumably the ion density
;For 3-D, plot the cross-section with the minimum density
if (ny gt 1) then !p.multi=[0,4,4]; if (nx ge 2*ny) then !p.multi=[0,2,4] 
densize=size(den2)

if (densize(0) ge 3) then iskip=round(densize(densize(0))/(4*8))+1

dmin=min(den2,idmin)
index_dmin=calc_index(den2,idmin)
print,'den2(',index_dmin,')= ', dmin

den2size=size(den2)
;Adjust nout to account for subcycling
den_subcycle = den_out_subcycle2

ps,'den2image-bw.ps',/landscape

;First cross-section, x-y
if (den2size(0) eq 4) then den=reform(den2(*,*,index_dmin(2),*)) $
  else den=den2 
dentitle='D2: x-y, z=' + strcompress(string(index_dmin(2)*dz*nout_avg))
.run den_image

;Add a second cross-section, x-z
if (den2size(0) eq 4) then den=reform(den2(*,index_dmin(1),*,*)) $
  else den=0
dentitle='D2: x-z, y=' + strcompress(string(index_dmin(1)*dy*nout_avg))
.run den_image

;Add a third cross-section, y-z
if (den2size(0) eq 4) then den=reform(den2(index_dmin(0),*,*,*)) $
  else den=0
dentitle='D2: y-z, x=' + strcompress(string(index_dmin(0)*dx*nout_avg))
.run den_image

stop_ps

