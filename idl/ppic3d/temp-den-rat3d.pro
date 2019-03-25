; This routine calculates the temp/density ratio to evaluate from
; thermal instabilities
@den3d
@flux3d
@nvsqr3d

vx0= fluxx0/((1+den0)*n0d0)
if n_elements(fluxy0) gt 2 then vy0= fluxy0/((1+den0)*n0d0)
if n_elements(fluxz0) gt 2 then vz0= fluxz0/((1+den0)*n0d0)

vsqrx0= nvsqrx0/((1+den0)*n0d0) - vx0^2
if n_elements(nvsqry0) gt 2 then vsqry0= nvsqry0/((1+den0)*n0d0) - vy0^2
if n_elements(nvsqrz0) gt 2 then vsqrz0= nvsqrz0/((1+den0)*n0d0) - vz0^2

if ndist ge 1 then begin
  vx1= fluxx1/((1+den1)*n0d1)
  if n_elements(fluxy1) gt 2 then vy1= fluxy1/((1+den1)*n0d1)
  if n_elements(fluxz1) gt 2 then vz1= fluxz1/((1+den1)*n0d1)

  vsqrx1= nvsqrx1/((1+den1)*n0d1) - vx1^2
  if n_elements(nvsqry1) gt 2 then vsqry1= nvsqry1/((1+den1)*n0d1) - vy1^2
  if n_elements(nvsqrz1) gt 2 then vsqrz1= nvsqrz1/((1+den1)*n0d1) - vz1^2
endif

; Plot Density and temperature plots on top of each other

densize=size(den0)
!p.multi=[0,2,2]
nplot=8
iskip=round( (densize[4] -1)/(nplot))+1
for i=0.,densize[4], iskip do begin
    image_plot,reform(den0(*,*,densize[3]/2,fix(i))),$
      xv,yv,/legend, /aspect, $
      title='den0 at t='+string(i*nout*dt)

    vsqr0=reform(vsqrx0(*,*,densize[3]/2,fix(i))+ $
                vsqry0(*,*,densize[3]/2,fix(i)))
    vsqr0=vsqr0/mean(vsqr0)-1
    vsqr1=reform(vsqrx1(*,*,densize[3]/2,fix(i))+ $
                vsqry1(*,*,densize[3]/2,fix(i)))
    vsqr1=(vsqr1-mean(vsqr1))/(mean(vsqr1)+mean(vsqr0))

;    ratio=-vsqr1/den0[*,*,densize[3]/2,fix(i)]
    image_plot,vsqr1,$
      xv,yv,/legend,/aspect, $ 
      title='-Ti/(mean(Ti)+mean(Te)) at t='+string(i*nout*dt)
    
    date_plot,title
endfor

end




