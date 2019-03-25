; This routine calculates pressure=nkT in space

@den3d
@flux3d
@nvsqr3d
if (n_elements(eps) eq 0) then eps=1.
;if (eps lt 1e-11) then kb=1.38E-23 else kb=1


px0= md0 * (nvsqrx0-fluxx0^2/((1+den0)*n0d0))
if n_elements(nvsqry0) gt 2 then py0= md0 * (nvsqry0-fluxy0^2/((1+den0)*n0d0))
if n_elements(nvsqrz0) gt 2 then pz0= md0 * (nvsqrz0-fluxz0^2/((1+den0)*n0d0))

if ndist ge 1 then begin
    px1= md1 * (nvsqrx1-fluxx1^2/( (1+den1)*n0d1) )
    if n_elements(nvsqry1) gt 2 then py1=md1*(nvsqry1-fluxy1^2/((1+den1)*n0d1))
    if n_elements(nvsqrz1) gt 2 then pz1=md1*(nvsqrz1-fluxz1^2/((1+den1)*n0d1))
endif

; Plot Density and Pressures

fluxx0size=size(fluxx0)
xv=findgen(fluxx0size(1))*dx*nout_avg
yv=findgen(fluxx0size(2))*dy*nout_avg
if fluxx0size(3) gt 1 then zv=findgen(fluxx0size(3))*dz*nout_avg
densize=size(den0)
!p.multi=[0,3,4]
nplot=16
iskip=round( (densize[4] -1)/(nplot))+1
for i=0.,densize[4], iskip do begin
    image_plot,reform(den0(*,*,densize[3]/2,fix(i))),$
      xv,yv,/legend, /aspect, $
      title='den0 at t='+string(i*nout*dt)

    P0=reform(sqrt(px0(*,*,densize[3]/2,fix(i))^2+ $
                   py0(*,*,densize[3]/2,fix(i))^2))
;    P0=(P0-mean(P0))/mean(P0)
    image_plot,P0,xv,yv,/legend,/aspect, $
      title='Pressure of Dist 0 at t='+string(i*nout*dt)

    P1=reform(sqrt(px1(*,*,densize[3]/2,fix(i))^2+ $
                   py1(*,*,densize[3]/2,fix(i))^2))
;    P1=(P1-mean(P1))/mean(P1)
    image_plot,P1,xv,yv,/legend,/aspect, $
      title='Pressure of Dist 1 at t='+string(i*nout*dt)

    date_plot,title
endfor

end




