; This routine calculates v^2 (temperatures) in space
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

; Only works in 2D & 3D:
if ny2 gt 1 then begin
    !p.multi=[0,3,4]
    nplot=16
    iskip=round( (densize[4] -1)/(nplot))+1
    for i=0.,densize[4], iskip do begin
        image_plot,reform(den0(*,*,densize[3]/2,fix(i))),$
          xv,yv,/legend, /aspect, $
          title='den0 at t='+string(i*nout*dt)
        
        vsqr=reform(vsqrx0(*,*,densize[3]/2,fix(i))+ $
                    vsqry0(*,*,densize[3]/2,fix(i)))
        vsqr=vsqr/mean(vsqr)-1
        image_plot,vsqr,$
          xv,yv,/legend,/aspect, $
          title='v0*v0 at t='+string(i*nout*dt)
        
        vsqr2=reform(vsqrx1(*,*,densize[3]/2,fix(i))+ $
                     vsqry1(*,*,densize[3]/2,fix(i)))
        vsqr2=vsqr2/mean(vsqr2)-1
        image_plot,vsqr2,$
          xv,yv,/legend,/aspect, $ 
          title='v1*v1 at t='+string(i*nout*dt)
        
        date_plot,title
    endfor

endif else begin
; In 1D
    
    !p.multi=[0,1,2]
    if n_elements(den_out_subcycle0) eq 0 then begin
        if n_elements(den_out_subcycle) gt 0 then $
          den_out_subcycle0 = den_out_subcycle else $
          den_out_subcycle0 = 1
    endif

    image_plot,reform(den0),$
          xv,nout*dt*den_out_subcycle0*findgen(densize[4]),/legend,  $
          title='den0',xtitle='den0',ytitle='time'

    if n_elements(vsqr_out_subcycle0) eq 0 then begin
        if n_elements(vsqr_out_subcycle) gt 0 then $
          vsqr_out_subcycle0 = vsqr_out_subcycle else $
          vsqr_out_subcycle0 = 1
    endif
    image_plot,reform(vsqrx0),$
          xv,nout*dt*vsqr_out_subcycle0*findgen((size(vsqrx0))[4]),/legend, $
          title='vxsqr',xtitle='vxsqr',ytitle='time'
        
        date_plot,title
    
endelse

end




