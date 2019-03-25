; Analyze the output of a ppic3d run, making a set of standardized
; outputs

if n_elements(color_on) eq 0 then color_on=1
if n_elements(iskip) eq 0 then iskip=1
if n_elements(Bz) eq 0 and n_elements(Bx) ne 0 then Bz=Bx

!p.multi=0

; Velocity distribution and Temperature calculation
iskip=1
!p.multi=[0,2,2]
ps,'temp-vel3d-bw.ps'
@moment_plot_ep.pro
stop_ps


;;Drift calculations
;!p.multi=[0,2,0]
;ps,'drifts.ps'
;.r drift3d.pro
;stop_ps

;Plot a sequence (or all) of the density values:
if (nx ge 2*ny) then nplot=8 else nplot=16
if (nx ge 2*ny) then !p.multi=[0,2,4] else !p.multi=[0,4,4] 
;@denep
;nden0=size(den0)
;if color_on eq 0 then ps,'den0_image-bw.ps',/landscape $
;else ps,'den0_image-c.ps',/landscape,/color
;for i=0.,(nden0(4)-1.),(nden0(4)-1.)/(nplot-1.) do $
;  image_plot, reform(den0(*,*,(size(den0))[3]/2,fix(i))),$
; xv,yv,/legend,$
;  title='Density 0 at t='+string(i*nout*dt),xtitle='x',ytitle='y' & $
;  date_plot,title
;stop_ps


if (nx ge 2*ny) then !p.multi=[0,2,4] else !p.multi=[0,4,4] 
nden0=size(den0)
if color_on eq 0 then ps,'den0_image-bw.ps',/landscape $
else ps,'den0_image-c.ps',/landscape,/color
for i=0.,nden0(4)-1.,(nden0(4)-1.)/(nplot-1.)*0.999999 do $
  image_plot, reform(den0(*,*,(size(den0))[3]/2,fix(i))),$
  xv,yv,/legend,/aspect,$
  title='Density 1 at t='+string(fix(i)*nout*dt),xtitle='x',ytitle='y' & $
  date_plot,title
if (nden0[3] gt 1) then $
for i=0.,nden0(4)-1.,(nden0(4)-1.)/(nplot-1.)*0.999999 do $
  image_plot, reform(den0(*,(size(den0))[2]/2,*,fix(i))),$
  xv,zv,/legend,/aspect,$
  title='Density 1 at t='+string(fix(i)*nout*dt),xtitle='x',ytitle='z' & $
  date_plot,title
if (nden0[3] gt 1) then $
for i=0.,nden0(4)-1.,(nden0(4)-1.)/(nplot-1.)*0.999999 do $
  image_plot, transpose(reform(den0((size(den0))[2]/2,*,*,fix(i)))),$
  zv,yv,/legend,/aspect,$
  title='Density 1 at t='+string(fix(i)*nout*dt),xtitle='z',ytitle='y' & $
  date_plot,title
stop_ps

!p.multi=0
ps,'denavg.ps'
.r denavg3d
stop_ps

;plot den0(w,|k|)
!p.multi=[0,2,2]
iavg=(densize[4]+3)/4
ps,'denkt_avg.ps',/color,/landscape
.r denkt_avg.pro
stop_ps

; Make some figures of density spectra at fixed k for a range of angles:
; Use only the second half of the simulation:
nden=densize[densize[0]]
den0=den0[*,*,*,nden/2:nden-1]
densize=size(den0)
nden=densize[densize[0]]
if (Bx ne 0) then den0=transpose(den0,[2,1,0,3]) & densize=size(den0)

;To save memory eliminate den1
den1=0

normalize_spectra=0
if color_on eq 0 then ps,'fixed_k_spectra_bw.ps',/landscape $
else ps,'fixed_k_spectra.ps',/landscape,/color
!p.multi=[0,2,2] 
wavelength=[0.75, 1.5, 3.0, 6.0]
if nz gt 1 then aspect_ang=[90,89,87,85,75]*!PI/180.
.r fixed_k_spectral_image

;Add a normalized figure.
normalize_spectra=1
.r fixed_k_spectral_image
date_plot, title+' renormalized hanning'
stop_ps

;Plot the spectral change with k:
!p.multi=0
ps,'k_spectra_bw.ps'
.r params_k_images
stop_ps

;plot den0(w,|k|)
!p.multi=0
ps,'denkw_image.ps',/color
.r den_k_w_images.pro
stop_ps



;Look at spacial distribution of pressures
;iskip=1
;if color_on eq 0 then ps,'pressure-bw.ps',/landscape $
;else ps,'pressure-c.ps',/landscape,/color
;.r pressure3d
;stop_ps

;Look at spacial distribution of temps 
;iskip=1
;if color_on eq 0 then ps,'vsqr3d-bw.ps',/landscape $
;else ps,'vsqr3d-c.ps',/landscape,/color
;.r vsqr3d
;stop_ps

print, 'AnalyzeFBep ended normally'
