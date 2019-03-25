; Analyze the output of a ppic3d run, making a set of standardized
; outputs

if n_elements(color_on) eq 0 then color_on=1
iskip=1
!p.multi=0

;Plot out mean and max of E^2 averaged over the simulation box:
ps,'E2avg.ps'
.run E2avg3d
stop_ps
;Plot a sequence (or all) of the electric potential (phi) values:
phisize=size(phi)

if (nx ge 2*ny) then !p.multi=[0,2,4] else !p.multi=[0,4,4]
if (nx ge 2*ny) then nplot=8 else nplot=16
iskip=round((phisize(4)-1)/(nplot))+1
if color_on eq 0 then ps,'phi_image-bw.ps',/landscape $
else ps,'phi_image-c.ps',/landscape,/color
.run phi_image3d
stop_ps

;Plot a sequence (or all) of the electric field energy (E^2) values:
if nx2 ge 2*ny2 then !p.multi=[0,2,4] else !p.multi=[0,4,4]
if ny2 eq 1 then !p.multi=0
if nx2 ge 2*ny2 then nplot=8 else nplot=16
iskip=round((phisize(4)-1)/(nplot))+1
if color_on eq 0 then ps,'E2_image-bw.ps',/landscape $
else ps,'E2_image-c.ps',/landscape,/color
.run E2image3d
stop_ps

iskip=round((phisize(4)-1)/27)+1
!p.multi=[0,3,3]
avg_3d=0
;Only plot the inner component of k space
xfrac=.3
yfrac=.3
zfrac=.3
if color_on eq 0 then ps,'E2kt-3x3-3d-bw.ps',/landscape $
else ps,'E2kt-3x3-3d-c.ps',/landscape,/color 
.run E2kt3d.pro
stop_ps

;!p.multi=[0,3,3]
;avg_3d=1
;ps,'E2kt-avg-3x3-3d-bw.ps',/landscape
;.run E2kt3d.pro
;stop_ps
;delvar, phi

delvar, phi, E

; Velocity distribution and Temperature calculation
iskip=1
!p.multi=[0,2,2]
ps,'temp-vel3d-bw.ps'
@moment_plots.pro
stop_ps

;if color_on eq 0 then ps,'vdist0-bw.ps',/landscape $
;else ps,'vdist0-c.ps',/landscadispe,/color
;!p.multi=[0,2,2] 
;nvdist0=(size(vdist0))[4] 
;for i=0.,nvdist0-1.,(nvdist0-1)/3. do $
;  image_plot,reform(vdist0(*,*,(size(vdist0))[3]/2,fix(i))),$
;  vdistvx0,vdistvy0,/zlog,decades=3,/legend,$
;  title='Vel. Dist. 0 at t='+string(i*nout*dt),xtitle='Vx',ytitle='Vy' & $
;  date_plot,title
;stop_ps
;
;if (n_elements(vdist1) ge 1) then $
;  if color_on eq 0 then ps,'vdist1-bw.ps',/landscape $
;  else ps,'vdist1-c.ps',/landscape,/color & $
;  nvdist1=(size(vdist1))[4]  & $
;  for i=0.,nvdist1-1.,(nvdist1-1.)/3. do $
;  image_plot,reform(vdist1(*,*,(size(vdist1))[3]/2,fix(i))),$
;  vdistvx1,vdistvy1,/zlog,decades=3,/legend, $
;  title='Vel. Dist. 1 at t='+string(i*nout*dt),xtitle='Vx',ytitle='Vy' & $
;  date_plot,title
;stop_ps

;;Drift calculations
;!p.multi=[0,2,0]
;ps,'drifts.ps'
;.r drift3d.pro
;stop_ps

;Plot a sequence (or all) of the density values:
if (nx ge 2*ny) then nplot=8 else nplot=16
if (nx ge 2*ny) then !p.multi=[0,2,4] else !p.multi=[0,4,4] 

@den3d
nden0=size(den0)
if color_on eq 0 then ps,'den0_image-bw.ps',/landscape $
else ps,'den0_image-c.ps',/landscape,/color
for i=0.,(nden0(4)-1.),(nden0(4)-1.)/(nplot-1.) do $
  image_plot, reform(den0(*,*,(size(den0))[3]/2,fix(i))),$
 xv,yv,/legend,$
  title='Density 0 at t='+string(i*nout*dt),xtitle='x',ytitle='y' & $
  date_plot,title
stop_ps

if (nx ge 2*ny) then !p.multi=[0,2,4] else p.multi=[0,4,4] 
nden1=size(den1)
if color_on eq 0 then ps,'den1_image-bw.ps',/landscape $
else ps,'den1_image-c.ps',/landscape,/color
for i=0.,nden1(4)-1.,(nden1(4)-1.)/(nplot-1.) do $
  image_plot, reform(den1(*,*,(size(den1))[3]/2,fix(i))),$
  xv,yv,/legend,$
  title='Density 1 at t='+string(i*nout*dt),xtitle='x',ytitle='y' & $
  date_plot,title
stop_ps

; Make some figures of density spectra at fixed k for a range of angles:
; Use only the second half of the simulation:

nden=(size(den1))[4]
den0=den1[*,*,*,nden/2:nden-1]
densize=size(den0)
nden=densize[densize[0]]

normalize_spectra=0
if color_on eq 0 then ps,'fixed_k_spectra_bw.ps',/landscape $
else ps,'fixed_k_spectra.ps',/landscape,/color
!p.multi=[0,2,2] 
wavelength=[0.75, 1.5, 3.0, 6.0]
.r fixed_k_spectral_image

;Add a normalized figure.
normalize_spectra=1
.r fixed_k_spectral_image
date_plot, title+' renormalized hanning'
stop_ps

;Plot the spectral change with k:
ps,'k_spectra_bw.ps'
.r params_k_images
stop_ps

;plot den0(w,|k|)
!p.multi=0
ps,'denkw_image.ps',/color
.r den_k_w_images.pro
stop_ps

;Look at spacial distribution of pressures
iskip=1
if color_on eq 0 then ps,'pressure-bw.ps',/landscape $
else ps,'pressure-c.ps',/landscape,/color
.r pressure3d
stop_ps

;Look at spacial distribution of temps 
iskip=1
if color_on eq 0 then ps,'vsqr3d-bw.ps',/landscape $
else ps,'vsqr3d-c.ps',/landscape,/color
.r vsqr3d
stop_ps


