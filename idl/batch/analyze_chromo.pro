; Analyze the output of an eppic run for E-region waves, 
; making a set of standardized outputs

@eppic.i
if n_elements(color_on) eq 0 then color_on=1
if n_elements(iskip) eq 0 then iskip=1

; Velocity distribution and Temperature calculation
iskip=1
!p.multi=[0,2,2]
ps,'temp-vel3d-bw.ps'
@moment_plot_ep.pro
stop_ps

!p.multi=0

;Plot out mean and max of E^2 averaged over the simulation box:
@phiep
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
decades=2
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

;if color_on eq 0 then ps,'vdist0-bw.ps',/landscape $
;else ps,'vdist0-c.ps',/landscadispe,/color
;!p.multi=[0,2,2] 
;nvdist0=(size(vdist0))[4] 
;for i=0.,nvdist0-1.,(nvdist0-1)/3. do $
;  image_plot,reform(vdist0(*,*,(size(vdist0))[3]/2,fix(i))),$
;  vdistvx0,vdistvy0,/zlog,decades=3,/legend,$
;  title='Vel. Dist. 0 at t='+string(long(i)*nout*dt),xtitle='Vx',ytitle='Vy' & $
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
;  title='Vel. Dist. 1 at t='+string(long(i)*nout*dt),xtitle='Vx',ytitle='Vy' & $
;  date_plot,title
;stop_ps

;;Drift calculations
;!p.multi=[0,2,0]
;ps,'drifts.ps'
;.r drift3d.pro
;stop_ps

iskip=1
@eppic.i
@params_in.pro
nx3=nx2/nsubdomains

if (n_elements(istart) eq 0) then istart=0
if (n_elements(iend) eq 0) then iend=-1

sizepertime=long64((nx3*1.0*ny2))
if(ndim_space eq 3) then sizepertime = sizepertime*nz2

den0=read_domains('den0.bin',[nx3,ny2,nz2], ndomains=nsubdomains, $
                  order=[2,1,0],/binary,skip=iskip,$
                  first=long64(istart*sizepertime),last=long64(iend*sizepertime))
ps,'den0_image-c.ps',/landscape,/color
plot_set,den0,xv,yv,title='Density 0 at t=',dt_scale=long(nout)*iskip*dt, label=title
stop_ps

if (ndim_space eq 2) then image_mpeg, den0,/selfnormalize,filename='denxy.mpg' 

if (ndim_space eq 3) then image_mpeg, den0(*,*,(size(den0))[3]/2,*),/selfnormalize,filename='denxy.mpg' 
if (ndim_space eq 3) then image_mpeg, den0(*,(size(den0))[2]/2,*,*),/selfnormalize,filename='denxz.mpg'
if (ndim_space eq 3) then image_mpeg, den0((size(den0))[1]/2,*,*,*),/selfnormalize,filename='denyz.mpg'

den0=0

den1=read_domains('den1.bin',[nx3,ny2,nz2], ndomains=nsubdomains, $
                  order=[2,1,0],/binary,skip=iskip,$
                  first=long64(istart*sizepertime),last=long64(iend*sizepertime))
ps,'den1_image-c.ps',/landscape,/color
plot_set,den1,xv,yv,title='Density 1 at t=',dt_scale=long(nout)*iskip*dt, label=title
stop_ps
den1=0

if ndist gt 2 then $
den2=read_domains('den2.bin',[nx3,ny2,nz2], ndomains=nsubdomains, $
                  order=[2,1,0],/binary,skip=iskip,$
                  first=long64(istart*sizepertime),last=long64(iend*sizepertime)) & $
   ps,'den2_image-c.ps',/landscape,/color & $
   plot_set,den2,xv,yv,title='Density 2 at t=',dt_scale=long(nout)*iskip*dt, label=title & $
   stop_ps
den2=0

if ndist gt 3 then $
den3=read_domains('den3.bin',[nx3,ny2,nz2], ndomains=nsubdomains, $
                  order=[2,1,0],/binary,skip=iskip,$
                  first=long64(istart*sizepertime),last=long64(iend*sizepertime)) & $
   ps,'den3_image-c.ps',/landscape,/color & $
   plot_set,den3,xv,yv,title='Density 3 at t=',dt_scale=long(nout)*iskip*dt, label=title & $
   stop_ps
den3=0

if ndist gt 4 then $
den4=read_domains('den4.bin',[nx3,ny2,nz2], ndomains=nsubdomains, $
                  order=[2,1,0],/binary,skip=iskip,$
                  first=long64(istart*sizepertime),last=long64(iend*sizepertime)) & $
   ps,'den4_image-c.ps',/landscape,/color & $
   plot_set,den4,xv,yv,title='Density 4 at t=',dt_scale=long(nout)*iskip*dt, label=title & $
   stop_ps
den4=0


nden0=size(den0)

!p.multi=0
ps,'denavg.ps'
.r denavg3d
stop_ps

densize=size(den0)
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
if n_elements(Bx) ne 0 then if (Bx ne 0) then begin den0=transpose(den0,[2,1,0,3]) & densize=size(den0)

;To save memory eliminate den0
den0=0

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
