; Analyze the output of an eppic run for meteors, 
; making a set of standardized outputs

;iskip=8

@eppic.i
if n_elements(color_on) eq 0 then color_on=1
if n_elements(iskip) eq 0 then iskip=1
iskip_bck=iskip

@params_in.pro
nx3=nx2/nsubdomains

if (n_elements(istart) eq 0) then istart=0
if (n_elements(iend) eq 0) then iend=-1

sizepertime=long64((nx3*1.0*ny2))
if(ndim_space eq 3) then sizepertime = sizepertime*nz2


; Velocity distribution and Temperature calculation
;iskip=1
!p.multi=[0,2,2]
ps,'temp-vel3d-bw.ps'
@moment_plot_ep.pro
stop_ps

!p.multi=0

 den1=read_domains('den1.bin',[nx3,ny2,nz2], ndomains=nsubdomains, $
                   order=[2,1,0],/binary,skip=iskip,$
                   first=long64(istart*sizepertime),last=long64(iend*sizepertime))
; Fix the nomalization
den1=(den1+1)*n0d1*param6_1

;Plot a sequence (or all) of the density values:
if (nx ge 2*ny) then nplot=8 else nplot=16
if (nx ge 2*ny) then !p.multi=[0,2,4] else !p.multi=[0,4,4] 

;nden0=size(den1)

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
nden1=size(den1)
if color_on eq 0 then ps,'den1_image-bw.ps',/landscape $
else ps,'den1_image-c.ps',/landscape,/color
;loadct,1

for i=0.,nden1(4)-1.,(nden1(4)-1.)/(nplot-1.)*0.999999 do $
  image_plot, reform(den1(*,*,(size(den1))[3]/2,fix(i))),$
  xv,yv,/legend,/aspect,$
  title='Density 1 at t='+string(fix(i)*long(nout)*iskip*dt),xtitle='x',ytitle='y' & $
  date_plot,title
if (nden1[3] gt 1) then $
for i=0.,nden1(4)-1.,(nden1(4)-1.)/(nplot-1.)*0.999999 do $
  image_plot, reform(den1(*,(size(den1))[2]/2,*,fix(i))),$
  xv,zv,/legend,/aspect,$
  title='Density 1 at t='+string(fix(i)*long(nout)*iskip*dt),xtitle='x',ytitle='z' & $
  date_plot,title
if (nden1[3] gt 1) then $
for i=0.,nden1(4)-1.,(nden1(4)-1.)/(nplot-1.)*0.999999 do $
  image_plot, transpose(reform(den1((size(den1))[1]/2,*,*,fix(i)))),$
  zv,yv,/legend,/aspect,$
  title='Density 1 at t='+string(fix(i)*long(nout)*iskip*dt),xtitle='z',ytitle='y' & $
  date_plot,title
stop_ps

; Make movies
image_mpeg, den1(*,*,(size(den1))[3]/2,*),/selfnormalize,filename='denxy.mpg'
image_mpeg, den1(*,(size(den1))[2]/2,*,*),/selfnormalize,filename='denxz.mpg'
image_mpeg, den1((size(den1))[1]/2,*,*,*),/selfnormalize,filename='denyz.mpg'

image_mpeg, den1(*,*,(size(den1))[3]/2,*),/zlog,filename='denxy_log.mpg'
image_mpeg, den1(*,(size(den1))[2]/2,*,*),/zlog,filename='denxz_log.mpg'
image_mpeg, den1((size(den1))[1]/2,*,*,*),/zlog,filename='denyz_log.mpg'

densize=size(den1)
;plot den0(w,|k|)
!p.multi=[0,2,2]
iavg=(densize[4]+3)/4
ps,'denkt_avg.ps',/color,/landscape
.r denkt_avg.pro
stop_ps

; Make some figures of density spectra at fixed k for a range of angles:
; Use only the second half of the simulation:
nden=densize[densize[0]]
den0=den1[*,*,*,nden/2:nden-1]
densize=size(den0)
nden=densize[densize[0]]
if n_elements(Bx) ne 0 then if (Bx ne 0) then begin den0=transpose(den0,[2,1,0,3]) & densize=size(den0)

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


@eppic.i

!p.multi=0

;Plot out mean and max of E^2 averaged over the simulation box:
iskip=iskip_bck
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

iskip=iskip_bck
;!p.multi=[0,3,3]
;avg_3d=1
;ps,'E2kt-avg-3x3-3d-bw.ps',/landscape
;.run E2kt3d.pro
;stop_ps
;delvar, phi

delvar, phi, E

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

print, 'Analyze_meteor ended normally'
