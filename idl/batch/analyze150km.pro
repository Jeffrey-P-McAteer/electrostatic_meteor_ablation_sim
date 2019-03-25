; Analyze the output of an eppic run for 150km echoes
; making a set of standardized outputs

@eppic.i
if n_elements(color_on) eq 0 then color_on=1
if n_elements(iskip) eq 0 then iskip=1
if (nx ge 2*ny) then nplot=8 else nplot=16

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

; Velocity distribution plots - dist 0+2
iskip=4
@params_in.pro
vdist_tmp=read_domains('vdist0.bin',[pnvx2,pnvy2,pnvz2], ndomains=nsubdomains, $
                       order=[2,1,0],/binary,skip=iskip)
; Each domain is now a piece of the matrix concatinated in the
; x-direction. Combine these:
vdist=vdist_tmp[0:pnvx0-1,*,*,*]
for i=1,nsubdomains-1 do vdist += vdist_tmp[i*pnvx0:(i+1)*pnvx0-1,*,*,*]

vdist_tmp=0
vdist_tmp=read_domains('vdist2.bin',[pnvx2,pnvy2,pnvz2], ndomains=nsubdomains, $
                       order=[2,1,0],/binary,skip=iskip)
; Each domain is now a piece of the matrix concatinated in the
; x-direction. Combine these:
for i=1,nsubdomains-1 do $
  vdist += vdist_tmp[i*pnvx2:(i+1)*pnvx2-1,*,*,*]
vdist_tmp=0

izsize=((size(vdist))(3))
vxs=findgen(pnvx0)/float(pnvx0)*(pvxmax0-pvxmin0) + pvxmin0
if n_elements(pnvy0) gt 0 then vys=findgen(pnvy0)/float(pnvy0)*(pvymax0-pvymin0) + pvymin0
if izsize gt 1 then vzs=findgen(pnvz0)/float(pnvz0)*(pvzmax0-pvzmin0) + pvzmin0 else vzs=0

ps,'vdist_e.ps',/landscape,/color
!p.multi=[0,2,2] 
nvdist=(size(vdist))[4] 
nimages=8.
for i=0.,nvdist-1.,(nvdist-1)/(nimages-1+1.0E-5) do $
  image_plot,reform(vdist(*,*,(size(vdist))[3]/2,fix(i))),$
  vxs,vys,/zlog,nlabels=10,/legend,$
  title='Vel. Dist. 0 at t='+string(long(i)*nout*dt*iskip*vdist_out_subcycle0),xtitle='Vx',ytitle='Vy' & $
  date_plot,title
stop_ps

;Make a movie of the evolving distributions
image_mpeg,vdist[*,*,(size(vdist))[3]/2,*],/zlog,decades=10,filename='vdist_e.mpg',quality=75

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

delvar, vdist

@eppic.i
@params_in.pro
iskip=1
nx3=nx2/nsubdomains

iend=169984

if (n_elements(istart) eq 0) then istart=0
if (n_elements(iend) eq 0) then iend=-1

sizepertime=long64((nx3*1.0*ny2))
if(ndim_space eq 3) then sizepertime = sizepertime*nz2

;print,"start should be",long64(istart*sizepertime)
;print,"end should be  ",long64(iend*sizepertime)

den0=read_domains('den0.bin',[nx3,ny2,nz2], ndomains=nsubdomains, $
                  order=[2,1,0],/binary,skip=iskip,$
                  first=long64(istart*sizepertime),last=long64(iend*sizepertime))

help,den0

;Plot a sequence (or all) of the density values:
if (nx ge 2*ny) then nplot=8 else nplot=16
if (nx ge 2*ny) then !p.multi=[0,2,4] else !p.multi=[0,4,4] 

nden0=size(den0)
densize=size(den0)

!p.multi=0
ps,'denavg.ps'
.r denavg3d
stop_ps

if (nx ge 2*ny) then !p.multi=[0,2,4] else !p.multi=[0,4,4] 
nden0=size(den0)
if color_on eq 0 then ps,'den0_image-bw.ps',/landscape $
else ps,'den0_image-c.ps',/landscape,/color
;loadct,1

for i=0.,nden0(4)-1.,(nden0(4)-1.)/(nplot-1.)*0.999999 do $
  image_plot, reform(den0(*,*,(size(den0))[3]/2,fix(i))),$
  xv,yv,/legend,/aspect,$
  title='Density 1 at t='+string(fix(i)*long(nout)*iskip*dt),xtitle='x',ytitle='y' & $
  date_plot,title
if (nden0[3] gt 1) then $
for i=0.,nden0(4)-1.,(nden0(4)-1.)/(nplot-1.)*0.999999 do $
  image_plot, reform(den0(*,(size(den0))[2]/2,*,fix(i))),$
  xv,zv,/legend,/aspect,$
  title='Density 1 at t='+string(fix(i)*long(nout)*iskip*dt),xtitle='x',ytitle='z' & $
  date_plot,title
if (nden0[3] gt 1) then $
for i=0.,nden0(4)-1.,(nden0(4)-1.)/(nplot-1.)*0.999999 do $
  image_plot, transpose(reform(den0((size(den0))[1]/2,*,*,fix(i)))),$
  zv,yv,/legend,/aspect,$
  title='Density 1 at t='+string(fix(i)*long(nout)*iskip*dt),xtitle='z',ytitle='y' & $
  date_plot,title
stop_ps

;Make movie of den0 (warning: very slow)
;image_mpeg,den0(*,*,0,*),filename='den0xy.mpg',quality=75

!p.multi=[0,2,2]
ps,'den0wk.ps',/color,/landscape

denkw_slice_images,den0,[0,5,45,85,88,89,89.5,90],[nx*nsubdomains*dx,ny*dy],dt*nout*iskip*(densize[4]-1),nlables=5,filename='den0'
date_plot,title 
denkw_slice_images,den0,[0,5,45,85,88,89,89.5,90],[nx*nsubdomains*dx,ny*dy],dt*nout*iskip*(densize[4]-1),zrange=[1E-3,1E-7],nlabels=5
date_plot,title 
stop_ps

;plot den0(w,|k|)
!p.multi=0
ps,'denkw_image.ps',/color
.r den_k_w_images.pro
stop_ps

; reformat den0 to the correct size
den0=reform(den0,[densize[1],densize[2],den0size[3],densize[4]])

normalize_spectra=0
if color_on eq 0 then ps,'fixed_k_spectra_bw.ps',/landscape $
else ps,'fixed_k_spectra.ps',/landscape,/color
!p.multi=[0,2,2] 
wavelength=[20, 30, 60, 120]*dx
if nz gt 1 then aspect_ang=[90,89,87,85,75]*!PI/180.
.r fixed_k_spectral_image

;Add a normalized figure.
normalize_spectra=1
.r fixed_k_spectral_image
date_plot, title+' renormalized hanning'
stop_ps

;;Plot the spectral change with k:
;!p.multi=0
;ps,'k_spectra_bw.ps'
;.r params_k_images
;stop_ps

den0=0

; Make ion plots
den1=read_domains('den1.bin',[nx3,ny2,nz2], ndomains=nsubdomains, $
                  order=[2,1,0],/binary,skip=iskip,$
                  first=long64(istart*sizepertime),last=long64(iend*sizepertime))

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

;Make movie of den1
;image_mpeg,den1(*,*,0,*),filename='den1xy.mpg',quality=75

densize=size(den1)
;plot denk(w,|k|)
!p.multi=[0,2,2]
iavg=(densize[4]+3)/4
ps,'denkt_avg.ps',/color,/landscape
.r denkt_avg.pro
stop_ps

;Make den1(w,k) figures for 2D system
den1size=size(den1)
nhalf=den1size[4]/2
;fden1=abs(shift(fft(reform(den1[*,*,0,nhalf:nhalf*2-1])),den1size[1]/2,den1size[2]/2,nhalf/2))
!p.multi=[0,2,2]

ps,'den1wk.ps',/color,/landscape
denkw_slice_images,den1,[0,5,45,85,88,89,89.5,90],[nx*nsubdomains*dx,ny*dy],dt*nout*iskip*(den1size[4]-1)*den_out_subcycle1,filename='den1'
date_plot,title 
denkw_slice_images,den1,[0,5,45,85,88,89,89.5,90],[nx*nsubdomains*dx,ny*dy],dt*nout*iskip*(den1size[4]-1)*den_out_subcycle1,zrange=[1E-6,1E-4]*0.2,wrange=[-.03,.03],krange=[0,.5]
date_plot,title 
stop_ps

; Make some figures of density spectra at fixed k for a range of angles:
; Use only the second half of the simulation:
nden=densize[densize[0]]
den1=den1[*,*,*,nden/2:nden-1]
densize=size(den1)
nden=densize[densize[0]]
if n_elements(Bx) ne 0 then if (Bx ne 0) then begin den1=transpose(den1,[2,1,0,3]) & densize=size(den1)

;To save memory eliminate den1
den1=0

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

print, 'Analyze150km ended normally'

cd,'../PES2D15VB.577CollBoost2xV160310'
; Do 2 at once
@eppic.i
if n_elements(color_on) eq 0 then color_on=1
if n_elements(iskip) eq 0 then iskip=1
if (nx ge 2*ny) then nplot=8 else nplot=16

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

; Velocity distribution plots - dist 0+2
iskip=4
@params_in.pro
vdist_tmp=read_domains('vdist0.bin',[pnvx2,pnvy2,pnvz2], ndomains=nsubdomains, $
                       order=[2,1,0],/binary,skip=iskip)
; Each domain is now a piece of the matrix concatinated in the
; x-direction. Combine these:
vdist=vdist_tmp[0:pnvx0-1,*,*,*]
for i=1,nsubdomains-1 do vdist += vdist_tmp[i*pnvx0:(i+1)*pnvx0-1,*,*,*]

vdist_tmp=0
vdist_tmp=read_domains('vdist2.bin',[pnvx2,pnvy2,pnvz2], ndomains=nsubdomains, $
                       order=[2,1,0],/binary,skip=iskip)
; Each domain is now a piece of the matrix concatinated in the
; x-direction. Combine these:
for i=1,nsubdomains-1 do $
  vdist += vdist_tmp[i*pnvx2:(i+1)*pnvx2-1,*,*,*]
vdist_tmp=0

izsize=((size(vdist))(3))
vxs=findgen(pnvx0)/float(pnvx0)*(pvxmax0-pvxmin0) + pvxmin0
if n_elements(pnvy0) gt 0 then vys=findgen(pnvy0)/float(pnvy0)*(pvymax0-pvymin0) + pvymin0
if izsize gt 1 then vzs=findgen(pnvz0)/float(pnvz0)*(pvzmax0-pvzmin0) + pvzmin0 else vzs=0

ps,'vdist_e.ps',/landscape,/color
!p.multi=[0,2,2] 
nvdist=(size(vdist))[4] 
nimages=8.
for i=0.,nvdist-1.,(nvdist-1)/(nimages-1+1.0E-5) do $
  image_plot,reform(vdist(*,*,(size(vdist))[3]/2,fix(i))),$
  vxs,vys,/zlog,nlabels=10,/legend,$
  title='Vel. Dist. 0 at t='+string(long(i)*nout*dt*iskip*vdist_out_subcycle0),xtitle='Vx',ytitle='Vy' & $
  date_plot,title
stop_ps

;Make a movie of the evolving distributions
image_mpeg,vdist[*,*,(size(vdist))[3]/2,*],/zlog,decades=10,filename='vdist_e.mpg',quality=75

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

delvar, vdist

@eppic.i
@params_in.pro
iskip=1
nx3=nx2/nsubdomains

iend=169984

if (n_elements(istart) eq 0) then istart=0
if (n_elements(iend) eq 0) then iend=-1

sizepertime=long64((nx3*1.0*ny2))
if(ndim_space eq 3) then sizepertime = sizepertime*nz2

;print,"start should be",long64(istart*sizepertime)
;print,"end should be  ",long64(iend*sizepertime)

den0=read_domains('den0.bin',[nx3,ny2,nz2], ndomains=nsubdomains, $
                  order=[2,1,0],/binary,skip=iskip,$
                  first=long64(istart*sizepertime),last=long64(iend*sizepertime))

help,den0

;Plot a sequence (or all) of the density values:
if (nx ge 2*ny) then nplot=8 else nplot=16
if (nx ge 2*ny) then !p.multi=[0,2,4] else !p.multi=[0,4,4] 

nden0=size(den0)
densize=size(den0)

!p.multi=0
ps,'denavg.ps'
.r denavg3d
stop_ps

if (nx ge 2*ny) then !p.multi=[0,2,4] else !p.multi=[0,4,4] 
nden0=size(den0)
if color_on eq 0 then ps,'den0_image-bw.ps',/landscape $
else ps,'den0_image-c.ps',/landscape,/color
;loadct,1

for i=0.,nden0(4)-1.,(nden0(4)-1.)/(nplot-1.)*0.999999 do $
  image_plot, reform(den0(*,*,(size(den0))[3]/2,fix(i))),$
  xv,yv,/legend,/aspect,$
  title='Density 1 at t='+string(fix(i)*long(nout)*iskip*dt),xtitle='x',ytitle='y' & $
  date_plot,title
if (nden0[3] gt 1) then $
for i=0.,nden0(4)-1.,(nden0(4)-1.)/(nplot-1.)*0.999999 do $
  image_plot, reform(den0(*,(size(den0))[2]/2,*,fix(i))),$
  xv,zv,/legend,/aspect,$
  title='Density 1 at t='+string(fix(i)*long(nout)*iskip*dt),xtitle='x',ytitle='z' & $
  date_plot,title
if (nden0[3] gt 1) then $
for i=0.,nden0(4)-1.,(nden0(4)-1.)/(nplot-1.)*0.999999 do $
  image_plot, transpose(reform(den0((size(den0))[1]/2,*,*,fix(i)))),$
  zv,yv,/legend,/aspect,$
  title='Density 1 at t='+string(fix(i)*long(nout)*iskip*dt),xtitle='z',ytitle='y' & $
  date_plot,title
stop_ps

;Make movie of den0 (warning: very slow)
;image_mpeg,den0(*,*,0,*),filename='den0xy.mpg',quality=75

!p.multi=[0,2,2]
ps,'den0wk.ps',/color,/landscape

denkw_slice_images,den0,[0,5,45,85,88,89,89.5,90],[nx*nsubdomains*dx,ny*dy],dt*nout*iskip*(densize[4]-1),nlables=5,filename='den0'
date_plot,title 
denkw_slice_images,den0,[0,5,45,85,88,89,89.5,90],[nx*nsubdomains*dx,ny*dy],dt*nout*iskip*(densize[4]-1),zrange=[1E-3,1E-7],nlabels=5
date_plot,title 
stop_ps

;plot den0(w,|k|)
!p.multi=0
ps,'denkw_image.ps',/color
.r den_k_w_images.pro
stop_ps

; reformat den0 to the correct size
den0=reform(den0,[densize[1],densize[2],den0size[3],densize[4]])

normalize_spectra=0
if color_on eq 0 then ps,'fixed_k_spectra_bw.ps',/landscape $
else ps,'fixed_k_spectra.ps',/landscape,/color
!p.multi=[0,2,2] 
wavelength=[20, 30, 60, 120]*dx
if nz gt 1 then aspect_ang=[90,89,87,85,75]*!PI/180.
.r fixed_k_spectral_image

;Add a normalized figure.
normalize_spectra=1
.r fixed_k_spectral_image
date_plot, title+' renormalized hanning'
stop_ps

;;Plot the spectral change with k:
;!p.multi=0
;ps,'k_spectra_bw.ps'
;.r params_k_images
;stop_ps

den0=0

; Make ion plots
den1=read_domains('den1.bin',[nx3,ny2,nz2], ndomains=nsubdomains, $
                  order=[2,1,0],/binary,skip=iskip,$
                  first=long64(istart*sizepertime),last=long64(iend*sizepertime))

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

;Make movie of den1
;image_mpeg,den1(*,*,0,*),filename='den1xy.mpg',quality=75

densize=size(den1)
;plot denk(w,|k|)
!p.multi=[0,2,2]
iavg=(densize[4]+3)/4
ps,'denkt_avg.ps',/color,/landscape
.r denkt_avg.pro
stop_ps

;Make den1(w,k) figures for 2D system
den1size=size(den1)
nhalf=den1size[4]/2
;fden1=abs(shift(fft(reform(den1[*,*,0,nhalf:nhalf*2-1])),den1size[1]/2,den1size[2]/2,nhalf/2))
!p.multi=[0,2,2]

ps,'den1wk.ps',/color,/landscape
denkw_slice_images,den1,[0,5,45,85,88,89,89.5,90],[nx*nsubdomains*dx,ny*dy],dt*nout*iskip*(den1size[4]-1)*den_out_subcycle1,filename='den1'
date_plot,title 
denkw_slice_images,den1,[0,5,45,85,88,89,89.5,90],[nx*nsubdomains*dx,ny*dy],dt*nout*iskip*(den1size[4]-1)*den_out_subcycle1,zrange=[1E-6,1E-4]*0.2,wrange=[-.03,.03],krange=[0,.5]
date_plot,title 
stop_ps

; Make some figures of density spectra at fixed k for a range of angles:
; Use only the second half of the simulation:
nden=densize[densize[0]]
den1=den1[*,*,*,nden/2:nden-1]
densize=size(den1)
nden=densize[densize[0]]
if n_elements(Bx) ne 0 then if (Bx ne 0) then begin den1=transpose(den1,[2,1,0,3]) & densize=size(den1)

;To save memory eliminate den1
den1=0

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

print, 'Analyze150km ended normally'
