; Analyze the output of an eppic run for 150km echoes
; making a set of standardized outputs

@eppic.i
if n_elements(color_on) eq 0 then color_on=1
if n_elements(iskip) eq 0 then iskip=1
if (nx ge 2*ny) then nplot=8 else nplot=16

; Velocity distribution and Temperature calculation
iskip=1
!p.multi=[0,2,2]
@params_in.pro
nx3=nx2/nsubdomains

if (n_elements(istart) eq 0) then istart=0
if (n_elements(iend) eq 0) then iend=-1

sizepertime=long64((nx3*1.0*ny2))
if(ndim_space eq 3) then sizepertime = sizepertime*nz2

den0=read_domains('den0.bin',[nx3,ny2,nz2], ndomains=nsubdomains, $
                  order=[2,1,0],/binary,skip=iskip,$
                  first=long64(istart*sizepertime),last=long64(iend*sizepertime))
densize=size(den0)
!p.multi=[0,2,2]
ps,'den0wk.ps',/color,/landscape

denkw_slice_images,den0,[0,90,5,85],[nx*nsubdomains*dx,ny*dy],dt*nout*iskip*(densize[4]-1)
stop_ps

den0=0

; Make ion plots
den1=read_domains('den1.bin',[nx3,ny2,nz2], ndomains=nsubdomains, $
                  order=[2,1,0],/binary,skip=iskip,$
                  first=long64(istart*sizepertime),last=long64(iend*sizepertime))
den1size=size(den1)

!p.multi=[0,2,2]

ps,'den1wk.ps',/color,/landscape
denkw_slice_images,den1,[0,90,5,85],[nx*nsubdomains*dx,ny*dy],dt*nout*iskip*(den1size[4]-1)*den_out_subcycle1
stop_ps

print, 'Ended normally'
