; Make a set of denwk plots with the same amplitude range to compare.  
; Note: the size must be the same for all (eppic.i only read for the
; first directory)

;'PES2D15VB.5halfCollBoostV150925','PES2D15VB.5BoostV151007',
dirs=['PES2D15VB.5Boost2xmi256V151012','PES2D15VB.5Boost2xV151009','PES2D0VB.5B2xV151017','PES2D10VB.5Boost2xV151021','PES2D15VB0Boost2xV151023']

ndirs=n_elements(dirs)
iskip=1

@$SCRATCH/PES2D15VB.5Boost2xmi256V151012/eppic.i 
if n_elements(ndim_space) eq 0 then ndim_space=1+(ny gt 1)+(nz gt 1)
if ndim_space lt 3 then nz=1
if ndim_space lt 3 then dz=0.
if ndim_space lt 2 then ny=1
if ndim_space lt 2 then dy=0.
if n_elements(nsubdomains) eq 0 then nsubdomains=1 
nx2=long(nx*nsubdomains/nout_avg) > 1
ny2=long(ny/nout_avg) > 1
nz2=long(nz/nout_avg) > 1

nx3=nx2/nsubdomains
sizepertime=long64((nx3*1.0*ny2))

!p.multi=[0,2,2]
ps,'den1wk_compare.ps',/color,/landscape
; Make ion plots

for i=0,ndirs-1 do begin
   cd,dirs[i]
   print,dirs[i]
   den1=read_domains('den1.bin',[nx3,ny2,nz2], ndomains=nsubdomains, $
                     order=[2,1,0],/binary)
                     
   den1size=size(den1) 
   nhalf=den1size[4]/2 
   denkw_slice_images,den1,[0,90,5,85],[nx*nsubdomains*dx,ny*dy],dt*nout*iskip*(den1size[4]-1)*den_out_subcycle1,zrange=[1E-6,1E-4]*0.4
   date_plot,dirs[i]+': '+title 
   cd,'..' 
endfor

stop_ps

end
