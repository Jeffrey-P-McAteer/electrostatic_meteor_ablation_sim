; A script to produce vtk files from den1 for meteor runs
@eppic.i
if n_elements(color_on) eq 0 then color_on=1
if n_elements(iskip) eq 0 then iskip=1
iskip_bck=iskip

@params_in.pro
nx3=nx2/nsubdomains

if (n_elements(istart) eq 0) then istart=0
if (n_elements(iend) eq 0) then iend=-1

sizepertime=long64((nx3*1.0*ny2))

 den1=read_domains('den1.bin',[nx3,ny2,nz2], ndomains=nsubdomains, $
                   order=[2,1,0],/binary,skip=iskip,$
                   first=long64(istart*sizepertime),last=long64(iend*sizepertime))
; Fix the nomalization                                                                         
den1=(den1+1)*n0d1*param6_1

;
den1 = filter3d(den1, fwidth=2, fsteep=3)

vtk_output, den1,'den1Filtered2VTK',[dx*nout_avg,dy*nout_avg,dz*nout_avg,dt*nout]


