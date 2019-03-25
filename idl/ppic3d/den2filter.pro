; Read den2 and filter out all modes not represented by at least 8 grid points
@ppic3d.i
nx2=long(nx/nout_avg) > 1
if (n_elements(ny) eq 0) then ny=1 
ny2=long(ny/nout_avg) > 1 
if (n_elements(nz) eq 0) then nz=1 
nz2=long(nz/nout_avg) > 1
nsize=nx2*ny2*nz2

if (ndist gt 2 and n_elements(den2) lt 2) then begin
    den2=readarray('den2.bin',[nx2,ny2,nz2],order=[2,1,0],/binary,skip=iskip)
endif

dsize=size(den2)

;xfilter=[findgen(dsize(1)/2),dsize(1)/2-findgen(dsize(1)/2)]/(dsize(1)/2)
;xfilter=1.0/(1+(xfilter/(dsize(1)/2-16))^10)

gap=6

xfilter=replicate(1,dsize(1))
xfilter[dsize(1)/2-gap:dsize(1)/2+gap-1]=0

yfilter=replicate(1,dsize(2))
yfilter[dsize(2)/2-gap:dsize(2)/2+gap-1]=0

zfilter=replicate(1,dsize(3))
zfilter[dsize(3)/2-gap:dsize(3)/2+gap-1]=0

filter=reform(den2(*,*,*,0))

for ix=0, dsize(1)-1 do $
  for iy=0, dsize(2)-1 do $
  for iz=0, dsize(3)-1 do $
  filter(ix,iy,iz)=xfilter(ix)*yfilter(iy)*zfilter(iz)

den2f=den2

for it=0,dsize(4)-1 do begin
    den2f(*,*,*,it)=abs(fft(fft(den2(*,*,*,it),1)*filter,-1))
endfor

end
