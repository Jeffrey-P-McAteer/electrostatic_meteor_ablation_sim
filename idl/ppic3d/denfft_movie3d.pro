densize=size(den1)
nt=densize(4)

if n_elements(den1f) lt 2 then begin
    den1f=den1
    for i=0,nt-1 do den1f(*,*,*,i)=shift(reform(abs(fft(den1(*,*,*,i))),nx2,ny2,nz2),nx2/2,ny2/2,nz2/2)
    for i=0,nt-1 do den1f(nx2/2,ny2/2,nz2/2,i)=0
endif

;Find the range over which there is interesting content:
den1f_sum=total(den1f,4)
nzm=densize[3]/2
range=image_range(den1f_sum[*,*,nzm] gt max(den1f_sum[*,*,0])/100.)

ixmin=nx2/2
ixmax=range[0,1]
iymin=range[1,0]
iymax=range[1,1]

image_movie,den1f(ixmin:ixmax,iymin:iymax,nzm,*),expand=4,/zlog,decades=2
end
