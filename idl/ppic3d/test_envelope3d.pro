; This routine tests envelope3d
@params_in.pro

snt=phisize(4)
phi=fltarr(nx2,ny2,nz2,nt)
divj=phi
fill=fltarr(ny2,nz2)+1.0
wp=1
kp=2*!PI/(dx*nx/4.)
dt2=dt*nout*iskip
for it=0, nt-1 do begin
    for ix=0, nx2-1 do begin
        phi(ix,*,*,it)=sin((kp*ix)*dx*nout_avg-wp*it*dt2)*fill
        divj(ix,*,*,it)=-kp^2*wp*cos((kp*ix)*dx*nout_avg-wp*it*dt2)*fill
    endfor
endfor
window,0
image_plot,reform(phi(*,*,0,nt-1)),/legend
window,2
image_plot,reform(divj(*,*,0,nt-1)),/legend

end
