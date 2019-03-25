; This function calculates the envelope of E assuming phi, divj and wp
; are known


@phi3d

nx2=long(nx/nout_avg) > 1
ny2=long(ny/nout_avg) > 1
nz2=long(nz/nout_avg) > 1

if (n_elements(divj) ne n_elements(phi)) then begin
    if (file_exist('divj.bin')) then $
      divj=readarray('divj.bin',[nx2,ny2,nz2],order=[2,1,0],/binary,skip=iskip)
endif

II=complex(0.,1.)
wp=1

; Calculate A
divAsize=size(divj)
A2=phi
time=fltarr(divAsize(4))
for i=0,divAsize(4)-1 do begin
    divE=-laplace(phi(*,*,*,i),dx=dx*nout_avg,dy=dy*nout_avg,dz=dz*nout_avg)
    t=i*dt*nout*iskip
    time(i)=t
    divA=(divE + II/wp/eps*divj(*,*,*,i))*exp(-II*wp*t)
    poten=poisson_invert(divA, dx=dx*nout_avg, dy=dy*nout_avg, dz=dz*nout_avg)
    A=gradient(poten,dx=dx*nout_avg, dy=dy*nout_avg, dz=dz*nout_avg)
    if (size(A))(0) gt 1 then $
      A2(*,*,*,i)=total(float(A*conj(A)),(size(A))(0)) else $
      A2(*,*,*,i)=float(A*conj(A))
endfor

end
