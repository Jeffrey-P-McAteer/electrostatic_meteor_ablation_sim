;Convert the electrostatic potential, phi, into electric field data, E^2

@phi3d

phisize=size(phi)
nt=phisize(phisize(0))
E2avg=fltarr(nt)
;Calculate E2 in real space:
E2=phi
for i=0,nt-1 do begin
    E=gradient(phi(*,*,*,i),dx=dx*nout_avg,dy=dy*nout_avg,dz=dz*nout_avg)
    if (nz2 gt 1) then begin
        E2(*,*,*,i)=E(*,*,*,0)^2+E(*,*,*,1)^2+E(*,*,*,2)^2 
    endif else begin
	if (ny2 gt 1) then E2(*,*,0,i)=E(*,*,0)^2+E(*,*,1)^2 $
        else E2(*,0,0,i)=E^2
    endelse
    E2avg(i)=mean(E2(*,*,*,i)) 
endfor

end
