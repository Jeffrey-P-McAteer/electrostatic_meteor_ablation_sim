;Convert the electrostatic potential, phi, into electric field data

@phi3d

; Set up E field arrays - size them
if (nx gt 1) then Ex=phi
if (ny gt 1) then Ey=phi
if (nz gt 1) then Ez=phi
;Calculate E2 in real space:
nt=(size(phi))(4)

for i=0,nt-1 do begin
    E=-gradient(phi[*,*,*,i],dx=dx*nout_avg,dy=dy*nout_avg,dz=dz*nout_avg)
    Esize=size(E)
    if (nz2 gt 1) then begin
        Ex(*,*,*,i)=E(*,*,*,0)
        Ey(*,*,*,i)=E(*,*,*,1)
        Ez(*,*,*,i)=E(*,*,*,2)
    endif else if (ny2 gt 1) then begin
        Ex(*,*,0,i)=E(*,*,0)
        Ey(*,*,0,i)=E(*,*,1)
    endif else Ex(*,0,0,i)=E
endfor



end
