;Convert the electrostatic potential, phi, into electric field data

if (n_elements(Ex) le 2 and n_elements(phi) le 2) then $
  message,' Please run E2 first'

;Calculate E2 in real space:
E2=phi
for i=0,nt-1 do begin
    E=gradient(phi(*,*,i),dx=dx*nout_avg,dy=dy*nout_avg)
    if ny gt 1 then E2(*,*,i)=E(*,*,0)^2+E(*,*,1)^2 $
    else E2(*,0,i)=E^2
endfor
delvar,E

;endif else if (E_in_k_space ne 'false') then begin
;
;Calculate E^2 in k space:
;    E2=fltarr(nx2,ny2,nt)
;    kxv=shift(findgen(Esize(1))*2*!PI/(nx*dx),nx2/2)
;    kxv(0:nx2/2-1)=kxv(0:nx2/2-1)-2*!PI/(dx*nout_avg)
;    if ny gt 1 then begin
;        kyv=shift(findgen(Esize(2))*2*!PI/(ny*dy),ny2/2)
;        kyv(0:ny2/2-1)=kyv(0:ny2/2-1)-2*!PI/(dy*nout_avg)
;        k2=((fltarr(ny2)+1.)##kxv^2+kyv^2##(fltarr(nx2)+1.))
;    endif else K2=kxv^2
;    
;    for iw=0,nt-1 do E2(*,*,iw)=(k2)*float(phi(*,*,iw)*conj(phi(*,*,iw)))
;
;endif 

end

