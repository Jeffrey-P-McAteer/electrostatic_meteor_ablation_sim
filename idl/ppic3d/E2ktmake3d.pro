; Make E2(k,t) from phi

@phi3d

Esize=size(phi)

if n_elements(ifirst) eq 0 then ifirst=0
if n_elements(ilast) eq 0 then ilast=Esize(4)-1
if n_elements(iskip) eq 0 then iskip=1

;Set up the k vector for the axes
;and construct k^2 array
kxv=shift(findgen(Esize(1))*2*!PI/(nx*dx),nx2/2)
kxv(0:nx2/2-1)=kxv(0:nx2/2-1)-2*!PI/(dx*nout_avg)
k2=kxv^2
if Esize(0) gt 2 then begin
    kyv=shift(findgen(Esize(2))*2*!PI/(ny*dy),ny2/2)
    kyv(0:ny2/2-1)=kyv(0:ny2/2-1)-2*!PI/(dy*nout_avg)
    k2=((fltarr(ny2)+1.)##kxv^2+kyv^2##(fltarr(nx2)+1.))
endif
if Esize(0) gt 3 and Esize(3) gt 1 then begin
    kzv=shift(findgen(Esize(3))*2*!PI/(nz*dz),nz2/2)
    kzv(0:nz2/2-1)=kzv(0:nz2/2-1)-2*!PI/(dz*nout_avg)
    k2store=k2
    k2=fltarr(nx2,ny2,nz2)
    for iz=0, nz2-1 do k2(*,*,iz)= k2store + kzv(iz)^2
endif


E2kt=phi
for i=ifirst,ilast,iskip do begin
    if Esize(0) gt 3 and Esize(3) gt 1 then begin
        phifft=shift(fft(phi(*,*,*,i)),nx2/2,ny2/2,nz2/2) 
        E2kt(*,*,*,i)=float(phifft*conj(phifft))*k2
    endif else begin
        phifft=shift(fft(phi(*,*,*,i)),nx2/2,ny2/2) 
        E2kt(*,*,0,i)=float(phifft*conj(phifft))*k2
    endelse

endfor

end
