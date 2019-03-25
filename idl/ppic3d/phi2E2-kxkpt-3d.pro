; This routine creates E2 in k space

@phi3d

if n_elements(ifirst) eq 0 then ifirst=0
if n_elements(ilast) eq 0 then ilast=phisize(4)-1
if n_elements(iskip) eq 0 then iskip=1
if n_elements(avg_3d) eq 0 then avg_3d=0

;Set up the k vector for the axes
;and construct k^2 array
kxv=shift(findgen(phisize(1))*2*!PI/(nx*dx),nx2/2)
kxv(0:nx2/2-1)=kxv(0:nx2/2-1)-2*!PI/(dx*nout_avg)
k2=kxv^2
if phisize(0) gt 2 then begin
    kyv=shift(findgen(phisize(2))*2*!PI/(ny*dy),ny2/2)
    kyv(0:ny2/2-1)=kyv(0:ny2/2-1)-2*!PI/(dy*nout_avg)
    k2=((fltarr(ny2)+1.)##kxv^2+kyv^2##(fltarr(nx2)+1.))
endif
if phisize(0) gt 3 and phisize(3) gt 1 then begin
    kzv=shift(findgen(phisize(3))*2*!PI/(nz*dz),nz2/2)
    kzv(0:nz2/2-1)=kzv(0:nz2/2-1)-2*!PI/(dz*nout_avg)
    k2store=k2
    k2=fltarr(nx2,ny2,nz2)
    for iz=0, nz2-1 do k2(*,*,iz)= k2store + kzv(iz)^2
endif

if phisize(0) eq 2 then E2kt=phi(*,ifirst:ilast/iskip)
if phisize(0) eq 3 then E2kt=phi(*,*,ifirst:ilast/iskip)
if phisize(0) eq 4 then begin
    E2kt_2d=fltarr(nx2,ny2,(ilast-ifirst+1)/iskip)
    nperp = ny2 > nz2
    bounds = Fltarr(4)
    if nz2 gt 0 then $
      dummy = int_theta(reform(kyv # kzv, ny2, nz2) ,kyv, kzv, /grid, $
                        points=nperp, bounds=bounds)
    k_perp_max = bounds[2]
endif

for i=ifirst,ilast,iskip do begin
    t=nout*dt*i
    if phisize(0) eq 2 then begin
        phifft=shift(fft(reform(phi(*,i))),nx2/2)
        E2kt(*,i)=float(phifft*conj(phifft))*k2
    endif else if phisize(0) eq 3 then begin
        phifft=shift(fft(reform(phi(*,*,i))),nx2/2,ny2/2)
        E2kt(*,*,i)=float(phifft*conj(phifft))*k2
    endif else if phisize(0) eq 4 then begin
        phifft=shift(fft(reform(phi(*,*,*,i))),nx2/2,ny2/2,nz2/2)
        E2k_3d=float(phifft*conj(phifft))*k2
        ; Integrate around theta to compress the perp dimension
        if nz2 gt 0 then $        
          FOR ix = 0, nx2-1 DO BEGIN
            E2kt_2d(ix,*,(i-ifirst)/iskip) = $
              int_theta(Reform(E2k_3d[ix,*,*]) $
                        ,kyv, kzv, /grid, /recall)
          ENDFOR
    endif
endfor    

end
