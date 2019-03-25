function poisson_invert, b, dx=dx, dy=dy, dz=dz
;
; NAME: POISSON-INVERT
;
; PURPOSE: To calculate p in the equation del^2(p)=b where b is given.
;          Assume periodic boundary conditions. 
;          Return a complex array
;

bsize=size(b)
nx=bsize(1)
if (n_elements(dx) eq 0) then dx=1.
if (n_elements(dy) eq 0) then dy=1.
if (n_elements(dz) eq 0) then dz=1.

kx=fltarr(nx)
kx(0:nx/2)=findgen(nx/2+1)*2*!PI/(dx*nx)
kx(nx/2+1:nx-1)=-(nx/2-1-findgen(nx/2-1))*2*!PI/(dx*nx)
if (bsize(0) eq 1) then k2=kx^2

if (bsize(0) ge 2) then begin
    ny=bsize(2)
    ky=fltarr(ny)
    ky(0:ny/2)=findgen(ny/2+1)*2*!PI/(dy*ny)
    ky(ny/2+1:ny-1)=-(ny/2-1-findgen(ny/2-1))*2*!PI/(dy*ny)
    if (bsize(0) eq 2) then begin
        k2=((fltarr(ny)+1.)##kx^2+ky^2##(fltarr(nx)+1.))
    endif
endif

if (bsize(0) eq 3) then begin
    nz=bsize(2)
    kz=fltarr(nz)
    kz(0:nz/2)=findgen(nz/2+1)*2*!PI/(dz*nz)
    kz(nz/2+1:nz-1)=-(nz/2-1-findgen(nz/2-1))*2*!PI/(dz*nz)

    k2store=((fltarr(ny)+1.)##kx^2+ky^2##(fltarr(nx)+1.))
    k2=fltarr(nx,ny,nz)
    for iz=0, nz-1 do k2(*,*,iz)= k2store + kz(iz)^2
        
endif

k2(0)=1.                        ;Avoid dividing by 0

bfft=fft(b)
bfft=bfft/(-k2)
bfft(0)=0.                      ;Set the DC component to 0
bfft=fft(bfft,/inverse,/overwrite)
    
    



return, bfft
end
