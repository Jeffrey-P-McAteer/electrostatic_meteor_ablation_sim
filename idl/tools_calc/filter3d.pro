Function filter3d, Ain, fwidth=fwidth, fsteep=fsteep
;
; filter3d filters Ain, removing the high frequency components
;          by applying the following damping function to Ain
;          after it has been fourier transformed into k-space:
;          
;          Aout=Ain*exp(-(fwidth*k)^(2*fsteep))
;
; If Ain is 4D then this process will be done sequentially to each
; successive frame.
; 
;ON_ERROR,2              ;Return to caller if an error occurs

Asize=size(Ain)
if (Asize(0) lt 3 or Asize(0) gt 4) then $
  message, ' FILTER3D requires a 3D or 4D array'

if (n_elements(fwidth) eq 0) then fwidth=2
if (n_elements(fsteep) eq 0) then fsteep=3


;Set up damping matrix:
nx=Asize(1)
kx=fltarr(nx)
kx(0:nx/2)=fwidth*findgen(nx/2+1)*2*!PI/(nx)
kx(nx/2+1:nx-1)=-fwidth*(nx/2-1-findgen(nx/2-1))*2*!PI/(nx)

ny=Asize(2)
ky=fltarr(ny)
ky(0:ny/2)=fwidth*findgen(ny/2+1)*2*!PI/(ny)
ky(ny/2+1:ny-1)=-fwidth*(ny/2-1-findgen(ny/2-1))*2*!PI/(ny)

nz=Asize(3)
kz=fltarr(nz)
kz(0:nz/2)=fwidth*findgen(nz/2+1)*2*!PI/(nz)
if nz gt 1 then kz(nz/2+1:nz-1)=-fwidth*(nz/2-1-findgen(nz/2-1))*2*!PI/(nz)

;Build 3D damping array:
kdamp=make_array(nx,ny,nz,/float)
kdamp2d=((fltarr(ny)+1.)##kx^2+ky^2##(fltarr(nx)+1.))
for iz=0,nz-1 do kdamp(*,*,iz)=kdamp2d+kz(iz)^2
kdamp=exp(-1*kdamp^fsteep)

if (Asize(0) eq 3) then begin
    Aout=fft(Ain,-1)*kdamp
    Aout=fft(Aout,1,/overwrite)
endif else begin
    Aout=Ain
    for i=0,Asize(4)-1 do begin
        Aout2=fft(Aout(*,*,*,i),-1)*kdamp
        Aout(*,*,*,i)=fft(Aout2,1)
    endfor
endelse

return, Aout
end


