function divergence, x, y, z, dx=dx, dy=dy, dz=dz
;
; NAME:
;       DIVERGENCE
;
; PURPOSE:
;       To calculate the divergence of a vector quantity on a uniform
;       mesh.  Works on up to a 3-D array.
;


if (n_elements(dx) eq 0) then dx=1.
if (n_elements(dy) eq 0) then dy=1.
if (n_elements(dz) eq 0) then dz=1.

xsize=size(x)

if xsize(0) eq 1 then begin
    div=(shift(x,-1)-shift(x,+1))/(2*dx)
endif

if xsize(0) eq 2 then begin
    div=(shift(x,-1,0)-shift(x,+1,0))/(2*dx)

    if (n_elements(y) eq n_elements(x)) then begin
        div=div+(shift(y,0,-1)-shift(y,0,+1))/(2*dy)
    endif else if (n_elements(y) ne 0) then message, $
      'Divergence error: y vector not same size as x vector'
endif

if xsize(0) eq 3 then begin
    div=(shift(x,-1,0,0)-shift(x,+1,0,0))/(2*dx)

    if (n_elements(y) eq n_elements(x)) then begin
        div=div+(shift(y,0,-1,0)-shift(y,0,+1,0))/(2*dy)
    endif else if (n_elements(y) ne 0) then message, $
      'Divergence error: y vector not same size as x vector'

    if (n_elements(z) eq n_elements(x)) then begin
        div=div+(shift(z,0,0,-1)-shift(z,0,0,+1))/(2*dz)
    endif else if (n_elements(z) ne 0) then message, $
      'Divergence error: z vector not same size as x vector'
endif

if xsize(0) ge 4 then message, 'Error: divergence function only works on 1D, 2D, or 3D arrays'

return, div
end
