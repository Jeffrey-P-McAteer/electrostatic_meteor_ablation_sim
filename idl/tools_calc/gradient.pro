function gradient, p, dx=dx, dy=dy, dz=dz
;
; NAME:
;       GRADIENT
;
; PURPOSE:
;       To calculate the GRADIENT of a scaler quantity on a uniform
;       mesh.  Works on up to a 3-D array.
;

psize=size(p)

if (n_elements(dx) eq 0) then dx=1.
if (n_elements(dy) eq 0) then dy=1.
if (n_elements(dz) eq 0) then dz=1.

if psize(0) eq 1 then begin
    grad=(shift(p,-1)-shift(p,+1))/(2*dx)
    return, grad
endif

if psize(0) eq 2 then begin
    gradx=(shift(p,-1,0)-shift(p,+1,0))/(2*dx)
    grady=(shift(p,0,-1)-shift(p,0,+1))/(2*dy)
    return, [ [[gradx]], [[grady]]]
endif

if psize(0) eq 3 then begin

    ; Problem: IDL appears to limit the concatination depth
    ;          so we'll use the following technique:
    psize=size(p)
    grad_vec=make_array(psize(1),psize(2),psize(3),psize(0), $
                        type=psize(psize(0)+1))
    grad_vec(*,*,*,0)=(shift(p,-1,0,0)-shift(p,+1,0,0))/(2*dx)
    grad_vec(*,*,*,1)=(shift(p,0,-1,0)-shift(p,0,+1,0))/(2*dy)
    grad_vec(*,*,*,2)=(shift(p,0,0,-1)-shift(p,0,0,+1))/(2*dz)
    return, grad_vec
endif

end
