function laplace, p, dx=dx, dy=dy, dz=dz
;
; NAME:
;       LAPLACE
;
; PURPOSE:
;       To calculate the LAPLACIAN of a scaler quantity on a uniform
;       mesh.  Works on up to a 3-D array.
;
if (n_elements(dx) eq 0) then dx=1.
if (n_elements(dy) eq 0) then dy=1.
if (n_elements(dz) eq 0) then dz=1.

psize=size(p)

if psize(0) eq 1 then begin
    v=(shift(p,-1)+shift(p,+1)-2*p)/(dx^2)
endif

if psize(0) eq 2 then begin
    v=(shift(p,-1,0)+shift(p,+1,0)-2*p)/(dx^2)+$
      (shift(p,0,-1)+shift(p,0,+1)-2*p)/(dy^2)
endif

if psize(0) eq 3 then begin
    v=(shift(p,-1,0,0)+shift(p,+1,0,0)-2*p)/(dx^2)
    v=v+(shift(p,0,-1,0)+shift(p,0,+1,0)-2*p)/(dy^2)
    v=v+(shift(p,0,0,-1)+shift(p,0,0,+1)-2*p)/(dz^2)
endif

return,v
end
