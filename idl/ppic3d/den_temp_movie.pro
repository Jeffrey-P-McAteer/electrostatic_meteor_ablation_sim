; This routine makes a movie of density adjacent to temp
; thermal instabilities
@den3d
@flux3d
@nvsqr3d

vx0= fluxx0/((1+den0)*n0d0)
if n_elements(fluxy0) gt 2 then vy0= fluxy0/((1+den0)*n0d0)
if n_elements(fluxz0) gt 2 then vz0= fluxz0/((1+den0)*n0d0)

vsqrx0= nvsqrx0/((1+den0)*n0d0) - vx0^2
if n_elements(nvsqry0) gt 2 then vsqry0= nvsqry0/((1+den0)*n0d0) - vy0^2
if n_elements(nvsqrz0) gt 2 then vsqrz0= nvsqrz0/((1+den0)*n0d0) - vz0^2
T0=vsqrx0(*,*,densize[3]/2,*)+vsqry0(*,*,densize[3]/2,*)

if ndist ge 1 then begin
  vx1= fluxx1/((1+den1)*n0d1)
  if n_elements(fluxy1) gt 2 then vy1= fluxy1/((1+den1)*n0d1)
  if n_elements(fluxz1) gt 2 then vz1= fluxz1/((1+den1)*n0d1)

  vsqrx1= nvsqrx1/((1+den1)*n0d1) - vx1^2
  if n_elements(nvsqry1) gt 2 then vsqry1= nvsqry1/((1+den1)*n0d1) - vy1^2
  if n_elements(nvsqrz1) gt 2 then vsqrz1= nvsqrz1/((1+den1)*n0d1) - vz1^2
  T1=vsqrx1(*,*,densize[3]/2,*)+vsqry1(*,*,densize[3]/2,*)
endif

image_movie, [(den1-min(den1))/max(den1),(T1-min(T1))/max(T1)]


end


