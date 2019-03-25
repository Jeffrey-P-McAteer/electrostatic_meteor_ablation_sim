; Make some figures of parameters which are a function of fixed k
; values for a range of angles:

if (n_elements(phi) le 2) then begin
@phi3d
endif

; We want only the 2nd half of the data:
nphi=phisize[phisize[0]]
phi=phi[*,*,*,nphi/2:nphi-1]
phisize=size(phi)
nphi=phisize[phisize[0]]

E2=phi
; We want only the 2nd half of the data:
E2size=size(E2)
nE2=E2size[E2size[0]]
E2=E2[*,*,*,nE2/2:nE2-1]
E2size=size(E2)
nE2=E2size[E2size[0]]

if (phisize(3) eq 1) then begin
    for i=0,nE2-1 do begin
        E=gradient(reform(phi[*,*,0,i]),dx=dx*nout_avg,dy=dy*nout_avg)
        E2(*,*,0,i)=total(abs(E*conj(E)),3)
    endfor
endif

; We need a vector of all prefered wavelengths:
npts=max([nx,ny,nz])/nout_avg/2-1
if nz eq 1 then dz=0
k_want=(findgen(npts)+1.)/float(npts)*2*!PI/(sqrt( (dx)^2+(dy)^2+(dz)^2 )*(nout_avg*2))

;k_max=params_k_spectra(phi, dx*nout_avg ,dy*nout_avg, dz*nout_avg, dt*nout, k_want,$
;                       k_max_angle=k_max_angle, k_max_omega=k_max_omega, $
;                       k_avg=k_avg)

k_max=params_k_spectra(sqrt(E2), dx*nout_avg ,dy*nout_avg, dz*nout_avg, dt*nout, k_want,$
                       k_max_angle=k_max_angle, k_max_omega=k_max_omega, $
                       k_avg=k_avg)

;Plot K_max
;plot, k_want, k_max, /ylog,  yrange=max(k_max)*[1E-3,1.],  $
;  title='Max Spectral Energy from phi', xtitle = 'k (1/m)',
;  ytitle='phi(k)'
plot, k_want, k_max, /ylog,  yrange=max(k_max)*[1E-3,1.],  $
  title='Max Spectral Energy from E', xtitle = 'k (1/m)', ytitle='|E(k)|'
oplot, k_want, k_avg
date_plot,title

end
