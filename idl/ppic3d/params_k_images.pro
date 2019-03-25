; Make some figures of parameters which are a function of fixed k
; values for a range of angles:

if (n_elements(den0) le 2) then begin
@den3d
den1=0
endif

densize=size(den0)
if n_elements(den_out_subcycle0) eq 0 then den_out_subcycle0=1

; We need a vector of all prefered wavelengths:
npts=max([nx,ny,nz])/nout_avg/2-1
if nz eq 1 then dz=0
k_want=(findgen(npts)+1.)/float(npts)*2*!PI/(sqrt( (dx)^2+(dy)^2+(dz)^2 )*(nout_avg*2))

k_max=params_k_spectra(den0, dx*nout_avg ,dy*nout_avg, dz*nout_avg, dt*nout, k_want,$
                       k_max_angle=k_max_angle, k_max_omega=k_max_omega, $
                       k_avg=k_avg)

;Plot K_max
plot, k_want, k_max, /ylog,  yrange=max(k_max)*[1E-3,1.],  $
  title='Max Spectral Energy all k', xtitle = 'k (1/m)', ytitle='|E|'
oplot, k_want, k_avg
date_plot,title



end
