;    Find the maximum amplitude of the spectra at a fixed k for each
;    angle and plot it.

normalize_spectra=0
densize=size(den0)
if n_elements(den_out_subcycle0) eq 0 then den_out_subcycle0=1

for i=0, n_elements(wavelength)-1 do begin
    kwave=2*!PI/wavelength[i]
    k_fixed=fixed_k_spectra( den0, dx*nout_avg, dy*nout_avg, kwave )


    nomega=(size(k_fixed))[2]
    nangle=(size(k_fixed))[1]
    angle_vec=360./nangle * findgen(nangle)

    max_k_fixed=max(k_fixed(*,nomega/2:nomega-1),dimension=2)
    
    max_max=max(max_k_fixed)

    plot, angle_vec, 10*alog10(max_k_fixed), $
      yrange=10*alog10([max([max_max*1e-3,min(max_k_fixed)]), max_max]), $
      title='Power of Spectra at fixed wavelength,'+Strcompress(String(wavelength[i]))+' m', $ 
      xtitle = 'Angle in Degrees', ytitle='Power n(k,w)^2'
endfor

date_plot,title
end
