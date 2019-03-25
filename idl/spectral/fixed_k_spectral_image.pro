;pro, fixed_k_specral_image, wavelength
; Make some figures of spectra at fixed k for a range of angles:

densize=size(den0)
nden=densize[densize[0]]
if n_elements(den_out_subcycle0) eq 0 then den_out_subcycle0=1

if densize[0] ne 4 then message, 'den0 is not 3D as required (3rd dim can = 1)'

if n_elements(aspect_ang) eq 0 then aspect_ang=0
if densize[3] gt 1 then begin
    if n_elements(aspect_angle) eq 1 and aspect_ang[0] eq 0 then begin
        den0bck=den0
        den0=total(den0,1)/densize[1]
        densize=size(den0)
        den0=reform(den0,[densize[1],densize[2],1,densize[3]])
        den0=transpose(den0,[1,0,2,3]) ; to make it look similar to 2D
        den0=reverse(den0,1)
    endif
endif


for i=0, n_elements(wavelength)-1 do begin
    kwave=2*!PI/wavelength[i]
    
    if densize[3] eq 1 then $
      k_fixed=fixed_k_spectra2d( den0, dx*nout_avg, dy*nout_avg, kwave ) $
    else $
      k_fixed=fixed_k_spectra3d( den0, dx*nout_avg, dy*nout_avg, dz*nout_avg, kwave, aspect_ang=aspect_ang ) 
    
    kf_size=size(k_fixed)
    if kf_size[0] eq 2 then k_fixed=reform(k_fixed,[1,kf_size[1],kf_size[2]]) &kf_size=size(k_fixed)

    nomega=0
    nangle=0
    if n_elements(k_fixed) gt 0 then begin
        nomega=(size(k_fixed))[3]
        nangle=(size(k_fixed))[2]
    endif

    if (n_elements(normalize_spectra) eq 0) then normalize_spectra=0
    if (normalize_spectra eq 1) then $
    for ia=0, nangle-1 do k_fixed(*,ia,*) /= max(k_fixed(*,ia,*))

    omega_min=2.*!pi/(nden*nout*dt*den_out_subcycle0)
    omega_vec=omega_min*(findgen(nomega)-nomega/2)
    angle_vec=360./nangle * findgen(nangle)

; Put a known through as a test case
;    k_fixed(nangle-1,*)=1E-3*exp(-(omega_vec-500.*kwave)^2/(2*(1000.*kwave)^2))

    for iphi=0,n_elements(aspect_ang)-1 do begin

; Calculate the mean Omega
        mean_omega=k_fixed(iphi,*,0)
;
        max_k_fixed=max(k_fixed)

;   Calculate a noise floor - the noise effects the mean and variance
;                             because of the limited spectral range.
;   Set the noise so that 50% of the signal is captured.
;       noise=(k_fixed(sort(k_fixed)))[n_elements(k_fixed)*0.50]
;   Set the noise so that any signal below 2% of the peak is ignored.
        noise=k_fixed(iphi,*,0) ;set array size
        stddev_omega=k_fixed(iphi,*,0) ;set array size
        for ia=0,nangle-1 do begin
            noise[ia]=(k_fixed(iphi,ia,sort(k_fixed(iphi,ia,*))))[nomega*0.50]
;        noise[ia]=max(k_fixed(iphi,ia,*))*0.01
            dist=reform((k_fixed(iphi,ia,*) > noise[ia]) -noise[ia])
            mean_omega(ia)=mean(dist*omega_vec)/mean(dist)
                                ; Calculate the spectral width
            stddev_omega(ia)=sqrt(total(dist*(omega_vec-mean_omega(ia))^2)/total(dist))
            Amag=total(dist)*omega_min/(sqrt(2*!PI)*stddev_omega[ia])
;   plot,omega_vec/kwave,Amag*exp(-(omega_vec-mean_omega[ia])^2/(2*stddev_omega[ia]^2))
        endfor
        
                                ; A better way to get the
                                ; standard_deviation is to use a
                                ; Gaussioan fit with quadratitic terms
                                ; (see gaussfit)
        stddev2=fltarr(nangle)
        mean2=fltarr(nangle)
        chisq=fltarr(nangle)
        yerror=fltarr(nangle)
        peak=fltarr(nangle)
        
        for ia=0,nangle-1 do begin
            estimates=[Amag,mean_omega[ia]/kwave,stddev_omega[ia]/kwave,noise[ia],0.,0.]
            gfit=gaussfit2((omega_vec/kwave),reform(k_fixed(iphi,ia,*)),gcoeff, $
                           estimates=estimates, chisq=chisq2,yerror=yerror2, $
                           status=status,itmax=100,tol=max_k_fixed*1E-6, iter=iter) 
            stddev2[ia]=abs(gcoeff[2]) ; Sometimes the fit finds a negative stddev
            mean2[ia]=gcoeff[1]
            chisq[ia]=chisq2
            yerror[ia]=yerror2
            peak[ia]=gcoeff[0]
;            print, ia, status, iter, yerror2, mean2[ia], stddev2[ia]
            if status ne 0 then yerror[ia]=-1
        endfor
        
                                ;	Calculate the range of frequencies to display
        if n_elements(Ey0_external) ne 0 and n_elements(Bz) ne 0 then $
          vph_max=min([abs(Ey0_external/Bz*1.5),max(omega_vec/kwave)]) $
        else if n_elements(E0z) ne 0 and n_elements(B0) ne 0 then $
          vph_max=min([abs(E0z/B0*1.5),max(omega_vec/kwave)]) $
        else vph_max=max(omega_vec/kwave)
        
        y_index=where(omega_vec/kwave ge (-vph_max) and omega_vec/kwave le vph_max)
        
        if n_elements(k_fixed) gt 2 then begin
            image_plot, k_fixed(iphi,0:nangle-1,y_index), angle_vec, (omega_vec/kwave)[y_index] $
              , /legend $
              , title='Spectra at '+Strcompress(String(wavelength[i]))+' m; p='+Strcompress(String(aspect_ang[iphi]*180/!PI)) $
              , xtitle = 'Angle in Degrees', ytitle='Phase Velocity (m/s)'
            oplot, angle_vec, mean_omega/kwave, thick=2.0
            xyouts, 0.,(omega_vec/kwave)[y_index(1)] , $
              'max Vph='+Strcompress(string(max(mean_omega/kwave)))
            
            if n_elements(plot_stddev) eq 0 then plot_stddev=1
            if (plot_stddev) then $
                plot, angle_vec, stddev_omega/kwave, $
                  title='std dev. at '+Strcompress(String(wavelength[i]))+' m' $
                  , xtitle = 'Angle in Degrees', ytitle='std dev (m/s)'
            ; Limit to "good" fits
            good_fits=where(yerror gt 0.00)
            if n_elements(good_fits) gt 1 then $
              oplot, angle_vec(good_fits), stddev2(good_fits), color=120 ;& $
              ;oplot, angle_vec(good_fits), mean2(good_fits), color=80
                
        endif
            
    endfor
endfor

if n_elements(den0bck) gt n_elements(den0) then den0=den0bck & den0bck=1

date_plot,title

end
