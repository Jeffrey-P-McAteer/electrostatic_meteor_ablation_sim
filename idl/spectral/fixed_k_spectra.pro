; This routine takes a 3-D data set and analyzes it at a fixed set
;  of wavelengths for:
;      Phase Velocity
;      Spectral Width
;      Spectral Amplitude
;
; INPUT
;      den      a 3-D array of data in real Space
;      k_want  	the wave-number to interpolated to
;		mode_want=k_want*k_min
;
; Output
;      fixed_k_spectra    a 2-D real array of size (nx+ny, nt) 
;                         normalized to the peak of the 3D fft of den
;

function fixed_k_spectra, den, dx, dy, dz, k_want

dsize=size(reform(den))
if dsize[0] eq 4 then begin
    denfft=reform(den[*,*,*,*])
    dsize=size(reform(denfft))
endif else if dsize[0] eq 3 then denfft=reform(den) $
           else message, 'Input array not required 3-D'

nx=dsize[1]
ny=dsize[2]
kx_min=2*!PI/(nx*dx)
ky_min=2*!PI/(ny*dy)
k_min=sqrt(kx_min^2+ky_min^2)

;Estimate the number of unique points around 360 degree circle
isize=fix(k_want/min([kx_min,ky_min]))*16

if isize gt 0 then begin
    theta = dindgen(isize)/(isize)*2*!pi
    ikx = cos(theta)*k_want/kx_min + nx/2
    iky = sin(theta)*k_want/ky_min + ny/2

;Put a Hanning window over the data
    H_time=Hanning(dsize[3]/2)
    for i=0,dsize[3]/2-1 do begin
        denfft(*,*,i)=denfft(*,*,i)*H_time(i)
    endfor

;Fourier transform the input matrix:
    denfft=shift(abs(fft(denfft)), dsize[1]/2, dsize[2]/2, dsize[3]/2)
    
;Since the time transform should be a negative transform, swap the
;time dimension:
    nt=dsize[dsize[0]]
    for i = 0, nt/2-1 do begin
        temp=denfft(*,*,i)
        denfft(*,*,i)=denfft(*,*,nt-i-1)
        denfft(*,*,nt-i-1)=temp
    endfor

;Zero the D.C. component
    denfft(dsize[1]/2, dsize[2]/2, dsize[3]/2)=0

;Normalize
    denfft=denfft/max(denfft)

;Loop through w (the third dim) and interpolate to the ikx, ikz grid
    
    k_fix=fltarr(isize,dsize[3])
    for iw=0, dsize[3] - 1 do begin
        k_fix[*,iw] = interpolate(denfft[*,*,iw], ikx, iky, missing=0.)
    end
    
endif else begin
    print, 'Requested an inappropriate wave number', k_want, ' from fixed_k_spectra'
    print, 'Continuing...'
    k_fix=0
endelse

return, k_fix

end

