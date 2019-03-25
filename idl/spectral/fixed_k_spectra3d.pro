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
;      aspect_ang = angles with respect to z axis desired.
;
; Output
;      fixed_k_spectra    a 2-D real array of size (nx+ny, nt) 
;                         normalized to the peak of the 3D fft of den
;

function fixed_k_spectra3d, den, dx, dy, dz, k_want, aspect_ang=aspect_ang

dsize=size(reform(den))
if dsize[0] eq 4 then begin
    denfft=reform(den[*,*,*,*])
    dsize=size(reform(denfft))
endif else if dsize[0] eq 3 then denfft=reform(den) $
           else message, 'Input array not required 3-D'

nx=dsize[1]
ny=dsize[2]
nz=dsize[3]
nt=dsize[4]
kx_min=2*!PI/(nx*dx)
ky_min=2*!PI/(ny*dy)
if nz gt 1 then kz_min=2*!PI/(nz*dz) else kz_min=0
k_min=sqrt(kx_min^2+ky_min^2+kz_min^2)

;Estimate the number of unique points around 360 degree circle
isize=fix(k_want/min([kx_min,ky_min,kz_min]))*16

if isize gt 0 then begin

;Put a Hanning window over the data
    H_time=Hanning(nt/2)
    for i=0,nt/2-1 do begin
        denfft(*,*,*,i)=denfft(*,*,*,i)*H_time(i)
    endfor

;Fourier transform the input matrix:
    denfft=shift(abs(fft(denfft)), dsize[1]/2, dsize[2]/2, dsize[3]/2, dsize[4]/2)
    
;Since the time transform should be a negative transform, swap the
;time dimension:
    nt=dsize[dsize[0]]
    for i = 0, nt/2-1 do begin
        temp=denfft(*,*,*,i)
        denfft(*,*,*,i)=denfft(*,*,*,nt-i-1)
        denfft(*,*,*,nt-i-1)=temp
    endfor

;Zero the D.C. component
    denfft(dsize[1]/2, dsize[2]/2, dsize[3]/2, dsize[4]/2)=0

;Normalize
    denfft=denfft/max(denfft)


    if n_elements(aspect_ang) eq 0 then aspect_ang=0

; IN 3D it will interpolate along circle of length k_want offset from
; the x-y plane by the angle aspect_ang.  If aspect_ang=0, then it
; will interpolate in the x-y plane (one can more efficiently do the
; aspect_ang=0 case by passing a 2D array where one averages along z).

    nphi=n_elements(aspect_ang)
    k_fix=fltarr(nphi,isize,nt)
    for iphi=0, nphi-1 do begin
        phi=aspect_ang[iphi]
    
        theta = dindgen(isize)/(isize)*2*!pi
        ikx = sin(phi)*sin(theta)*k_want/kx_min + nx/2
        iky = sin(phi)*cos(theta)*k_want/ky_min + ny/2
        ikz = cos(phi)*(dblarr(isize)+1.0)*k_want/kz_min+nz/2

; Loop through w (the fourth dim) and interpolate to the ikx, ikz grid
    
        for iw=0, nt - 1 do begin
            k_fix[iphi,*,iw] = interpolate(denfft[*,*,*,iw], ikx, iky, ikz, missing=0.)
        endfor

    endfor
    
endif else begin
    print, 'Requested an inappropriate wave number', k_want, ' from fixed_k_spectra'
    print, 'Continuing...'
    k_fix=0
endelse

return, k_fix

end

