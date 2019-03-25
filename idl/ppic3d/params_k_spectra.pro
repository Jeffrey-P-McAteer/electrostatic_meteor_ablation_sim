; This routine takes a 3-D data set and returns statistical values for
; fixed k
;
; INPUT
;      den      a 3-D array of data in real Space
;      k_want  	the wave-numbers to interpolated to
;		mode_want=k_want*k_min
;
; Output
;      k_max    a 1-D real array of size n_element(k_want)


function params_k_spectra, den, dx, dy, dz, dt, k_want, $
                           k_max_angle=k_max_angle, k_max_omega=k_max_omega, $
                           k_avg=k_avg

;Set up the output arrays
k_max=fltarr(n_elements(k_want))
k_max_angle=fltarr(n_elements(k_want))
k_max_omega=fltarr(n_elements(k_want))
k_avg=fltarr(n_elements(k_want))

dsize=size(reform(den))
if dsize[0] eq 4 then begin
    denfft=reform(den[*,*,0,*])
    dsize=size(reform(denfft))
endif else if dsize[0] eq 3 then denfft=reform(den) $
else message, 'Input array not required 3-D'

;clear memory

nx=dsize[1]
ny=dsize[2]
kx_min=2*!PI/(nx*dx)
ky_min=2*!PI/(ny*dy)
k_min=sqrt(kx_min^2+ky_min^2)

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

;Loop through each chosen k:

for i=0,n_elements(k_want)-1 do begin
;Estimate the number of unique points around 360 degree circle
    isize=fix(k_want[i]/min([kx_min,ky_min]))*16
        
    if isize gt 0 then begin
        theta = dindgen(isize)/(isize)*2*!pi
        ikx = cos(theta)*k_want[i]/kx_min + nx/2
        iky = sin(theta)*k_want[i]/ky_min + ny/2
        
;Loop through w (the third dim) and interpolate to the ikx, ikz grid
        
        k_fix=fltarr(isize,dsize[3])
        for iw=0, dsize[3] - 1 do begin
            k_fix[*,iw] = interpolate(denfft[*,*,iw], ikx, iky, missing=0.)
        end
        
                                ; Extract parameters
        k_max[i]=max(k_fix,index)
        k_max_angle[i]=theta[(calc_index(k_fix,index))[0]]
        k_max_omega[i]=(calc_index(k_fix,index))[1]*(!PI/(nt*dt)/dsize[3])
        k_avg[i]=sqrt(mean(k_fix^2))

    endif else begin
        print, 'Requested an inappropriate wave number', k_want[i], ' from fixed_k_spectra'
        print, 'Continuing...'
        k_fix=0
    endelse
    
endfor

return, k_max

end

