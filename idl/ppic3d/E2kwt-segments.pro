; Plot E^2(k,w) in time nseg segments.

st = systime(1)
rect_time = 0.0
plot_time = 0.0
fft_time  = 0.0

@params_in.pro

ndecades = 4
; plot at most that many pages of spectra for fixed t
npages_max=3
; plot only spectra for fixed k_para with a maximum greater than minmax
minmax_para = 0.5
; plot only spectra for fixed k_perp with a maximum greater than minmax
minmax_perp = 0.5

nperp = ny2 > nz2

;== open and assoc the data file ======================================
Openr, unit, 'phi.bin', /GET_LUN
phi_a = Assoc(unit, Fltarr(nz2,ny2,nx2))     
nphi = (Fstat(phi_a)).size / (4L * nz2*ny2*nx2 )            
    
;nphi = 32 < nphi
;Print, '### nphi changed to', nphi


;== some plot settings ================================================
px0   = 0.06
py0   = 0.74
pxs   = 0.43
pys   = 0.16
pxoff = 0.48
pyoff = -0.22

Cd,curr=curr
!P.Charsize=0.8

spawn, 'ls -l phi.bin', res
ls = str_sep(strcompress(res[0]), ' ')
spawn, 'date;uname -n;pwd', res
info = 'printed at ' + res[0] $
     + ', data created at ' + ls[5] + ' ' + ls[6]+' ' + ls[7] $
     + ' in ' + res[1] + ':' + res[2]

print, info

Loadct, 0
Tvlct, r, g, b, /GET
Tvlct, Reverse(r), Reverse(g), Reverse(b)

IF !D.Name EQ "PS" THEN BEGIN
    !P.Font=0
    Device, /LANDSCAPE, /HELVETICA, $
        /COLOR, BITS_PER_PIXEL=8, $
        FILENAME=('E2wk-t-segments' + '.ps')
ENDIF
IF !D.Name EQ "X" THEN BEGIN
    !P.Font=-1
ENDIF

IF !D.Name EQ 'PS' THEN BEGIN
	st_omega = '!9w!3'
	st_perp  = '!9^!3'
	st_par   = '||'
ENDIF ELSE BEGIN
	st_omega = '!7x!3'
	st_perp  = '!9x!3'
	st_par   = '!9#!3'
ENDELSE

;== prepare the transformation ========================================
pi2dx = !PI / (dx * NOUT_AVG)
pi2dy = !PI / (dy * NOUT_AVG)
pi2dz = !PI / (dz * NOUT_AVG)

kx = FIndgen(nx2) / nx2 * 2.*pi2dx
greater_pi = Where(kx GT pi2dx, count)
IF count GT 0 THEN kx[greater_pi] = kx[greater_pi] - 2.*pi2dx

ky = FIndgen(ny2) / ny2 * 2.*pi2dy
greater_pi = Where(ky GT pi2dy, count)
IF count GT 0 THEN ky[greater_pi] = ky[greater_pi] - 2.*pi2dy

kz = FIndgen(nz2) / nz2 * 2.*pi2dz
greater_pi = Where(kz GT pi2dz, count)
IF count GT 0 THEN kz[greater_pi] = kz[greater_pi] - 2.*pi2dz

;== initialize the transformation =====================================
bounds = Fltarr(4)
dummy = int_theta(reform(ky # kz, ny2, nz2) ,ky, kz, /grid, $
                     points=nperp, bounds=bounds)
k_perp_max = bounds[2]

k2 = Fltarr(nx2, ny2, nz2)
FOR ix=0, nx2-1 DO BEGIN
    FOR iy=0, ny2-1 DO BEGIN
	FOR iz=0, nz2-1 DO BEGIN
	    k2[ix,iy,iz] = kx[ix]^2 + ky[iy]^2 + kz[iz]^2
	ENDFOR
    ENDFOR
ENDFOR

;Set the default number of time segments

if (n_elements(nseg) eq 0) then nseg=8

Print, 'time spent for initialization:', SysTime(1) - st

;== compute the transformations =======================================

; iskip_phi indicates the number of phi components to skip (1 means no skipping)
if n_elements(iskip_phi) le 1 THEN begin 
    phisize=long(nx2)*ny2*nz2*nphi*8/nseg
    iskip_phi=long(phisize/6e8)+1
endif
if iskip_phi gt 1 then print, 'Skipping every ',iskip_phi,' phi arrays for space'
nphi2=nphi/iskip_phi/nseg

st_tot = SysTime(1)

st = SysTime(1)

;Each segment extends overlaps the next segment:
nphi_extent=nphi2*2

print,'Allocating a ComplexArr(nx2,ny2,nz2,nphi) of size', nx2*ny2*nz2*nphi_extent
phik = ComplexArr(nx2,ny2,nz2,nphi_extent, /nozero)

hann = hanning(nphi_extent)

; Loop through the segments
for iseg= 0, nseg-2 do begin

    istart=iseg*nphi2*iskip_phi
    FOR i=0, nphi_extent-1 DO BEGIN
        phik[*,*,*,i] = Transpose(convert_endian(phi_a[istart+i*iskip_phi]), [2,1,0])
        phik[*,*,*,i] = FFT(phik[*,*,*,i],/overwrite)
    ENDFOR
    Print, 'time spent for FFT(x):        ', SysTime(1) - st
    
    st = SysTime(1)
    FOR ix = 0, nx2-1 DO BEGIN
        FOR iy = 0, ny2-1 DO BEGIN
            FOR iz = 0, nz2-1 DO BEGIN
                phik[ix,iy,iz,*] = FFT(hann * phik[ix,iy,iz,*], /INVERSE)
            ENDFOR
        ENDFOR
    ENDFOR
    Print, 'time spent for FFT(t):        ', SysTime(1) - st
    
    st = SysTime(1)
; ignore negative frequencies
    E2kppw = FltArr(nx2, nperp, nphi_extent, /nozero)
    FOR it = 0, nphi_extent-1 DO BEGIN
        E2kw = Float( k2 * phik[*,*,*,it] * Conj(phik[*,*,*,it]) )
        FOR ix = 0, nx2-1 DO BEGIN
            E2kppw[ix,*,it] = int_theta(Reform(E2kw[ix,*,*], ny2, nz2) $
                                        ,ky, kz, /grid, /recall)
        ENDFOR
    ENDFOR
    Print, 'time spent for computing E^2: ', SysTime(1) - st
    
    Print, '---- total computing time ----', SysTime(1) - st_tot
    
    help, /memory

    gmax = max(E2kppw)

;== plot E^2(w, k_perp, t) for some k_par ===========================
    IF !D.Name EQ 'X' THEN Window, 2

    i = 0

    FOR il = 0, nx2-1 DO BEGIN

        id = (il + nx2/2+1) mod nx2
        
        data =  Reform( E2kppw[id,*,*] )
        maxi = Max( data ) / gmax

        n_par = ROUND(kx[id] / kx[1])

        IF (maxi GT minmax_para OR abs(n_par) LE 2) THEN BEGIN

            IF i EQ 0 THEN BEGIN
                erase
                XYOutS, 0.5, 0.97, alignment=0.5,/normal, charsize=1.2, $
                  color=255, $
                  'E!U2!N, ' + (Reverse(Str_Sep(curr,'/')))[0] $
                  + ', logarithmic scale, ' + StrTrim(String(ndecades),2) $
                  + ' decades'
                XYOutS, 0.02, 0.0, /normal, charsize=0.4, color=255, info
            ENDIF

            IF i  /  2 EQ 3 THEN BEGIN
                xtitle='k!D' + st_perp + '!N [m!U-1!N]'
            ENDIF ELSE BEGIN
                xtitle=''
            ENDELSE
            IF i Mod 2 EQ 0 THEN BEGIN
                ytitle = st_omega + ' [rad/s]'
            ENDIF ELSE BEGIN
                ytitle=''
            ENDELSE

            title = StrCompress('T= ' + String(iseg*nphi2*iskip_phi*dt*nout,format='(f5.0)')+ $
                                ' - ' + String((iseg+1)*nphi2*iskip_phi*dt*nout,format='(f5.0)') $
                                + 'k!D'+st_par+'!N =' + String( kx[id] ) $
                                + ', n!D' + st_par + '!N =' + String(n_par) $
                                + ', Max =' + String(maxi,FORMAT='(g30.5)') )
	
            dkperp = k_perp_max / nperp
            domega = 2.*!pi / (nphi*nout*dt)

            imageaxis, $
              data, /zlog, decades=ndecades, /zoom, title=title, $
              x0 = 0.0, dx = dkperp, xtitle = xtitle, $
              y0 = 0.0, dy = domega, ytitle = ytitle, $
              position = [px0+(i mod 2)*pxoff,     py0+(i/2)*pyoff, $
                          px0+(i mod 2)*pxoff+pxs, py0+(i/2)*pyoff+pys], $
              color=255

            i = (i + 1) mod 8
            
        ENDIF
        
    ENDFOR

;== plot E^2(w, k_par, t) for some k_perp ===========================
    IF !D.Name EQ 'X' THEN Window, 0

    i = 0

    FOR id = 0, nperp-1 DO BEGIN

        data =  Shift( Reform( E2kppw[*,id,*] ), nx2/2-1, 0)
        maxi = Max( data ) / gmax

        IF (maxi GT minmax_perp OR id LE 2) THEN BEGIN
            
            IF i EQ 0 THEN BEGIN
                erase
                Xyouts, 0.5, 0.97, alignment=0.5,/normal, charsize=1.2, $
                  color=255, $
                  'E!U2!N, ' + (Reverse(Str_Sep(curr,'/')))[0] $
                  + ', logarithmic scale, ' + StrTrim(String(ndecades),2) $
                  + ' decades'
                XYOutS, 0.02, 0.0, /normal, charsize=0.4, color=255, info
            ENDIF

            IF i  /  2 EQ 3 THEN BEGIN
                xtitle='k!D' + st_par + '!N [m!U-1!N]'
            ENDIF ELSE BEGIN
                xtitle=''
            ENDELSE
            IF i Mod 2 EQ 0 THEN BEGIN
                ytitle = st_omega + ' [rad/s]'
            ENDIF ELSE BEGIN
                ytitle=''
            ENDELSE
            
            title = StrCompress('T=' + String(iseg*nphi2*iskip_phi*dt*nout)+ $
                                ' - ' + String((iseg+1)*nphi2*iskip_phi*dt*nout) $
                                + 'k!D' + st_perp + '!N =' $
                                + String(id * k_perp_max / nperp) $
                                + ', n!D' + st_perp + '!N =' $
                                + String(id) $
                                + ', Max =' $
                                + String(maxi,FORMAT='(g30.5)') )
            
            dkpara = (max(kx) - min(kx)) / nx2
            domega = 2.*!pi / (nphi*nout*dt)
            
            imageaxis, $
              data, /zlog, decades=ndecades, /zoom, title=title, $
              x0 = min(kx), dx = dkpara, xtitle = xtitle, $
              y0 = 0.0,     dy = domega, ytitle = ytitle, $
              position = [px0+(i mod 2)*pxoff,     py0+(i/2)*pyoff, $
                          px0+(i mod 2)*pxoff+pxs, py0+(i/2)*pyoff+pys], $
              color=255
            
            i = (i + 1) mod 8
            
        ENDIF
        
    ENDFOR

ENDFOR


;== done with plotting ================================================

    IF !D.Name EQ "PS" THEN Device, /CLOSE

    Free_Lun, unit

END
