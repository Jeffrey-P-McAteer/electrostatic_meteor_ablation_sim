; Plot E^2(k)

st = systime(1)
rect_time = 0.0
plot_time = 0.0
fft_time  = 0.0

@params_in.pro

ndecades = 3
; plot at most that many pages of spectra for fixed t
npages_max=3
; plot only spectra for fixed k_para with a maximum greater than minmax
minmax_para = 0.0025
; plot only spectra for fixed k_perp with a maximum greater than minmax
minmax_perp = 0.0025

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
if n_elements(fdir) eq 0 then fdir=''

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
        FILENAME=fdir + (Reverse(Str_Sep(curr,'/')))[0] + $
	'.ppw.' + StrTrim(String(nphi),2) + '.ps'
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

;Calculate the simulation volume
vol=nx2*dx
if ny gt 1 then vol=vol*ny2*dy
if nz gt 1 then vol=vol*nz2*dz

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

Print, 'time spent for initialization:', SysTime(1) - st

;== compute the transformations =======================================

; iskip_phi indicates the number of phi components to skip (1 means no skipping)
if n_elements(iskip_phi) EQ 0 THEN iskip_phi=1
nphi2=nphi/iskip_phi

hann = hanning(nphi2)

st_tot = SysTime(1)

st = SysTime(1)

print,'Allocating a ComplexArr(nx2,ny2,nz2,nphi) of size', nx2*ny2*nz2*nphi2*2

phik = ComplexArr(nx2,ny2,nz2,nphi2, /nozero)
FOR i=0, nphi2-1 DO BEGIN
    phik[*,*,*,i] = $
      FFT(Transpose(convert_endian(phi_a[i*iskip_phi]),[2,1,0]))/vol
ENDFOR
Print, 'time spent for FFT(x):        ', SysTime(1) - st

st = SysTime(1)
FOR ix = 0, nx2-1 DO BEGIN
	 FOR iy = 0, ny2-1 DO BEGIN
	     FOR iz = 0, nz2-1 DO BEGIN
		 phik[ix,iy,iz,*] = $
                   FFT(hann * phik[ix,iy,iz,*], /INVERSE)/ $
                   (nphi2*dt*nout*iskip_phi)
	     ENDFOR
	 ENDFOR
ENDFOR
Print, 'time spent for FFT(t):        ', SysTime(1) - st

st = SysTime(1)
; ignore negative frequencies
E2kppw = FltArr(nx2, nperp, nphi2/2, /nozero)
FOR it = 0, nphi2/2-1 DO BEGIN
	 E2kw = Float( k2 * phik[*,*,*,it] * Conj(phik[*,*,*,it]) )
	 FOR ix = 0, nx2-1 DO BEGIN
	     E2kppw[ix,*,it] = int_theta($
				   Reform(E2kw[ix,*,*], ny2, nz2)$
			       ,ky, kz, /grid, /recall)
	 ENDFOR
ENDFOR
Print, 'time spent for computing E^2: ', SysTime(1) - st

Print, '---- total computing time ----', SysTime(1) - st_tot

help, /memory

gmax = max(E2kppw)

;== plot E^2(k_par, kerp, t) for some k_par ===========================
IF !D.Name EQ 'X' THEN Window, 2

;For a 1D data set:
if (size(reform(e2kppw)))[0] eq 2 then begin
    maxi = Max( E2kppw ) / gmax
    title = StrCompress('k!D'+st_par+'!N =' $
                        + String( kx[0] ) + ', n!D' $
                        + ', Maximum =' $
                        + String(maxi,FORMAT='(g30.5)') )
    
    imageaxis, reform(e2kppw), /zlog, decades=ndecades, /zoom, title=title, $
      x0 = 0.0, dx = dkperp, xtitle = xtitle, $
      y0 = 0.0, dy = domega, ytitle = ytitle, $
      position = [px0+(i mod 2)*pxoff,     py0+(i/2)*pyoff, $
                  px0+(i mod 2)*pxoff+pxs, py0+(i/2)*pyoff+pys], $
      color=255

endif else begin

i = 0

FOR il = 0, nx2-1 DO BEGIN

    id = (il + nx2/2+1) mod nx2

    data =  Reform( E2kppw[id,*,*] )
    maxi = Max( data ) / gmax

    IF (maxi GT minmax_para OR kx[id] EQ 0.) THEN BEGIN

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

	n_par = ROUND(kx[id] / kx[1])
	
	title = StrCompress('k!D'+st_par+'!N =' $
			+ String( kx[id] ) $
			+ ', n!D' + st_par + '!N =' $
			+ String(n_par) $
			+ ', Maximum =' $
			+ String(maxi,FORMAT='(g30.5)') )
	
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


;== plot E^2(k_par, kerp, t) for some k_perp ===========================
IF !D.Name EQ 'X' THEN Window, 0

i = 0

FOR id = 0, nperp-1 DO BEGIN

    data =  Shift( Reform( E2kppw[*,id,*] ), nx2/2-1, 0)
    maxi = Max( data ) / gmax

    IF (maxi GT minmax_perp OR kx[id] EQ 0.) THEN BEGIN

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
	
	title = StrCompress('k!D' + st_perp + '!N =' $
			+ String(id * k_perp_max / nperp) $
			+ ', n!D' + st_perp + '!N =' $
			+ String(id) $
			+ ', Maximum =' $
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

;== compute the fft in x again ========================================
FOR i=0, nphi2-1 DO BEGIN
    phik[*,*,*,i] = $
      FFT(Transpose(convert_endian(phi_a[i*iskip_phi]),[2,1,0]))/vol
ENDFOR


;== plot E^2(k_par, kerp, t) for some t ===============================
IF !D.Name EQ 'X' THEN Window, 1
npages=Min([npages_max, nphi])
to_plot_long = Round(Reform(FIndgen(8*npages)/(8*npages-1)*(nphi2-1),8,npages))

data = Fltarr(nx2,nperp)

FOR p=0, N_Elements(to_plot_long)/8-1 DO BEGIN

    to_plot = to_plot_long[*,p]

    Erase

    FOR i=0, 7 DO BEGIN
	IF i  /  2 EQ 3 THEN BEGIN
		xtitle='k!D'+st_par+'!N [m!U-1!N]'
	ENDIF ELSE BEGIN
		xtitle=''
	ENDELSE
	IF i Mod 2 EQ 0 THEN BEGIN
		ytitle = 'k!D' + st_perp + '!N [m!U-1!N]'
	ENDIF ELSE BEGIN
		ytitle=''
	ENDELSE
	
	it = to_plot[i]
	E2k  = Float( k2 * phik[*,*,*,it] * Conj(phik[*,*,*,it]) )
	FOR ix = 0, nx2-1 DO BEGIN
	data[ix,*] = int_theta(Reform(E2k[ix,*,*], ny2, nz2) $
                               ,ky, kz, /grid, /recall)
	ENDFOR

	data = Shift( data, nx2/2-1, 0 )

	maxi = Max( data )
	title = StrCompress('t =' $
			+ String( to_plot[i] * dt * nout * iskip_phi) $
			+ ', Max =' $
			+ String(maxi/gmax,FORMAT='(g30.5)') )
	
	dkperp = k_perp_max / nperp
	dkpara = 2.*!pi / (dx * nx)

	imageaxis, $
		data, /zlog, decades=ndecades, /zoom, title=title, $
		x0 = min(kx), dx = dkpara, xtitle = xtitle, $
		y0 = 0.0,     dy = dkperp, ytitle = ytitle, $
		position = [px0+(i mod 2)*pxoff,     py0+(i/2)*pyoff, $
			    px0+(i mod 2)*pxoff+pxs, py0+(i/2)*pyoff+pys], $
		color=255
	
	ENDFOR

	Xyouts, 0.5, 0.97, alignment=0.5,/normal, charsize=1.2, color=255, $
		'E!U2!N, ' + (Reverse(Str_Sep(curr,'/')))[0] $
		+ ', logarithmic scale, ' + StrTrim(String(ndecades),2) $
		+ ' decades'
	XYOutS, 0.02, 0.0, /normal, charsize=0.4, color=255, info

ENDFOR

ENDELSE

;== done with plotting ================================================

IF !D.Name EQ "PS" THEN Device, /CLOSE

Free_Lun, unit

END
