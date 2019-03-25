; A routine to fft den2 into w,k space, integrate around k_perp and
; plot it


; Plotting

ndecades = 4
; plot at most that many pages of spectra for fixed t
npages_max=3
; plot only spectra for fixed k_para with a maximum greater than minmax
minmax_para = 0.01
; plot only spectra for fixed k_perp with a maximum greater than minmax
minmax_perp = 0.01


;== some plot settings ================================================
px0   = 0.06
py0   = 0.74
pxs   = 0.43
pys   = 0.16
pxoff = 0.48
pyoff = -0.22

Cd,curr=curr
!P.Charsize=0.8

spawn, 'ls -l den2.bin', res
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
        FILENAME=(Reverse(Str_Sep(curr,'/')))[0] + $
	'.den2.' + StrTrim(String(nden2),2) + '.ps'
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

i = 0

FOR il = 0, nx2-1 DO BEGIN

    id = (il + nx2/2+1) mod nx2

    data =  Reform( den2ppw[id,*,*] )
    maxi = Max( data ) / gmax

    IF (maxi GT minmax_para OR kx[id] EQ 0.) THEN BEGIN

	IF i EQ 0 THEN BEGIN
	    erase
	    XYOutS, 0.5, 0.97, alignment=0.5,/normal, charsize=1.2, $
		color=255, $
		'den!U2!N, ' + (Reverse(Str_Sep(curr,'/')))[0] $
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
	domega = 2.*!pi / (nden2*nout*dt)

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

END
