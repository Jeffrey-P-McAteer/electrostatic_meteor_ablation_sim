; A routine to fft den2 into w,k space, integrate around k_perp and
; plot it

@ppic3d.i
nx2=long(nx/nout_avg) > 1
if (n_elements(ny) eq 0) then ny=1 
ny2=long(ny/nout_avg) > 1 
if (n_elements(nz) eq 0) then nz=1 
nz2=long(nz/nout_avg) > 1
nsize=nx2*ny2*nz2

if (n_elements(den2) lt 2) then begin
    den2=readarray('den2.bin',[nx2,ny2,nz2],order=[2,1,0],/binary,skip=iskip)
endif

den2size=size(den2)
nden2=den2size(den2size(0))

;== prepare the transformation ========================================
pi2dx = !PI / (dx * NOUT_AVG)
pi2dy = !PI / (dy * NOUT_AVG)
if (nz eq 1 ) then  dz=1 
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

hann = hanning(nden2)

st_tot = SysTime(1)

st = SysTime(1)

print,'Allocating a ComplexArr(nx2,ny2,nz2,nden) having', nx2*ny2*nz2*nden2*2,$
  ' elements'

denk = ComplexArr(nx2,ny2,nz2,nden2, /nozero)
FOR i=0, nden2-1 DO BEGIN
    denk[*,*,*,i] = FFT(den2[*,*,*,i])
ENDFOR
Print, 'time spent for FFT(x):        ', SysTime(1) - st

st = SysTime(1)
FOR ix = 0, nx2-1 DO BEGIN
	 FOR iy = 0, ny2-1 DO BEGIN
	     FOR iz = 0, nz2-1 DO BEGIN
		 denk[ix,iy,iz,*] = FFT(hann * denk[ix,iy,iz,*], /INVERSE)
	     ENDFOR
	 ENDFOR
ENDFOR
Print, 'time spent for FFT(t):        ', SysTime(1) - st

st = SysTime(1)
nperp = ny2 > nz2
if nz gt 1 then begin
; ignore negative frequencies
    den2ppw = FltArr(nx2, nperp, nden2/2, /nozero)
    FOR it = 0, nden2/2-1 DO BEGIN
        den2kw = Float( denk[*,*,*,it] * Conj(denk[*,*,*,it]) )
        FOR ix = 0, nx2-1 DO BEGIN
            den2ppw[ix,*,it] = int_theta(Reform(den2kw[ix,*,*], ny2, nz2) $
                                         ,ky, kz, /grid, /recall)
        ENDFOR
    ENDFOR
endif else begin
    den2ppw=reform(Float( denk * Conj(denk) ) )
endelse

Print, 'time spent for computing E^2: ', SysTime(1) - st

Print, '---- total computing time ----', SysTime(1) - st_tot

help, /memory

gmax = max(den2ppw)

; Plotting

ndecades = 4
; plot at most that many pages of spectra for fixed t
npages_max=3
; plot only spectra for fixed k_para with a maximum greater than minmax
minmax_para = 0.0025
; plot only spectra for fixed k_perp with a maximum greater than minmax
minmax_perp = 0.0025

IF !D.Name EQ 'PS' THEN BEGIN
	st_omega = '!9w!3'
	st_perp  = '!9^!3'
	st_par   = '||'
ENDIF ELSE BEGIN
	st_omega = '!7x!3'
	st_perp  = '!9x!3'
	st_par   = '!9#!3'
ENDELSE

Cd, curr=curr
spawn, 'ls -l den2.bin', res
ls = str_sep(strcompress(res[0]), ' ')
spawn, 'date;uname -n;pwd', res
info = 'printed at ' + res[0] $
     + ', data created at ' + ls[5] + ' ' + ls[6]+' ' + ls[7] $
     + ' in ' + res[1] + ':' + res[2]


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
	domega = 2.*!pi / (nden2*nout*dt*den_out_subcycle2)

        px0   = 0.06
        py0   = 0.74
        pxs   = 0.43
        pys   = 0.16
        pxoff = 0.48
        pyoff = -0.22

	image_plot, $
		data, /zlog, decades=ndecades, /zoom, title=title, $
		x0 = 0.0, dx = dkperp, xtitle = xtitle, $
		y0 = 0.0, dy = domega, ytitle = ytitle
;		position = [px0+(i mod 2)*pxoff,     py0+(i/2)*pyoff, $
;			    px0+(i mod 2)*pxoff+pxs, py0+(i/2)*pyoff+pys], $
;		color=255

	i = (i + 1) mod 8

    ENDIF

ENDFOR

END
