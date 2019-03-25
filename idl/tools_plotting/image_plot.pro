pro image_plot, a, x, y, WINDOW_SCALE = window_scale, ASPECT = aspect, $
                INTERP = interp,  nlevels = nlevels, LEGEND=legend, $
                ZRANGE = zrange, ZLOG=zlog, NLABELS=nlabels, DECADES=decades, $
                ZOOM=zoom, TITLE=title, _EXTRA = e
;
; NAME:
;	IMAGE_PLOT
;
; PURPOSE:
;	Make an image and plot it with an axis, To add contour
;	lines, use "nlevels=n".  This is an enhanced version of
;	image_cont, which also produces reversed images in postscript. 
;
; CATEGORY:
;	General graphics.
;
; CALLING SEQUENCE:
;	IMAGE_PLOT, A, x, y
;
; INPUTS:
;	A:	The two-dimensional array to display.
;       x:      A vector of x values
;       y:      A vector of y values
;
; KEYWORD PARAMETERS:
; WINDOW_SCALE:	Set this keyword to scale the window size to the image size.
;		Otherwise, the image size is scaled to the window size.
;		This keyword is ignored when outputting to devices with 
;		scalable pixels (e.g., PostScript).
;
;	ASPECT:	Set this keyword to retain the image's aspect ratio.
;		Square pixels are assumed.  If WINDOW_SCALE is set, the 
;		aspect ratio is automatically retained.
;
;	INTERP:	If this keyword is set, bilinear interpolation is used if 
;		the image is resized.
;       LEGEND: If this keyword is set a colormap legend is created on
;               the right of the image.
;       ZLOG:   Take the log of the amplitude (only 12 orders of
;                magnitude in range allowed)
;       ZOOM:   Use only the interesting part of the data
;
; OUTPUTS:
;	No explicit outputs.
;
; COMMON BLOCKS:
;	None.
;
; SIDE EFFECTS:
;	The currently selected display is affected.
;
; RESTRICTIONS:
;	None.
;
; PROCEDURE:
;	If the device has scalable pixels, then the image is written over
;	the plot window.
;
; MODIFICATION HISTORY:
;	DMS, May, 1988.
;-

;on_error,2                      ;Return to caller if an error occurs
a=reform(a)
sz = size(a)			;Size of image
if sz(0) ne 2  then message, 'Parameter not 2D'
; Test A for NaN 
if total(finite(a)) ne n_elements(a) then begin
    if total(finite(a)) ge n_elements(a)-1 then $
      message, 'Error: Only 1 or less finite value in array'
    print, 'Warning: Array contains non-finite values - ignoring them'
endif

if n_elements(x) eq 0 then x=findgen(sz(1))
if n_elements(y) eq 0 then y=findgen(sz(2))

xmargin=!X.MARGIN
;leave extra room for a legend if requested:
if keyword_set(legend) then xmargin=[xmargin(0),8]


;nlabels=ndecades+1 are equivelent - nlabels takes precedence.
if (n_elements(nlabels) eq 0) then $
  if (n_elements(ndecades) ne 0) then nlabels=ndecades+1 else nlabels = 4
ndecades=nlabels-1

if n_elements(zrange) eq 0 then begin
    amax=max(a,/nan)
    amin=min(a,/nan) 
endif else if n_elements(zrange) eq 1 then begin
    amax=zrange & amin=min(a,/nan) 
endif else if max(zrange) eq min(zrange) then begin
    amax=max(a,/nan)
    amin=min(a,/nan) 
endif else begin
    amax=max(zrange) & amin=min(zrange) 
endelse

; Using a logrithmic scale:
if not keyword_set(zlog) then begin
    b=a 
    bmax=amax
    bmin=amin
endif else begin
    b=alog(a > max([amax/10.^ndecades,amin],/nan))
    bmax=alog(amax)
    bmin=alog(max([amax/10.^ndecades,amin],/nan))
endelse

maxsize=!d.table_size
image=bytscl(b,top=maxsize-1,min=bmin,max=bmax) 

   ; Zoom into the interesting part of the data

IF Keyword_Set(zoom) THEN BEGIN
    index = where(b GT bmin, count)
    IF count GT 0 THEN BEGIN
        xmin   = min( index mod sz(1) )
        xmax   = max( index mod sz(1) )
        ymin   = min( index  /  sz(1) )
        ymax   = max( index  /  sz(1) )
        
        image  = image[xmin:xmax, ymin:ymax]
        b = b[xmin:xmax, ymin:ymax]
        x = x[xmin:xmax]
        y = y[ymin:ymax]
        sz=size(image)
    ENDIF
ENDIF

;Since X and Y point to the beginning of each element, adjust these
;for the contour routine:
dx=x(1)-x(0)
xrange=[x(0), x(n_elements(x)-1)+dx]
dy=y(1)-y(0)
yrange=[y(0), y(n_elements(y)-1)+dy]
xv=findgen(n_elements(x))/(n_elements(x)-1)*(xrange(1)-xrange(0))+xrange(0)
yv=findgen(n_elements(y))/(n_elements(y)-1)*(yrange(1)-yrange(0))+yrange(0)

;set window used by contour
contour,[[0,0],[1,1]], xrange, yrange $
  ,/nodata, xstyle=4, ystyle = 4, xmargin=xmargin, _EXTRA = e 

px = !x.window * !d.x_size	;Get size of window in device units
py = !y.window * !d.y_size
swx = px(1)-px(0)		;Size in x in device units
swy = py(1)-py(0)		;Size in Y
six = float(sz(1))		;Image sizes
siy = float(sz(2))
aspi = six*dx / (siy*dy)	;Image aspect ratio
aspw = swx / swy		;Window aspect ratio
f = aspi / aspw			;Ratio of aspect ratios

; If Postscript, reverse image
if !d.name eq 'PS' then image=!d.table_size-image-1

if (!d.flags and 1) ne 0 then begin ;Scalable pixels?
    if keyword_set(aspect) then begin ;Retain aspect ratio?
                                ;Adjust window size
        if f ge 1.0 then swy = swy / f else swx = swx * f
    endif

    tv,image,px(0),py(0),xsize = swx, ysize = swy, /device

endif else begin                ;Not scalable pixels	
    if keyword_set(window_scale) then begin ;Scale window to image?
        tv,image,px(0),py(0)    ;Output image
        swx = six		;Set window size from image
        swy = siy
    endif else begin		;Scale window
        if keyword_set(aspect) then begin
            if f ge 1.0 then swy = swy / f else swx = swx * f
        endif                   ;aspect
        tv,poly_2d(image,$      ;Have to resample image
                   [[0,0],[six/swx,0]], [[0,siy/swy],[0,0]],$
                   keyword_set(interp),swx,swy), $
          px(0),py(0)
    endelse			;window_scale
endelse                         ;scalable pixels

; Contours:
mx = !d.n_colors-1		;Brightest color
colors = [mx,mx,mx,0,0,0]	;color vectors
if !d.name eq 'PS' then colors = mx - colors ;invert line colors for pstscrp
if n_elements(nlevels) eq 0 then nlevels=0
if nlevels ne 0 then $
  contour,b, xv, yv,/noerase,/xst,/yst,$ ;Do the contour
  pos = [px(0),py(0), px(0)+swx,py(0)+swy],/dev,$
  c_color =  colors , title=title, _Extra = e $
else $
  contour,b, xv, yv,/noerase,/nodata,/xst,/yst,$ ;Do the contour
  pos = [px(0),py(0), px(0)+swx,py(0)+swy],/dev,$
  c_color =  colors, title=title, _Extra = e

; Legend
if !p.charsize eq 0 then !p.charsize = 1

if keyword_set(legend) then begin
    rot=1
    if !d.name eq 'PS' then rot=3
    xloc=!x.region(1)*!d.x_size+!d.x_ch_size*3
    color_bar_image, xpos=xloc, thick =1/2., ypos=py(0), rot=rot
    ywin=!y.window*!d.y_size

    for i=0, nlabels-1 do begin
        v=(bmax-bmin)*float(i)/(nlabels-1)+bmin
        if keyword_set(zlog) then v=exp(v)
        if (v ne 0) then vexp=fix(alog10(abs(v))) else vexp=0
        if vexp ge 0 and vexp lt 2 then begin
            str=strcompress(string(v,format='(F6.2)'),/remove_all)
            xyouts,xloc - !D.X_CH_SIZE, $
              (ywin(1)-ywin(0))*float(i)/(nlabels-1)+ywin(0),str,$
              /device, alignment=1.0, charsize=!p.charsize*0.75
        endif else begin
            str=[strcompress(string(v/10^(vexp*1.), $
                                    format='(F6.2)'), /remove_all), $
                 'E'+strcompress(string(vexp, $
                                           format='(I3)'), /remove_all)]
            xyouts,xloc - !D.X_CH_SIZE, $
              (ywin(1)-ywin(0))*float(i)/(nlabels-1)+ywin(0),str(0),$
              /device, alignment=1.0, charsize=!p.charsize*0.75
            xyouts,xloc - !D.X_CH_SIZE, $
              (ywin(1)-ywin(0))*float(i)/(nlabels-1)+ywin(0)-!D.Y_CH_SIZE, $
              str(1), /device, alignment=1.0, charsize=!p.charsize*0.75
        endelse
    endfor 
endif

return
end
