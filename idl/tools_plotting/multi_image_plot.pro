pro multi_image_plot, a, x, y, $              
                FIRST = first, LAST = last, SKIP = skip, $
                SUBTITLES = title_list, LEGEND=legend, TITLE=title, $
                NLABELS=nlabels, ZLOG=zlog, decades=decades, ZRANGE=zrange, $
                _EXTRA = e
;
; NAME:
;	IMAGE_PLOT
;
; PURPOSE:
;	Make a series of images and plot them with an axis, To add contour
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
;	A:	The tree-dimensional array to display. 
;               Displays a series of A(*,*,i) for all i.
;       x:      A vector of x values
;       y:      A vector of y values
;
; KEYWORD PARAMETERS:
;       LEGEND: If this keyword is set a colormap legend is created on
;               the right of the image.
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

on_error,2                      ;Return to caller if an error occurs
sz = size(a)			;Size of image
if sz(0) ne 3 then message, 'Parameter not 3D'

if n_elements(x) eq 0 then x=findgen(sz(1))
if n_elements(y) eq 0 then y=findgen(sz(2))

;Since X and Y point to the beginning of each element, adjust these
;for the contour routine:
dx=x(1)-x(0)
xrange=[x(0), x(n_elements(x)-1)+dx]
dy=y(1)-y(0)
yrange=[y(0), y(n_elements(y)-1)+dy]
xv=findgen(n_elements(x))/(n_elements(x)-1)*(xrange(1)-xrange(0))+xrange(0)
yv=findgen(n_elements(y))/(n_elements(y)-1)*(yrange(1)-yrange(0))+yrange(0)

if (n_elements(nlabels) eq 0) then nlabels=6
if (nlabels gt 1) then ndecades=nlabels-1 else ndecades=5
if n_elements(decades) eq 1 then begin
    ndecades=decades
    nlabels=ndecades+1
endif 

if n_elements(zrange) eq 0 then begin
    zrange=dblarr(2)
    zrange(0)=min(a)
    zrange(1)=max(a)
endif
if keyword_set(zlog) then begin
    zrange(0)=max([zrange(1)/10.^ndecades,zrange(0)])
endif

if n_elements(first) eq 0 then first=0
if n_elements(last) eq 0 then last=sz(3)-1
if n_elements(skip) eq 0 then skip=1

;subtitles
if n_elements(title_list) eq 0 then title_list=strarr(last+1)

;Set the first image on a new page
!p.multi(0)=0
npage=0

;Loop through each of the images:

for image=first,last,skip do begin

;First image on the page
    if (!p.multi(0) mod (!p.multi(1)*!p.multi(2)) eq 0) then begin
        npage=npage+1

; Set up some plot-wide fields
        p_store=!p.multi
        !p.multi=0
        !x.omargin=[0,0]
        !y.omargin=[0,0]
        xmargin=!X.MARGIN
;leave extra room for a legend if requested:
        charfactor=1
        if max([p_store,p_store(2)]) gt 2 then charfactor=2
        if keyword_set(legend) then xmargin=[xmargin(0),10*charfactor]

;set window used by contour
        contour,[[0,0],[1,1]], xrange, yrange ,/nodata, $
          xstyle=4, ystyle = 4, xmargin=xmargin, title=title, _EXTRA = e 
    
        if n_elements(title) ne 0 then !y.omargin=[0,2*charfactor]
        
        px = !x.window * !d.x_vsize ;Get size of window in device units
        py = !y.window * !d.y_vsize
        swx = px(1)-px(0)       ;Size in x in device units
        swy = py(1)-py(0)       ;Size in Y
        six = float(sz(1))      ;Image sizes
        siy = float(sz(2))
        aspi = six / siy        ;Image aspect ratio
        aspw = swx / swy        ;Window aspect ratio
        f = aspi / aspw         ;Ratio of aspect ratios
; Legend
        if keyword_set(legend) then begin
            rot=1
            if !d.name eq 'PS' then rot=3
            xloc=!x.region(1)*!d.x_size*(0.975)
            color_bar_image, xpos=xloc, thick =1/2., ypos=py(0), rot=rot
            for i=0, nlabels-1 do begin
                v=(zrange(1)-zrange(0))*float(i)/(nlabels-1)+zrange(0)
                if keyword_set(zlog) then $
                  v=exp((alog(zrange(1))-alog(zrange(0)))*float(i)/ $
                        (nlabels-1)+alog(zrange(0)))
                vexp=fix(alog10(v))
                str=[strcompress(string(v/10^(vexp*1.-1), $
                                 format='(F6.4)'), /remove_all), $
                    'X10^'+strcompress(string(vexp-1, $
                                format='(I3)'), /remove_all)]
                xyouts,xloc - !D.X_CH_SIZE/2., $
                  (py(1)-py(0))*float(i)/(nlabels-1)+py(0),str(0),$
                  /device, alignment=1.0, charsize=!p.charsize/2
                xyouts,xloc - !D.X_CH_SIZE/2., $
                  (py(1)-py(0))*float(i)/(nlabels-1)+py(0)-!D.Y_CH_SIZE, $
                  str(1), /device, alignment=1.0, charsize=!p.charsize/2
            endfor 
            !x.omargin=[0,10*charfactor]
        endif
        
        date_plot,' page '+strcompress(string(npage))

        !p.multi=p_store
        !p.multi(0)=max([!p.multi(1)*!p.multi(2),1])


    endif

    ;Generate image
    image_plot, reform(a(*,*,image)), x, y, title=title_list(image), $
      ZRANGE=zrange, NLABELS=nlabels, charsize=1.0, ZLOG=zlog, _EXTRA=e   
endfor

!x.omargin=[0,0]
!y.omargin=[0,0]
return
end
