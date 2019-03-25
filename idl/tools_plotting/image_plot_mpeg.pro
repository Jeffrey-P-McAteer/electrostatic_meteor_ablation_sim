pro image_plot_mpeg, xx, xv, yv, expand=expand, bytelimit=bytelimit, $
                 selfnormalize=selfnormalize, skip=skip, $
                 zlog=zlog, decades=decades, filename=filename, $
                 deleteimages=deleteimages, saveimages=saveimages, delay=delay, $
                     _EXTRA = e
;
; NAME:
;       TEST_IMAGE_MPEG
;
; PURPOSE:
;       Make an mpeg movie using image_plot features.
;
; CATEGORY:
;       General Graphics.
;
; CALLING SEQUENCE:
;       TEST_IMAGE_MPEG, xx
;
; INPUTS:
;       xx:    The three-dimensional array to make a movie out of.
;
; KEYWORD PARAMETERS:
;   EXPAND:        ?
;
;   BYTELIMIT:     ?
;
;   SELFNORMALIZE: If this keyword is set, each image will selfnormalize
;
;   SKIP:          The number of images you want to skip when creating the
;                  movie.
;
;   ZLOG:          Take the log of the amplitude (only 12 orders of
;                  magnitude in range allowed).
;
;   DECADES:       ? 
;
;   FILENAME:      The filename of the mpeg file, the defult is 'idl.mpg'
;
;   DELETEIMAGES:  If this keyword is set, the postscript images that
;                  are created to make the movie are deleted. (DEFAULT)
;
;   SAVEIMAGES:  If this keyword is set, the postscript images that
;                  are created to make the movie are saved. 
;
;   DELAY:         Puts x/100 seconds in between frames in the mpeg.
;
;   LEGEND:        If this keyword is set a colormap legend is created on
;                  the right of the image.
;
;   ZOOM:   Use only the interesting part of the data
;
; OUTPUTS:
;       Mpeg file in local directory and postscript files in
;       './mpeg_images'
;
; COMMON BLOCKS:
;	None.
;
; SIDE EFFECTS:
;       None.
;
; RESTRICTIONS:
;	None.

xx = reform(xx)
xxsize=size(xx)
if (xxsize(0) ne 3) then $
  message, 'The input to image_movie must be a 3D array'

if n_elements(filename) eq 0 then filename='idl.mpg'

if n_elements(skip) eq 0 then skip=1
;if n_elements(skip) eq 1 then skip=skip + 1
if n_elements(bytelimit) eq 0 then bytelimit=(1024.)^2*256.
if n_elements(expand) eq 0 then $
  if max([xxsize(1),xxsize(2)]) lt 128 then $
  expand=128/max([xxsize(1),xxsize(2)]) else expand=1

xsize=xxsize(1)*expand
ysize=xxsize(2)*expand
max_nimages=bytelimit/(xsize*ysize)
mend=xxsize(3)-1
mskip=skip
if xxsize(3)/mskip gt max_nimages then begin
    mskip=fix((mend-mstart)/max_nimages)+1 
    print, ' Plotting every ',mskip, ' images'
endif

if (keyword_set(zlog) and n_elements(decades) eq 0) then begin
  decades=4
  print, 'Movie using '+ strcompress(string(decades))+' decades of log(data)'
endif

if (n_elements(selfnormalize) eq 0) then begin
  xmin=min(xx)
  xmax=max(xx)
  if keyword_set(zlog) then begin
    xmin=alog(xmax*10.^(-1.*decades)) 
    xmax=alog(xmax) 
  endif
endif

;delete existing images in directory
existing_images = file_search('mpeg_images/*.ps')
n_existing_images =  n_elements(existing_images)
if n_existing_images gt 1 then FILE_DELETE, 'mpeg_images', /RECURSIVE

;make directory to store images in
FILE_MKDIR, 'mpeg_images'

range=[min(xx),max(xx)]

;write postscript images
for i = 0, (xxsize(3)-1), mskip do begin
    ps, 'mpeg_images/image' + string(FORMAT='(I05)', i) + '.ps',/color
    image_plot, xx(*, *, i), xv, yv, /REVERSE, ZLOG=zlog, DECADES=decades, $
      ZRANGE=range, SELFNORMALIZE=selfnormalize, _EXTRA = e
    stop_ps
endfor

;create movie
if (keyword_set(delay) or keyword_set(zlog)) then begin
    if not keyword_set(delay) then delay = 10
    spawn, 'convert -rotate 270 -delay ' + strcompress(delay, /REMOVE_ALL) + $
      ' mpeg_images/image*.ps ' + filename
endif else begin
    spawn, 'convert  -rotate 270 mpeg_images/image*.ps ' + filename
endelse

;delete image files
if keyword_set(images) then begin
    for i = 0, xxsize(3)-1 do begin
        FILE_DELETE, 'mpeg_images/image' + string(FORMAT='(I05)', i) + $
        '.ps'
    endfor
    FILE_DELETE, 'mpeg_images'
endif
 
end
