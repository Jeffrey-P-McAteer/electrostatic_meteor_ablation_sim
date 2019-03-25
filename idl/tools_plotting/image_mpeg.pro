pro image_mpeg, xx, expand=expand, bytelimit=bytelimit, $
                 selfnormalize=selfnormalize, skip=skip, $
                 zlog=zlog, decades=decades, filename=filename, _EXTRA=extra_args
; This program takes xx (a 3D array and turns it into a movie).

; We need the color tables because the mpeg routines don't handle color
common colors, r_orig, g_orig, b_orig, r_curr, g_curr, b_curr
xx=reform(xx)
xxsize=size(xx)
if (xxsize(0) ne 3) then $
  message, 'The input to image_movie must be a 3D array'

if n_elements(filename) eq 0 then filename='idl.mpg'

if n_elements(skip) eq 0 then skip=1
if n_elements(bytelimit) eq 0 then bytelimit=(1024.)^3*4
if n_elements(expand) eq 0 then $
  if max([xxsize(1),xxsize(2)]) lt 128 then $
  expand=128/max([xxsize(1),xxsize(2)]) else expand=1

mstart=0
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

;xinteranimate, set=[xsize,ysize,(mend-mstart)/mskip+1],
;mpeg_filename=filename, /mpeg_open
mpegID=mpeg_open([xsize,ysize], _EXTRA=extra_args)
d3=bytarr(3,xsize,ysize)

for i=mstart,mend,mskip do begin
    frame=xx(*,*,i)
    if keyword_set(zlog) then frame=alog(frame>exp(xmin)) 
    if (n_elements(selfnormalize) eq 0) then $
      d2=bytscl(rebin(frame,xsize,ysize) $
                ,min=xmin,max=xmax,top=!d.table_size) $
    else $
      d2=bytscl(rebin(frame,xsize,ysize),top=!d.table_size)

;I shouldn't have to manually set the color tables like this but ...
    d3(0,*,*)=r_curr(d2)
    d3(1,*,*)=g_curr(d2)
    d3(2,*,*)=b_curr(d2)

;    xinteranimate, frame=i/mskip, image=d2 
    mpeg_put, mpegID, frame=i/mskip, image=d3, /order
endfor
;xinteranimate, /mpeg_close, /keep_pixmaps
mpeg_save,mpegID, filename=filename
mpeg_close,mpegID

; To close type:
;xinteranimate, /close

end



