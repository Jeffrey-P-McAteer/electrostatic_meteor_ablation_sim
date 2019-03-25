pro image_movie, xx, expand=expand, bytelimit=bytelimit, $
                 selfnormalize=selfnormalize, skip=skip, $
                 zlog=zlog, decades=decades, zrange=zrange
; This program takes xx (a 3D array and turns it into a movie).

;on_error,2                      ;Return to caller if an error occurs

xxsize=size(reform(xx))
if (xxsize(0) ne 3) then $
  message, 'The input to image_movie must be a 3D array'

if n_elements(skip) eq 0 then skip=1
if n_elements(bytelimit) eq 0 then bytelimit=1024.^2*128.
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

if (n_elements(zrange) eq 2) then begin
	zmin=zrange(0)
	zmax=zrange(1)
endif else if (n_elements(selfnormalize) eq 0) then begin
  zmin=min(xx)
  zmax=max(xx)
  if keyword_set(zlog) then begin
    zmin=alog(zmax*10.^(-1.*decades)) 
    zmax=alog(zmax) 
  endif
endif

;-- first kill old processes
if xregistered('XInterAnimate',/noshow) ne 0 then xinteranimate,/close


xinteranimate, set=[xsize,ysize,(mend-mstart)/mskip+1],/showload,/track
for i=mstart,mend,mskip do begin
    frame=(reform(xx))[*,*,i]
    if keyword_set(zlog) then begin
        if (n_elements(selfnormalize) ne 0) then zmin=max(frame)*10.^(-1.*decades)
        frame=alog(frame>exp(zmin)) 
    endif
    if (n_elements(selfnormalize) eq 0) then $
      d2=bytscl(rebin(frame,xsize,ysize) $
                ,min=zmin,max=zmax,top=!d.table_size) $
    else $
      d2=bytscl(rebin(frame,xsize,ysize),top=!d.table_size)
        
    xinteranimate, frame=i/mskip, image=d2 
endfor
xinteranimate, /keep_pixmaps

; To close type:
;xinteranimate, /close

end



