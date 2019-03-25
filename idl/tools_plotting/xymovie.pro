pro xymovie, xin, yin, wait = wait, xrange=xrange, yrange=yrange, $
             selfnormalize = selfnormalize, $
             color=color, crange=crange, $
             skip=skip, label=label, title=title0, _EXTRA = e
; Animate a set of line data

on_error,2

if (n_elements(yin) eq 0) then begin
    y=xin
    ysize=size(y)
    x=findgen(ysize(1))
endif else begin
    y=yin
    x=xin
endelse

ysize=size(y)
if (ysize(0) ne 2) then $
  message, 'The input to image_movie must be a 2D array'

xsize=size(x)
if (ysize(1) ne xsize(1)) then $
  message, 'The x array must match the y array in size'

if n_elements(wait) eq 0 then wait = 0.1

if n_elements(yrange) ne 2 then begin
    ymin=min(y)
    ymax=max(y)
endif else begin
    ymin=yrange(0)
    ymax=yrange(1)
endelse

if n_elements(xrange) ne 2 then begin
    xmin=min(x)
    xmax=max(x)
endif else begin
    xmin=xrange(0)
    xmax=xrange(1)
endelse

if ( n_elements(skip) eq 0) then skip=1
if ( n_elements(title) eq 0) then title=''

if ( n_elements(color) eq 0 ) then color=0 
if ( n_elements(color) eq 1 ) then begin
    for i=0, ysize(2)-1, skip do begin 
        if keyword_set(label) then title=title0 + ' ' + strcompress(string(i))
        if n_elements(selfnormalize) eq 1 then begin
            ymin=min(y[*,i])
            ymax=max(y[*,i])
        endif
        if (xsize(0) eq 1) then $
          plot, x, y(*,i), xrange=[xmin,xmax], yrange=[ymin,ymax], xstyle=1, title=title,  _EXTRA = e $
        else $
          plot, x(*,i), y(*,i), xrange=[xmin,xmax], yrange=[ymin,ymax], xstyle=1,  title=title, _EXTRA = e
        wait, wait
    endfor
endif else begin $
                                ;Colored movies
    if n_elements(crange) ne 2 then begin
        cmin=min(color)
        cmax=max(color)
    endif else begin
        cmin=crange(0)
        cmax=crange(1)
    endelse

    for i=0, ysize(2)-1, skip do begin 
        if keyword_set(label) then title=title0 + ' ' + string(i)
                                ; Define Axis
        plot, [xmin,xmax], [ymin,ymax], xrange=[xmin,xmax], yrange=[ymin,ymax], $
          xstyle=1, title=title, /nodata
                                ; Designate Colors
        if (n_elements(color) eq n_elements(y)) then $
            colobj=fix(((((color(*,i)-cmin)/(cmax-cmin)) > 0.) < 1.)*!D.N_COLORS) $
            else colobj=fix(((((color-cmin)/(cmax-cmin)) > 0.) < 1.)*!D.N_COLORS) 
            
        if (xsize(0) eq 1) then $
          plots, x, y(*,i), color=colobj, _EXTRA = e $
        else $
          plots, x(*,i), y(*,i), color=colobj, _EXTRA = e
        wait, wait

        t_colbar, length=4, upper=!D.N_COLORS, xpos=!D.X_SIZE-12, $
          ypos=!D.Y_SIZE/2 - !D.N_COLORS*4/2., rotation=1

        wait, wait
    endfor
endelse


end
