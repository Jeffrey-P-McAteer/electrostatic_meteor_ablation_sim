pro veldist3d, plottype = plottype, color_on=color_on
;Generate images of the velocity distributions

@vdist3d

if n_elements(plottype) eq 0 then plottype='X'
if plottype eq 'ps' then begin
    if color_on eq 0 then ps,'vdist0-bw.ps',/landscape $
    else ps,'vdist0-c.ps',/landscape,/color
endif

nvdist0=(size(vdist0))[4] 
if n_elements(vdist_out_subcycle) eq 0 then vdist_out_subcycle=1
dvdist=vdist_out_subcycle*nout*dt

;1D case:
if (size(vdist0))[2] eq 1 then begin
    image_plot,reform(vdist0), $
      vdistvx0,dvdist*findgen(nvdist0), $
      /zlog,decades=3,/legend,$
      title='Vel. Dist. 0', xtitle='Vx',ytitle='time' 
  date_plot,titlew
endif else begin
; 2D & 3D case
    for i=0.,nvdist0-1.,(nvdist0-1)/3. do begin
        image_plot,reform(vdist0(*,*,(size(vdist0))[3]/2,fix(i))),$
          vdistvx0,vdistvy0,/zlog,decades=3,/legend,$
          title='Vel. Dist. 0 at t='+string(i*dvdist),xtitle='Vx',ytitle='Vy' 
        date_plot,titlew
    endfor
    if plottype eq 'ps' then stop_ps
endelse

if (n_elements(vdist1) ge 1) then begin
    if plottype eq 'ps' then begin
        if color_on eq 0 then ps,'vdist1-bw.ps',/landscape $
        else ps,'vdist1-c.ps',/landscape,/color 
    endif
    nvdist1=(size(vdist1))[4]  

                                ;1D case:
    if (size(vdist1))[2] eq 1 then begin
        image_plot,reform(vdist1), $
          vdistvx1,dvdist*findgen(nvdist0), $
          /zlog,decades=3,/legend,$
          title='Vel. Dist. 1', xtitle='Vx',ytitle='time' 
        date_plot,titlew
    endif else begin
        for i=0.,nvdist1-1.,(nvdist1-1.)/3. do begin
            image_plot,reform(vdist1(*,*,(size(vdist1))[3]/2,fix(i))), $
              vdistvx1,vdistvy1,/zlog,decades=3,/legend, $
              title='Vel. Dist. 1 at t='+string(i*dvdist),xtitle='Vx',ytitle='Vy' 
            date_plot,title
        endfor
    endelse
    if plottype eq 'ps' then stop_ps
endif


if (n_elements(vdist2) ge 2) then begin
    if plottype eq 'ps' then begin
        if color_on eq 0 then ps,'vdist2-bw.ps',/landscape $
        else ps,'vdist2-c.ps',/landscape,/color 
    endif
    nvdist2=(size(vdist2))[4]  
                                ;1D case:
    if (size(vdist2))[2] eq 1 then begin
        image_plot,reform(vdist2), $
          vdistvx2,dvdist*findgen(nvdist2), $
          /zlog,decades=3,/legend,$
          title='Vel. Dist. 2', xtitle='Vx',ytitle='time' 
        date_plot,titlew
    endif else begin
        for i=0.,nvdist2-1.,(nvdist2-1.)/3. do begin
            image_plot,reform(vdist2(*,*,(size(vdist2))[3]/2,fix(i))),$
              vdistvx2,vdistvy2,/zlog,decades=3,/legend, $
              title='Vel. Dist. 2 at t='+string(i*dvdist),xtitle='Vx',ytitle='Vy'
            date_plot,title
        endfor
    endelse
    if plottype eq 'ps' then stop_ps
endif

end
