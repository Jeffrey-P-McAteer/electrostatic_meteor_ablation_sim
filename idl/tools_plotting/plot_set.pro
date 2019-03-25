pro plot_set, a0, xv, yv, TITLE=title, dt_scale=dt_scale, label=label
;
; NAME:
;	PLOT_SET
;
; PURPOSE:
;	Make a set of summary plots of 2-D, 3-D or 4-D where the last
;	dimension is time and the previous dimensions are spatial.
;
; CATEGORY:
;	General graphics.
;
; CALLING SEQUENCE:
;	PLOT_SET, A0, xv, yv
;
; INPUTS:
;	A0:	The two-dimensional array to display.
;       xv:      A vector of x values
;       yv:      A vector of y values
;
; KEYWORD PARAMETERS:
;       title:    A string
;       dt_scale: A scaler to multiply by the last dimensions'
;                 counter to put on the lable
;
nA0=size(A0)
if n_elements(title) eq 0 then title=''
if n_elements(label) eq 0 then label=''
if n_elements(dt_scale) eq 0 then dt_scale=1


nx=nA0(1)
ny=nA0(2)
; Estimate the number to put on a page based on the ratio of nx to ny
!p.multi=0
if (nx ge 2*ny) then nplot=8 else nplot=16
if (nx ge 2*ny) then !p.multi=[0,2,4] else !p.multi=[0,4,4] 


for i=0.,nA0(4)-1.,(nA0(4)-1.)/(nplot-1.)*0.999999 do begin
   title2=title + ' ' + string(fix(i)*dt_scale)
   image_plot, reform(A0(*,*,(size(A0))[3]/2,fix(i))),$
               xv,yv,/legend,/aspect,$
               title=title2,xtitle='x',ytitle='y' 
   date_plot,label
endfor
if (nA0[3] gt 1) then begin
   for i=0.,nA0(4)-1.,(nA0(4)-1.)/(nplot-1.)*0.999999 do begin
      title2=title + ' ' + string(fix(i)*dt_scale)
      image_plot, reform(A0(*,(size(A0))[2]/2,*,fix(i))),$
                  xv,zv,/legend,/aspect,$
                  title=title2,xtitle='x',ytitle='z' 
      date_plot, label
   endfor
endif

if (nA0[3] gt 1) then begin
   for i=0.,nA0(4)-1.,(nA0(4)-1.)/(nplot-1.)*0.999999 do begin
      title2=title + ' ' + string(fix(i)*dt_scale)
      image_plot, transpose(reform(A0((size(A0))[1]/2,*,*,fix(i)))),$
                  zv,yv,/legend,/aspect,$
                  title=title2,xtitle='z',ytitle='y' 
      date_plot,label
   endfor
endif

!p.multi=1

return
end

