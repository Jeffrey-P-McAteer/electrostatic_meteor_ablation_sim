;Calculate the average |den| and plot it

if n_elements(den0) and n_elements(den1) le 2 then begin
@den3d
endif

;if den0 is zero, try to sub den 1 for den 0

if n_elements(den0) lt 2 then begin
    den0=den1
    den0swapped=1
    endif else den0swapped=0


den0size=size(den0)
den0avg=fltarr(den0size[4])
den0max=fltarr(den0size[4])
for i=0,den0size[4]-1 do begin
    den0avg(i)=sqrt(mean(den0(*,*,*,i)^2))
    den0max(i)=max(abs(den0(*,*,*,i)))
endfor

if n_elements(den1) gt 2 then begin
    den1size=size(den1)
    den1avg=fltarr(den1size[4])
    den1max=fltarr(den1size[4])
    for i=0,den1size[4]-1 do begin
        den1avg(i)=sqrt(mean(den1(*,*,*,i)^2))
        den1max(i)=max(abs(den1(*,*,*,i)))
    endfor
endif

if n_elements(den2) gt 2 then begin
    den2size=size(den2)
    den2avg=fltarr(den2size[4])
    den2max=fltarr(den2size[4])
    for i=1,den2size[4]-1 do begin
        den2avg(i)=sqrt(mean(den2(*,*,*,i)^2))
        den2max(i)=max(abs(den2(*,*,*,i)))
    endfor
endif

if n_elements(den3) gt 3 then begin
    den3size=size(den3)
    den3avg=fltarr(den3size[4])
    den3max=fltarr(den3size[4])
    for i=1,den3size[4]-1 do begin
        den3avg(i)=sqrt(mean(den3(*,*,*,i)^2))
        den3max(i)=max(abs(den3(*,*,*,i)))
    endfor
endif

;Make the plot
if n_elements(den0) gt 2 then nt=den0size(4) 

if n_elements(den_out_subcycle0) eq 0 then den_out_subcycle0=1
tv=findgen(nt)*dt*nout*iskip*den_out_subcycle0

; Determine yrange
yrange=[min(den0avg[2:nt-1]),max(den0avg[2:nt-1])] 

plot,tv[2:nt-1],den0avg[2:nt-1], yrange=yrange,$
  title='Perturbed Density Average', $
  xtitle='time (s)', ytitle='|den/n0| (SI units)'
;oplot,tv,den0avg
xyouts,tv(den0size(4)-1),den0avg(den0size(4)-1),' D0'

if (n_elements(den1) gt 2) then begin
;    oplot,tv,den1max,linestyle=1
    oplot,tv,den1avg,linestyle=1
    xyouts,tv(den1size(4)-1),den1avg(den1size(4)-1),' D1'
endif

if (n_elements(den2) gt 2) then  begin
;    oplot,tv,den2max,linestyle=2
    oplot,tv,den2avg,linestyle=2
    xyouts,tv(den2size(4)-1),den2avg(den2size(4)-1),' D2'
endif

if (n_elements(den3) gt 2) then  begin
;    oplot,tv,den3max,linestyle=3
    oplot,tv,den3avg,linestyle=3
    xyouts,tv(den3size(4)-1),den3avg(den3size(4)-1),' D3'
endif

if den0swapped then den0=0

date_plot,title
end


