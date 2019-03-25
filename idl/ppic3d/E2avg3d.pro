
;Calculate the average E^2 and plot it

;Calculate E2 in real space:
if n_elements(phi) le 2 then begin
@phiep
endif
if n_elements(Bx) eq 0 then Bx=0
if n_elements(Bz) eq 0 then Bz=0

phisize=size(phi)
E2avg=fltarr(phisize[4])
E2max=fltarr(phisize[4])
; Each component may be useful:
Ex2avg=fltarr(phisize[4])
Ey2avg=fltarr(phisize[4])
Ez2avg=fltarr(phisize[4])

if (n_elements(boundary_type) eq 0) then boundary_type=0

for i=0L,phisize[4]-1 do begin
    E=gradient(phi[*,*,*,i],dx=dx*nout_avg,dy=dy*nout_avg,dz=dz*nout_avg)
    if (boundary_type eq 1) then E[0,*,*,*]=0
    if (boundary_type eq 1) then E[phisize[1]-1,*,*,*]=0
    if (nz2 gt 1) then begin

       ;Add a filtering capability to filter E along B
       if (fwidth gt 0 and fsteep gt 0) then begin
          E(*,*,*,0)=filter3d(E(*,*,*,0),fwidth=fwidth,fsteep=fsteep)
       endif

        Ex2avg(i)=mean(E(*,*,*,0)^2)
        Ey2avg(i)=mean(E(*,*,*,1)^2)
        Ez2avg(i)=mean(E(*,*,*,2)^2)


        E2avg(i)=Ex2avg(i)+Ey2avg(i)+Ez2avg(i)
        E2max(i)=max(E(*,*,*,0)^2+E(*,*,*,1)^2+E(*,*,*,2)^2)
    endif else begin
        if (ny2 gt 1) then begin
            Ex2avg(i)=mean(E(*,*,0)^2)
            Ey2avg(i)=mean(E(*,*,1)^2)
            E2avg(i)=Ex2avg(i)+Ey2avg(i)
            E2max(i)=max(E(*,*,0)^2+E(*,*,1)^2)
        endif else begin
            E2avg(i)=mean(E^2)
            E2max(i)=max(E^2)
        endelse
    endelse
endfor

plot,tv,E2max,/ylog,linestyle=1, yrange=[max([min(E2avg[1:*]),max(E2max)*1E-4]) ,max(E2max)], $
  title='E!E2!N Max and Avg', $
  xtitle='time (s)', $
  ytitle='E!E2!N (SI units)'
oplot,tv,E2avg
date_plot,title

;plot,tv,Ez2avg/E2avg,  title='Ez!E2!N/E!E2!N Avg',  xtitle='time
;(s)', ytitle='Ez!E2!N/E!E2!N',/ylog,yrange=[.001,1]*max(Ez2avg/E2avg)

plot,tv,Ex2avg/E2avg,  title='Ex!E2!N/E!E2!N Avg',  xtitle='time (s)', ytitle='Ex!E2!N/E!E2!N',/ylog,yrange=[.001,1]*max(Ex2avg/E2avg)

end


