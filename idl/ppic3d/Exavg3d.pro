;Calculate the average Ex and plot it

@phi3d

;Calculate Ex in real space:

phisize=size(phi)
Eymax=fltarr(phisize[4])
Exmax=fltarr(phisize[4])
for i=0,phisize[4]-1 do begin
    E=gradient(phi(*,*,*,i),dx=dx*nout_avg,dy=dy*nout_avg,dz=dz*nout_avg)
    if (nz2 gt 1) then begin
        Exmax(i)=max(E(*,*,*,0))
        Eymax(i)=max([E(*,*,*,1),E(*,*,*,2)])
    endif else begin
        if (ny2 gt 1) then begin
            Exmax(i)=max(E(*,*,0))
            Eymax(i)=max(E(*,*,1))
        endif else begin
            Exmax(i)=max(Ex)
            Eymax(i)=0
        endelse
    endelse
endfor


Emax=max([max(Exmax),max(Eymax)])
plot,tv,Exmax, linestyle=1, title='E!Ex!N max and E!Ey!N max', $
  yrange=[0,Emax], xtitle='time (s)', ytitle='E (SI units)'
oplot,tv,Eymax
date_plot,title
end


