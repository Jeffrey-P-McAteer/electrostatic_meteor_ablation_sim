;This is a routine to calculate the total electron flux in a 2-stream
;simulation.
;Lars Dyrud
;BU 3/25/02

@flux3d

; Read in and calculate the electro-static energy
@phi3d
!P.font=0
!P.charsize=2.7
!y.margin=[.5,1.5]
!x.margin=[11,11]

tmax=4000  ; maximum time to run the plot out to in plasma periods
if (n_elements(E_in_k_space) eq 0) then begin
    phisize=size(phi)

    E2avg=fltarr(phisize[4])
    E2max=fltarr(phisize[4])
    for i=0,phisize[4]-1 do begin
        E=gradient(phi(*,*,*,i),dx=dx*nout_avg,dy=dy*nout_avg,dz=dz*nout_avg)
        if (nz2 gt 1) then begin
            E2avg(i)=mean(E(*,*,*,0)^2+E(*,*,*,1)^2+E(*,*,*,2)^2)
            E2max(i)=max(E(*,*,*,0)^2+E(*,*,*,1)^2+E(*,*,*,2)^2)
        endif else begin
            if (ny2 gt 1) then begin
                E2avg(i)=mean(E(*,*,0)^2+E(*,*,1)^2)
                E2max(i)=max(E(*,*,0)^2+E(*,*,1)^2)
            endif else begin
                E2avg(i)=mean(E^2)
                E2max(i)=max(E^2)
            endelse
        endelse
    endfor

endif else message, 'phi2E2_3D not implemented in k space'



ionflux=reform(fluxx2)
flux=reform(fluxx0+fluxx1)
fluxsize=size(flux)


fluxtot =total(flux,1)
fluxtotsq=total(flux^2,1)
ionfluxtot=total(ionflux,1)

if fluxsize[0] eq 4 then begin 
    fluxtot =total(fluxtot,1)
    fluxtotsq=total(fluxtotsq,1)
    fluxtot =total(fluxtot,1)
    fluxtotsq=total(fluxtotsq,1)
    
    ionfluxtot =total(ionfluxtot,1)
    
    ionfluxtot =total(ionfluxtot,1)

    tim=dt*nout*VEL_OUT_SUBCYCLE0*findgen(fluxsize[4])
endif

if fluxsize[0] eq 3 then begin 
    fluxtot =total(fluxtot,1)
    fluxtotsq=total(fluxtotsq,1)
    
    ionfluxtot =total(ionfluxtot,1)
    tim=dt*nout*VEL_OUT_SUBCYCLE0*findgen(fluxsize[3])
endif else begin	

    tim=dt*nout*VEL_OUT_SUBCYCLE0*findgen(fluxsize[2])
endelse

;calculate Efield from momentum change
Efield=gradient(fluxtot,dx=dt*nout)/(fluxsize[1]*fluxsize[2])*5.68e-12*100.0*(1e5^2) 
Efield[0]=0
Efield[fluxsize[3]-1]=0


;calculate the effective collision frequency
coll_freq=Efield
coll_freq=gradient(fluxtot,dx=dt*nout)/(fluxtot)
coll_freq[0]=0
coll_freq[fluxsize[3]-1]=0

!p.multi=[0,1,3]

flux_norm= fluxtot/fluxtot(0)
plot, tim, flux_norm , ytitle='!4j!De!N',yrange=[min(flux_norm),max(flux_norm)],title='',xtickname=replicate(' ',30), xrange=[0,tmax]
;xyouts, tim[fluxsize[2]*0.7], flux[fluxsize[2]*0.7], 'Electron Beam'

plot, tim, (fluxtotsq/(fluxsize[1]*fluxsize[2])), ytitle='Avg. Energy',title='',xtickname=replicate(' ',30), xrange=[0,tmax]
oplot, tv, e2avg/20., linestyle=2
oplot,tim ((ionfluxtot/(fluxsize[1]*fluxsize[2]))^2)*md2, linestyle=3

plot, tim, Efield, ytitle='Est. E-Field (V/m)', yrange=[min(Efield), max(Efield)],ystyle=4,xtitle='Time (!9w!3!E-1!N!Dpe!N)', yticks=5, xrange=[0,tmax]
axis,yaxis=0, yrange=[min(Efield), max(Efield)],  ytitle='Est. E (V/m)', yticks=5
axis,yaxis=1, yrange=[min(coll_freq), max(coll_freq)],ytitle='Eff. Coll. Rate per (!9w!3!E-1!N!Dpe!N)',yticks=5

;plot, tim, coll_freq,xtitle='Time (inverse plasma periods)', ytitle='Eff. Coll. Rate per inMp)'
stop
end


