;pro drift3d
;This program will read in the neccesary data files and produce plots
;of Hall, Pederson drifts for both ions and electrons.
;
;This program is designed to work with the ppic3d data set.
;Written by Lars Dyrud Sep 2 2001
;
;Modified by Meers Oppenheim Nov 20,2001

; Read in Fluxes & input parameters, if necessary
@flux3d
; Read in density data if necessary
@den3d
; If the coll rates are 0 set them to 0
if n_elements(coll_rate0) eq 0 then coll_rate0=0
if n_elements(coll_rate1) eq 0 then coll_rate1=0
coll_rate0=coll_rate0*1.0
coll_rate1=coll_rate1*1.0

; Find out the size of the arrays
flux_size=size(fluxx0)
nt=min([flux_size(4),(size(den0))(4)])

;Theoretical values:
; Cyclotron Frequency:
	wce  = qd0*Bz/md0
	wci  = qd1*Bz/md1

; Parallel Mobilities:
	if (coll_rate0 ne 0) then mu0=qd0/coll_rate0/md0 else mu0 = 0
	if (coll_rate1 ne 0) then mu1=qd1/coll_rate1/md1 else mu0 = 0

; Pedersen mobility:
      ped0 = mu0/(1+wce^2/coll_rate0^2)
      ped1 = mu1/(1+wci^2/coll_rate1^2)

; Hall mobility:
      hall0 = qd0*wce/(md0*(wce^2+coll_rate0^2))
      hall1 = qd1*wci/(md1*(wci^2+coll_rate1^2))

; Calculate the expected drift rates
      v_hall0 = Ey0_external*hall0
      v_hall1 = Ey0_external*hall1

      v_ped0= Ey0_external*ped0
      v_ped1= Ey0_external*ped1

;declare arrays
vxavg0 = fltarr(nt)
vyavg0 = fltarr(nt)
vzavg0 =fltarr(nt)
vxavg1 = fltarr(nt)
vyavg1 = fltarr(nt)
vzavg1 = fltarr(nt)

; Loop through the data and calculate average drift rates as a
; function of time

for i=0, nt-1 do begin

      vxavg0(i)= mean(fluxx0(*,*,*,i)/((1+den0(*,*,*,i))*n0d0))
      vyavg0(i)= mean(fluxy0(*,*,*,i)/((1+den0(*,*,*,i))*n0d0))

      vxavg1(i)= mean(fluxx1(*,*,*,i)/((1+den1(*,*,*,i))*n0d1))
      vyavg1(i)= mean(fluxy1(*,*,*,i)/((1+den1(*,*,*,i))*n0d1))

      if n_elements(fluxz0) gt 2 then begin
          vzavg0(i)= mean(fluxz0(*,*,*,i)/((1+den0(*,*,*,i))*n0d0))
          vzavg1(i)= mean(fluxz1(*,*,*,i)/((1+den1(*,*,*,i))*n0d1))
      endif

endfor

time=dt*nout*findgen(nt)*iskip

if (!D.NAME eq 'X')  then window, 0, retain=2
vrange0=[min([min(vxavg0),min(vyavg0),min(vzavg0),v_hall0,v_ped0]), $
	 max([max(vxavg0),max(vyavg0),max(vzavg0),v_hall0,v_ped0])]
plot, time, vxavg0, title = 'Dist 0 velocities', yrange=vrange0, $
	xtitle='time (s)',ytitle='velocity (m/s)'
oplot, time, vyavg0
oplot, time, vzavg0
oplot, [time(0),time(nt-1)],[v_hall0,v_hall0]
xyouts,time(nt-1),v_hall0,' V hall 0 ',alignment=1.0
oplot, [time(0),time(nt-1)],[v_ped0,v_ped0]
xyouts,time(0),v_ped0,' V perdersen 0 ',alignment=0.0
date_plot,title

if (!D.NAME eq 'X')  then window, 1, retain=2
vrange1=[min([min(vxavg1),min(vyavg1),v_hall1,v_ped1]),$
	 max([max(vxavg1),max(vyavg1),v_hall1,v_ped1])]
plot, time, vxavg1, title = 'Dist 1 velocities', yrange=vrange1, $
	xtitle='time (s)',ytitle='velocity (m/s)'
oplot, time, vyavg1
oplot, [time(0),time(nt-1)],[v_hall1,v_hall1]
xyouts,time(nt-1),v_hall1,' V hall 1 ',alignment=1.0
oplot, [time(0),time(nt-1)],[v_ped1,v_ped1]
xyouts,time(0),v_ped1,' V perdersen 1 ',alignment=0.0
date_plot,title


end
