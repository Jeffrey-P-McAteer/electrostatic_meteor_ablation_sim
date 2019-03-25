;pro drift_test
;This program will read in the neccesary data files and produce plots
;of Hall, Pederson drifts for both ions and electrons.
;
;This program is designed to work with the ppic3d data set.
;Written by Lars Dyrud Sep 2 2001
;
;Modified by Meers Oppenheim Nov 20,2001


;First check to make sure the neccesary data files have been read in

if ( n_elements(velx0) le 2)  then message, 'Please run vel3d first'
@den3d

if n_elements(coll_rate0) eq 0 then coll_rate0=0
if n_elements(coll_rate1) eq 0 then coll_rate1=0
coll_rate0=coll_rate0*1.0
coll_rate1=coll_rate1*1.0

; Find out the size of the arrays
vel_size=size(velx0)
nt=vel_size(4)
nx=vel_size(1)
nz=vel_size(2)


;Theoretical values:
; Cyclotron Frequency:
	wce  = qd0*Bz/md0
	wci  = qd1*Bz/md1

; Parallel Mobilities:
	if (coll_rate0 ne 0) then mu0=qd0/coll_rate0/md0 else mu0 = 0
	if (coll_rate1 ne 0) then mu1=qd1/coll_rate1/md1 else mu0 = 0

; Pedersen mobility:
      ped0 = mu0/(1+wce^2/coll_rate0^2)
      ped1 = mu0/(1+wce^2/coll_rate0^2)

; Hall mobility:
      hall0 = qd0*wce/(md0*(wce^2+coll_rate0^2))
      hall1 = qd1*wce/(md1*(wce^2+coll_rate1^2))

; Calculate the expected drift rates
      v_hall0 = Ermag*hall0
      v_hall1 = Ermag*hall1

      v_ped0= Ermag*ped0
      v_ped1= Ermag*ped1

;declare arrays
vxavg0 = fltarr(nt)
vyavg0 = fltarr(nt)
vxavg1 = fltarr(nt)
vyavg1 = fltarr(nt)

; Loop through the data and calculate average drift rates as a
; function of time

for i=0, nt-1 do begin

      vxavg0(i)= mean(velx0(*,*,0,i)/((1+den0(*,*,0,i))*n0d0))
      vyavg0(i)= mean(vely0(*,*,0,i)/((1+den0(*,*,0,i))*n0d0))

      vxavg1(i)= mean(velx1(*,*,0,i)/((1+den1(*,*,0,i))*n0d1))
      vyavg1(i)= mean(vely1(*,*,0,i)/((1+den1(*,*,0,i))*n0d1))

endfor

time=dt*nout*findgen(nt)

vrange0=[min([min(vxavg0),min(vyavg0)]), max([max(vxavg0),max(vyavg0)])]
plot, time, vxavg0, title = 'Dist 0 velocities', yrange=vrange0
oplot, time, vyavg0
oplot, [time(0),time(nt-1)],[v_hall0,v_hall0]
oplot, [time(0),time(nt-1)],[v_ped0,v_ped0]

vrange1=[min([min(vxavg1),min(vyavg1)]), max([max(vxavg1),max(vyavg1)])]
plot, time, vxavg1, title = 'Dist 1 velocities', yrange=vrange1
oplot, time, vyavg1
oplot, [time(0),time(nt-1)],[v_hall1,v_hall1]
oplot, [time(0),time(nt-1)],[v_ped1,v_ped1]





end