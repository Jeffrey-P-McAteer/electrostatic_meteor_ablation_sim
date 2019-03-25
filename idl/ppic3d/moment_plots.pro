;pro moment_plots
;This program will read in the neccesary data files and produce plots
;of Hall, Pederson drifts for both ions and electrons
;from moment information.
;
;written by Meers Oppenheim May 2006

;ps,'param_evol.ps'
@ppic3d.i
@moments3d

if n_elements(Bz) eq 0 and n_elements(Bx) ne 0 then Bz=Bx


; We need the Boltzmann constant will usually be mks or 1
if (n_elements(kb) eq 0) then if (md0 lt 1e-8) then kb=1.38e-23 else kb=1

;Theoretical values:
; Cyclotron Frequency:
	wce  = qd0*Bz/md0
	wci  = qd1*Bz/md1

; Parallel Mobilities:
        if (n_elements(coll_rate0) eq 0) then coll_rate0=0.
        if (n_elements(coll_rate1) eq 0) then coll_rate1=0.
	if (coll_rate0 ne 0) then mu0=qd0/coll_rate0/md0 else mu0 = 0
	if (coll_rate1 ne 0) then mu1=qd1/coll_rate1/md1 else mu0 = 0

; Pedersen mobility:
      ped0 = mu0/(1+wce^2/(coll_rate0*1.0)^2)
      ped1 = mu1/(1+wci^2/(coll_rate1*1.0)^2)

; Hall mobility:
      hall0 = qd0*wce/(md0*(wce^2+coll_rate0^2))
      hall1 = qd1*wci/(md1*(wci^2+coll_rate1^2))

; Calculate the expected drift rates
      if (n_elements(Ey0_external) eq 0) then Ey0_external=0.
      v_hall0 = Ey0_external*hall0
      v_hall1 = Ey0_external*hall1

      v_ped0= Ey0_external*ped0
      v_ped1= Ey0_external*ped1


mom_size=size(moments0)
imomax=mom_size(2)-1

maxtemp0=max([moments0[2,*],moments0[6,*],moments0[10,*]],/nan)*md0/kb
mintemp0=min([moments0[2,*],moments0[6,*],moments0[10,*]],/nan)*md0/kb

maxtemp1=max([moments1[2,*],moments1[6,*],moments1[10,*]],/nan)*md1/kb
mintemp1=min([moments1[2,*],moments1[6,*],moments1[10,*]],/nan)*md1/kb

;Plot distribution 0 temps
plot, moments0[0,*]*dt, moments0[2,*]*md0/kb, $
  ystyle=1, yrange=[mintemp0,maxtemp0], $
  title='Temperatures 0: ', ytitle='T (K)', xtitle='time (s)'
xyouts,moments0[0,imomax]*dt,moments0[2,(imomax)]*md0/kb, ' x0'
oplot, moments1[0,*]*dt, moments0[6,*]*md0/kb
xyouts,moments0[0,imomax]*dt,moments0[6,(imomax)]*md0/kb, ' y0'
oplot, moments1[0,*]*dt, moments0[10,*]*md0/kb
xyouts,moments0[0,imomax]*dt,moments0[10,(imomax)]*md0/kb, ' z0'
date_plot,title

;Plot distribution 1 temps
plot, moments1[0,*]*dt, moments1[2,*]*md1/kb, $
  ystyle=1, yrange=[mintemp1,maxtemp1],$
  title='Temperatures 1: ', ytitle='T (K)', xtitle='time (s)'
xyouts,moments1[0,imomax]*dt,moments1[2,(imomax)]*md1/kb, ' x1'
oplot, moments1[0,*]*dt, moments1[6,*]*md1/kb
xyouts,moments1[0,imomax]*dt,moments1[6,(imomax)]*md1/kb, ' y1'
oplot, moments1[0,*]*dt, moments1[10,*]*md1/kb
xyouts,moments1[0,imomax]*dt,moments1[10,(imomax)]*md1/kb, ' z1'
maxvel0=max([moments0[1,*],moments0[5,*],moments0[9,*], $
           moments1[1,*],moments1[5,*],moments1[9,*]],/nan)
minvel0=min([moments0[1,*],moments0[5,*],moments0[9,*], $
           moments1[1,*],moments1[5,*],moments1[9,*]],/nan)
date_plot,title

t0tfin=[moments0[0,0],moments0[0,imomax]]*dt
;Plot distribution 0 velocities
plot, moments0[0,*]*dt, moments0[1,*], ystyle=1, yrange=[minvel0,maxvel0], $
  title='Velocities 0: ', ytitle='V (m/s)', xtitle='time (s)'
xyouts,moments0[0,imomax]*dt,moments0[1,(imomax)], ' Vx0'
oplot, moments1[0,*]*dt, moments0[5,*]
xyouts,moments0[0,imomax]*dt,moments0[5,(imomax)], ' Vy0'
oplot, moments1[0,*]*dt, moments0[9,*]
xyouts,moments0[0,imomax]*dt,moments0[9,(imomax)], ' Vz0'
;Draw lines with the expected velocities:
oplot,t0tfin,[v_hall0,v_hall0]
xyouts,moments0[0,imomax]*dt,v_hall0, ' V_hall0'
oplot,t0tfin,[v_ped0,v_ped0]
xyouts,moments0[0,imomax]*dt,v_ped0, ' V_ped0'
date_plot,title

;Plot distribution 1 velocities
maxvel1=max([moments1[1,*],moments1[5,*],moments1[9,*]],/nan)
minvel1=min([moments1[1,*],moments1[5,*],moments1[9,*]],/nan)
plot, moments1[0,*]*dt, moments1[1,*], ystyle=1, yrange=[minvel1,maxvel1], $
  title='Velocities 1: ', ytitle='V (m/s)', xtitle='time (s)'
xyouts,moments1[0,imomax]*dt,moments1[1,(imomax)], ' Vx1'
oplot, moments1[0,*]*dt, moments1[5,*]
xyouts,moments1[0,imomax]*dt,moments1[5,(imomax)], ' Vy1'
oplot, moments1[0,*]*dt, moments1[9,*]
xyouts,moments1[0,imomax]*dt,moments1[9,(imomax)], ' Vz1'
;Draw lines with the expected velocities:
oplot,t0tfin,[v_hall1,v_hall1]
xyouts,moments1[0,imomax]*dt,v_hall1, ' V_hall1'
oplot,t0tfin,[v_ped1,v_ped1]
xyouts,moments1[0,imomax]*dt,v_ped1, ' V_ped1'
date_plot,title

; Calculate some important parameters
Te=(moments0[2,*] + moments0[6,*] + moments0[10,*])/3. * md0/kb
Ti=(moments1[2,*] + moments1[6,*] + moments1[10,*])/3. * md1/kb

Cs=sqrt(kb*(Te+Ti)/md1)
Te_start=0.5*md0/kb*((vxthd0*1.0)^2+(vythd0^2*1.0)+(vxthd0*1.0)^2)
Ti_start=0.5*md1/kb*((vxthd1*1.0)^2+(vythd1^2*1.0)+(vxthd1*1.0)^2)

Cs_start=sqrt(kb*(Te_start + Ti_start)/md1)

; electron collision rate can be calculted from the EXB-direction drift speed
; or the pederson drift speed. The best method may be to use the ratio of
; these speeds, essentially using the drift angle:

nue=moments0[5,*]/moments0[1,*]*wce
nue2=wce*sqrt(Ey0_external/Bz/moments0[1,*] -1)
nui=Ey0_external/(moments1(5,*))*(qd1/md1)

Psi=abs(nue*nui/(wce*wci))

driver=moments0(1,*)/(1+Psi)

imomin=8
plot,moments1[0,imomin:imomax]*dt,Cs[imomin:imomax], $
ystyle=1, yrange=[Cs_start,max([Cs[imomin:imomax],driver[imomin:imomax]])*1.02], $
  title='Acoustic speed:', ytitle='Cs (m/s)', xtitle='time (s)'
xyouts,moments1[0,imomax]*dt,Cs[(imomax)], ' Cs'
oplot,t0tfin,[Cs_start,Cs_start]
xyouts,moments1[0,imomax]*dt,Cs_start, ' Cs_start'
oplot,moments1[0,imomin:imomax]*dt, Driver[imomin:imomax]
xyouts,moments1[0,imomax]*dt,Driver[(imomax)], ' Driver'
date_plot,title

; Collision rates
plot,moments1[0,imomin:imomax]*dt,nue[imomin:imomax], $
  ystyle=1, yrange=[min([nue[imomin:imomax],nui[imomin:imomax]])>0.,max([nue[imomin:imomax],nue[imomin:imomax]])*1.1], $
  title='collision rates:', ytitle='freq (1/s)', xtitle='time (s)'
xyouts,moments1[0,imomax]*dt,nue[(imomax)], ' nue'
oplot,moments1[0,imomin:imomax]*dt,nui[imomin:imomax]
xyouts,moments1[0,imomax]*dt,nui[(imomax)], ' nui'
date_plot,title

; Psi parameter
plot,moments1[0,imomin:imomax]*dt,Psi[imomin:imomax], $
  ystyle=1, yrange=[0,max(Psi[imomin:imomax])*1.1], $
  title='Psi:', ytitle='Psi', xtitle='time (s)'
date_plot,title


;stop_ps


