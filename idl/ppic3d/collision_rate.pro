;pro collision_rate.pro
;This program will read in the neccesary data files and produce plots
;of collision rates assuming a given field, Ey0_external
;
;This program is designed to work with the ppic3d data set and the
;drifts3d program
;
;Written by Meers Oppenheim Oct, 2002

if n_elements(vxavg0) eq 0 then message,' Run drift3d.pro first'

; Use Pederson drifts to calculate collision rates:
nu0_mu0=double(qd0)*Ey0_external/(md0*vyavg0)
nu0_ped0=(nu0_mu0-sqrt(nu0_mu0^2-4*wce^2) )/2
nu0_ped1=(nu0_mu0+sqrt(nu0_mu0^2-4*wce^2) )/2
;nu0_ped0a=wce^2*md0*vyavg0/(qd0*Ey0_external)

; Use Hall drifts to calculate collision rates:
nu0_h0=double(qd0)*Ey0_external*wce/(md0*vxavg0)
nu0_hall=sqrt(nu0_h0-wce^2)

nu1_mu1=double(qd1)*Ey0_external/(md1*vyavg1)
nu1_ped0=(nu1_mu1+sqrt(nu1_mu1^2-4*wci^2) )/2
nu1_ped1=(nu1_mu1-sqrt(nu1_mu1^2-4*wci^2) )/2


plot, time(5:*), nu0_ped0(5:*), yrange=[0,max(nu0_ped0(5:*))], $
  title = 'Dist 0 effective coll_rates', $
  xtitle='time (s)',ytitle='collision rates (1/s)'
oplot, time(5:*), nu0_ped1(5:*)
oplot, time(5:*), nu0_hall(5:*)

date_plot,title

plot, time(5:*), nu1_ped0(5:*), title = 'Dist 1 effective coll_rates', $
	xtitle='time (s)',ytitle='collision rates (1/s)'
oplot, time(5:*), nu1_ped1(5:*)
date_plot,title



end
