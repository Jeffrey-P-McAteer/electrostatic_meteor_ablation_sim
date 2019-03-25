;pro temp_predict.pro
;This program will read in the neccesary data files to predict temperatures
;given appropriate parameters.
;
; by Meers Oppenheim Nov 20, 2002

@params_in.pro

; If the coll rates are 0 set them to 0
if n_elements(coll_rate0) eq 0 then coll_rate0=0
if n_elements(coll_rate1) eq 0 then coll_rate1=0
; Fix the collision rates by hand:
coll_rate0=coll_rate0/1.56
coll_rate1=coll_rate1*106./175.

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
      v_hall0 = double(Ey0_external*hall0)
      v_hall1 = double(Ey0_external*hall1)

      v_ped0= double(Ey0_external*ped0)
      v_ped1= double(Ey0_external*ped1)

if (n_elements(kb) eq 0) then if (md0<1e-8) then kb=1.38e-23 else kb=1

; Calculate the expected temperatures
      if n_elements(massd_neutral0) eq 0 then massd_neutral0=m_neutral
      if n_elements(vthd_neutral0) eq 0 then vthd_neutral0=vth_neutral
      vthd_neutral0=double(vthd_neutral0)
      mu0=double(massd_neutral0)*md0/(md0+massd_neutral0)
      delta0=2*md0/(md0+massd_neutral0)
      T0=( (2/3.)*(v_hall0^2+v_ped0^2)*mu0/delta0 + vthd_neutral0^2*massd_neutral0 )/kb
      print, 'Expected equilibrium T0=',T0

      if n_elements(massd_neutral1) eq 0 then massd_neutral1=m_neutral
      if n_elements(vthd_neutral1) eq 0 then vthd_neutral1=vth_neutral
      vthd_neutral1=double(vthd_neutral1)
      mu1=double(massd_neutral1)*md1/(md1+massd_neutral1)
      delta1=2*md1/(md1+massd_neutral1)
      T1=( (2/3.)*(v_hall1^2+v_ped1^2)*mu1/delta1 + vthd_neutral1^2*massd_neutral1 )/kb
      print, 'Expected equilibrium T1=',T1

end
