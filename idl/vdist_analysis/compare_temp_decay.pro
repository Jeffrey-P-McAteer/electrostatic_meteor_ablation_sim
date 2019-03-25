; This script compares temperature decays from FB eppic to theoretical
; decay curves

!p.multi=[0,2,2]

; Start by gathering data from run
@moment_plot_ep.pro

!p.multi=[0,2,1]

t=reform(moments1[0,*]*dt)

max_range=min([n_elements(t)-1,32])
t=t[0:max_range-1]

Tn= massd_neutral0*(1d0*Vthd_neutral0)^2/kb
T0= md0*(1d0*vxthd0)^2/kb
;nu=coll_rate0/9.2
nu=coll_rate0/3.

Tdif=(T0^.5-Tn^.5)

;For nu proportional to T^(1/2)
Thalf=Tn*((T0^.5+Tn^.5+Tdif*exp(-(Tn/T0)^.5*nu*t))/ $
           (T0^.5+Tn^.5-Tdif*exp(-(Tn/T0)^.5*nu*t)))^2

;For nu proportional to T
Tone=Tn*T0/(T0*(1.-exp(-(Tn/T0)*nu*t))+Tn*(exp(-(Tn/T0)*nu*t)))

;For nu not proportional to T
Tconst=Tn+(T0-Tn)*exp(-nu*t)

plot, t, moments0[2,0:max_range-1]*md0/kb-Tn, title='e- temp evolution comparison',yrange=[0.001,1]*(T0-Tn), ytitle='Temp - Tn (K)'
oplot, t, Tone-Tn, color= 40
oplot, t, Tconst-Tn, color= 120
oplot, t, thalf-Tn , color=160
xyouts, t[max_range-1]*3./4.,max(moments0[2,*]*md0/kb-Tn)*.95, ' Simulation'
xyouts, t[max_range-1]*3./4.,max(moments0[2,*]*md0/kb-Tn)*.90, ' Thalf',color=40
xyouts, t[max_range-1]*3./4.,max(moments0[2,*]*md0/kb-Tn)*.85, ' Tone',color=120
xyouts, t[max_range-1]*3./4.,max(moments0[2,*]*md0/kb-Tn)*.80, ' Tconst',color=160

; Ions ...

nu=coll_rate1/0.9

Tn= m_neutral*(1d0*Vth_neutral)^2/kb
T0= md1*(1d0*vxthd1)^2/kb


Tdif=(T0^.5-Tn^.5)

;For nu proportional to T^(1/2)
Thalf=Tn*((T0^.5+Tn^.5+Tdif*exp(-(Tn/T0)^.5*nu*t))/ $
           (T0^.5+Tn^.5-Tdif*exp(-(Tn/T0)^.5*nu*t)))^2

;For nu proportional to T
Tone=Tn*T0/(T0*(1.-exp(-(Tn/T0)*nu*t))+Tn*(exp(-(Tn/T0)*nu*t)))

;For nu not proportional to T
Tconst=Tn+(T0-Tn)*exp(-nu*t)

plot, t, moments1[2,0:max_range-1]*md1/kb-Tn, title='Ion Temp evolution comparison', ytitle='Temp - Tn (K)',ystyle=1
oplot, t, Tone-Tn, color= 40
oplot, t, Tconst-Tn, color= 120
oplot, t, thalf-Tn, color=160
xyouts, t[max_range-1]*3./4.,max(moments1[2,*]*md0/kb-Tn)*.95, ' Simulation'
xyouts, t[max_range-1]*3./4.,max(moments1[2,*]*md0/kb-Tn)*.90, ' Thalf',color=40
xyouts, t[max_range-1]*3./4.,max(moments1[2,*]*md0/kb-Tn)*.85, ' Tone',color=120
xyouts, t[max_range-1]*3./4.,max(moments1[2,*]*md0/kb-Tn)*.80, ' Tconst',color=160
date_plot
