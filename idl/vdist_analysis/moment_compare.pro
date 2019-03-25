; This script plots a number of temps on the same plot

cd,'/project/eregion/meerso/eppic-archive/FBI/'
dirs=['FBI3D512x512x1024psi.3n04dx.1E140nofilter','FBI3D2048E140nu.5V111210','FBI2D512x512psi.3n04dx.1E140-T2','FBI2D2048x2048E140nu.5']
labels=['3D:1.0','3D:0.5','2D:1.0','2D:0.5']

ndirs=n_elements(dirs)

id=0
cd,dirs[id]
@moment_plot_ep

; Make the arrays
Temp_e=fltarr(ndirs,mom_size(2)*8)
Temp_i=fltarr(ndirs,mom_size(2)*8)
time=fltarr(ndirs,mom_size(2)*8)
imax=intarr(ndirs)

;Fill arrays
Time[id,0:imomax]=moments0[0,*]*dt
Temp_e[id,0:imomax]=(moments0[2,*]+moments0[6,*]+moments0[10,*])/3.*md0/kb
Temp_i[id,0:imomax]=(moments1[2,*]+moments1[6,*]+moments1[10,*])/3.*md1/kb
imax(id)=imomax
cd,'..'

id=1
cd,dirs[id]
@moment_plot_ep
;Fill arrays
Time[id,0:imomax]=moments0[0,*]*dt
Temp_e[id,0:imomax]=(moments0[2,*]+moments0[6,*]+moments0[10,*])/3.*md0/kb
Temp_i[id,0:imomax]=(moments1[2,*]+moments1[6,*]+moments1[10,*])/3.*md1/kb
imax(id)=imomax
cd,'..'

id=2
cd,dirs[id]
@moment_plot_ep
;Fill arrays
Time[id,0:imomax]=moments0[0,*]*dt
Temp_e[id,0:imomax]=(moments0[2,*]+moments0[6,*]+moments0[10,*])/3.*md0/kb
Temp_i[id,0:imomax]=(moments1[2,*]+moments1[6,*]+moments1[10,*])/3.*md1/kb
imax(id)=imomax
cd,'..'

id=3
cd,dirs[id]
@moment_plot_ep
;Fill arrays
Time[id,0:imomax]=moments0[0,*]*dt
Temp_e[id,0:imomax]=(moments0[2,*]+moments0[6,*]+moments0[10,*])/3.*md0/kb
Temp_i[id,0:imomax]=(moments1[2,*]+moments1[6,*]+moments1[10,*])/3.*md1/kb
imax(id)=imomax
cd,'..'

;id=4
;cd,dirs[id]
;@moment_plot_ep
;;Fill arrays
;Time[id,0:imomax]=moments0[0,*]*dt
;Temp_e[id,0:imomax]=(moments0[2,*]+moments0[6,*]+moments0[10,*])/3.*md0/kb
;Temp_i[id,0:imomax]=(moments1[2,*]+moments1[6,*]+moments1[10,*])/3.*md1/kb
;cd,'..'

;Plot distribution 0 temps against each other:
; Find the range of temps
rangeT=[min(Temp_e),max(Temp_e)*1.05]
range_time=[min(Time),max(Time)]
id=0
ps,'eTE140.ps'
!p.charsize=1.5
plot, time[id,0:imax(id)-1], Temp_e[id,0:imax(id)-1],  $
  xstyle=1, xrange=range_time,  ystyle=1, yrange=rangeT, $
  title='Electron moments for E=140mV dif. coll. rates ', ytitle='T (K)', xtitle='time (s)' 
  xyouts, range_time[1],Temp_e[id,imax(id)-1], labels[id], alignment=1.0
for id=1, ndirs-1 do $
  oplot, time[id,0:imax(id)-1], Temp_e[id,0:imax(id)-1],linestyle=id
for id=1, ndirs-1 do xyouts, range_time[1],Temp_e[id,imax(id)-1], labels[id], alignment=1.0

date_plot,title
stop_ps


