; Analyze the output of a ppic3d run, making a set of standardized
; outputs

color_on=1
if n_elements(color_on) eq 0 then color_on=1 

;Plot out mean and max of E^2 averaged over the simulation box:
ps,'E2avg.ps'
.run E2avg3d
stop_ps

;Plot a sequence (or all) of the electric field energy (E^2) values:
if nx2 ge 2*ny2 then !p.multi=[0,2,4] else !p.multi=[0,4,4]
if ny2 eq 1 then !p.multi=0
if nx2 ge 2*ny2 then nplot=8 else nplot=16
iskip=round((phisize(4)-1)/(nplot))+1
if color_on eq 0 then ps,'E2_image-bw.ps',/landscape $
else ps,'E2_image-c.ps',/landscape,/color
.run E2image3d
stop_ps

iskip=round(phisize(4)/(4*8))+1
if ny2 gt 1 then !p.multi=[0,3,3]
avg_3d=0
ps,'E2kt-3x3-3d-bw.ps',/landscape
.run E2kt3d.pro
stop_ps

if ny2 gt 1 then !p.multi=[0,3,3]
avg_3d=1
ps,'E2kt-avg-3x3-3d-bw.ps',/landscape
.run E2kt3d.pro
stop_ps

;Look at spacial distribution of temps 
if color_on eq 0 then ps,'vsqr3d-bw.ps',/landscape $
else ps,'vsqr3d-c.ps',/landscape,/color
.r vsqr3d
stop_ps

ps
.r e2kwt
stop_ps

!p.multi=[0,2,2] 
veldist3d, 'ps'

@analyze_den3d

loadct,5
delvar,E
iskip=1
.run phi2E3d
if Esize(0) gt 1 then Eimage=[[reform(Ex(*,*,0,*))],[reform(Ey(*,*,0,*))]]
if Esize(0) gt 1 then image_mpeg,Eimage,/zlog,decades=3,expand=2,filename='ExEy.mpg'
delvar,Ex,Ey,Ez




