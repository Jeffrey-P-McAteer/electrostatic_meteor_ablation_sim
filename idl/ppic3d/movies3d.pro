; Analyze the output of a ppic3d run, making a set of standardized
; outputs

if n_elements(color_on) eq 0 then color_on=1 

;Read in the potential
@phi3d
if phisize(0) eq 4 then image_mpeg,reform(phi[*,*,0,*]),expand=4,filename='phi.mpg' &
else if phisize(0) eq 3 then image_mpeg,reform(phi),expand=4,filename='phi.mpg'

loadct,5
delvar,E
iskip=1
.run phi2E3d
if phisize(0) gt 1 then Eimage=[[reform(Ex(*,*,0,*))],[reform(Ey(*,*,0,*))]]
if phisize(0) gt 1 then image_mpeg,Eimage,/zlog,decades=3,expand=2,filename='ExEy.mpg'
delvar,phi,Ex,Ey,Ez



