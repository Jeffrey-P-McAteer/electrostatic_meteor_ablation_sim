; Analyze the output of a ppic3d run, assuming no B0 field
; making a set of standardized outputs

color_on=1
if n_elements(color_on) eq 0 then color_on=1 

;Calculate envelope amplitude
.run envelope3d

; Plot an x vs t plot for the largest y,z crossection
!p.multi=0
ps,'A2iymax.ps',/color
.run A2iymax3d
stop_ps

;Plot a sequence (or all) of the electric field envelope energy (A^2) values:

phisize=size(phi)
iskip=round(phisize(4)/(4*8))+1
ps,'A2image-bw.ps',/landscape
!p.multi=[0,3,3]
.run A2image3d
stop_ps

; Make A2 movie:

loadct,5
iskip=1
if (size(reform(A2)))(0) gt 2 then $
  image_mpeg,reform(A2(*,*,0,*)),/zlog,decades=3,expand=2,filename='A2.mpg'

exit

