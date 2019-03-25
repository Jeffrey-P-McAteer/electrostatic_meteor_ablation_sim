pro altitude, Z, x, y , log=log, _EXTRA = e
;This creates a contour plot of any complex Matrix

Zs=size(Z)
if Zs(0) ne 2 then begin
	print, 'Error: Argument to altitude must be a complex 2D array'
	return
endif
if not (Zs(3) eq 6 or Zs(3) eq 9) then  begin
	print, 'Error: Argument to altitude must be a complex 2D array'
	return
endif

if (n_elements(x) eq 0) then x=findgen(Zs(1))
if (n_elements(y) eq 0) then y=findgen(Zs(1))
	
if n_elements(log) eq 0 then $
  contour, abs(Z), x, y, /FOLLOW, NLEVELS=16., _EXTRA = e $
else $
  contour, alog(abs(Z)), x, y, /FOLLOW, NLEVELS=16., _EXTRA = e 


Zphase= (180/!PI)*atan(imaginary(Z),float(Z))
contour, Zphase,x,y, LEVELS=findgen(24)/24.*360.-180. $
 ,/FOLLOW,/OVERPLOT, c_linestyle=1

return
end