function rk4c, y, dydx, x, h, func
; 
; NAME:
;       RK4C
;
; PURPOSE:
;       A fourth order Runge-Kutta method allowing complex (or any
;       other) type arguments.
;


hh=0.5*h
h6=h/6.0
xh=x+hh
yt=y+hh*dydx
dyt=call_function(func,xh,yt)
yt=y+hh*dyt
dym=call_function(func,xh,yt)
yt=y+h*dym
dym=dym+dyt
dyt=call_function(func,x+h,yt)
return,y+h6*(dydx+dyt+2.0*dym)

end
