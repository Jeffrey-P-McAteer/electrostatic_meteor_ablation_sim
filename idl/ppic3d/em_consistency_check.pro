; This routine runs checks pempic for divB=0 and div E =p/eps0

Bsize=size(Bx)
; Let's calculate the sumsquare of this quantity across the mesh
nB=Bsize[Bsize[0]]
divBtot=dblarr(nB)
divBmax=dblarr(nB)
for i=0,nB-1 do begin
    divB=divergence(Bx[*,*,*,i],By[*,*,*,i],Bz[*,*,*,i],dx=dx,dy=dy,dz=dz)
    divBtot(i)=total(divB^2)
    divBmax(i)=max((divB^2))
endfor
divBtot /= nx*dx
if (ndim_space ge 2) then divBtot /= ny*dy
if (ndim_space ge 3) then divBtot /= nz*dz

divBmax /= dx
if (ndim_space ge 2) then divBmax /= dy
if (ndim_space ge 3) then divBmax /= dz

;Assume ions and electrons are the only species
den0=0
den1=0
@den3d

divEtot=dblarr(nB)
divEmax=dblarr(nB)
for i=0,nB-1 do begin
    divE=divergence(Ex[*,*,*,i],Ey[*,*,*,i],Ez[*,*,*,i],dx=dx,dy=dy,dz=dz)
    rho=qd0*den0[*,*,*,i]+qd1*den1[*,*,*,i]
    divEtot(i)=total((divE-rho/eps)^2)
    divEmax(i)=max(((divE-rho/eps)^2))
endfor
divEtot /= nx*dx
if (ndim_space ge 2) then divEtot /= ny*dy
if (ndim_space ge 3) then divEtot /= nz*dz

divEmax /= dx
if (ndim_space ge 2) then divEmax /= dy
if (ndim_space ge 3) then divEmax /= dz

if !D.Name EQ "X" THEN window,0
plot,tv,divEtot,/ylog,$
  title='Avg of divB^2 and (divE-rho/eps)^2 per volume', xtitle='time'
xyouts,tv(nB-1),divEtot(nB-1),'avg divE'

oplot,tv,divBtot,linestyle=2
xyouts,tv(nB-1),divBtot(nB-1),'avg divB'

if !D.Name EQ "X" THEN window,1
plot,tv,divEmax,/ylog, $
  title='Max of divB^2 and (divE-rho/eps)^2 per volume', xtitle='time'
xyouts,tv(nB-1),divEmax(nB-1),'max divE'


oplot,tv,divBmax,linestyle=2
xyouts,tv(nB-1),divBmax(nB-1),'max divB'


end
