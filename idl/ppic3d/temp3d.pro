;Calculate and plot the temperatures (only works in 2D)

if n_elements(vdist0) lt 2 then begin
@vdist3d
endif

; We need the Boltzmann constant will usually be mks or 1
if (n_elements(kb) eq 0) then if (md0<1e-8) then kb=1.38e-23 else kb=1

izsize=((size(vdist0))(3))

vxs=findgen(pnvx0)/float(pnvx0)*(pvxmax0-pvxmin0) + pvxmin0
if n_elements(pnvy0) gt 0 then vys=findgen(pnvy0)/float(pnvy0)*(pvymax0-pvymin0) + pvymin0
if izsize gt 1 then vzs=findgen(pnvz0)/float(pnvz0)*(pvzmax0-pvzmin0) + pvzmin0 else vzs=0

vxmat2d=vxs#(fltarr(pnvy0)+1.0)
vymat2d=vys##(fltarr(pnvx0)+1.0)

vxmat=vdist0(*,*,*,0)
vymat=vdist0(*,*,*,0)
vzmat=vdist0(*,*,*,0)
for iz= 0, izsize -1 do begin
    vxmat(*,*,iz)=vxmat2d
    vymat(*,*,iz)=vymat2d
    vzmat(*,*,iz)=vzs(iz)*(fltarr(pnvx0,pnvy0)+1.0)
endfor

vsize=size(vdist0)
time_iter=vsize(4)
Tempx0=fltarr(time_iter)
Tempy0=fltarr(time_iter)
Tempz0=fltarr(time_iter)
time0=findgen(time_iter)*dt*nout*iskip*vdist_out_subcycle0
for i=0,time_iter-1 do begin
    avgvx=mean(vdist0(*,*,*,i)*vxmat)/n0d0
    Tempx0(i)=mean(vdist0(*,*,*,i)*(vxmat-avgvx)^2)/n0d0*md0/kb
    avgvy=mean(vdist0(*,*,*,i)*vymat)/n0d0
    Tempy0(i)=mean(vdist0(*,*,*,i)*(vymat-avgvy)^2)/n0d0*md0/kb
    avgvz=mean(vdist0(*,*,*,i)*vzmat)/n0d0
    Tempz0(i)=mean(vdist0(*,*,*,i)*(vzmat-avgvz)^2)/n0d0*md0/kb
endfor

if (ndist ge 2) then begin
    vxs=findgen(pnvx1)/float(pnvx1)*(pvxmax1-pvxmin1) + pvxmin1
    vys=findgen(pnvy1)/float(pnvy1)*(pvymax1-pvymin1) + pvymin1
    if izsize gt 1 then vzs=findgen(pnvz1)/float(pnvz1)*(pvzmax1-pvzmin1) + pvzmin1 else vzs=0
    vxmat2d=vxs#(fltarr(pnvy1)+1.0)
    vymat2d=vys##(fltarr(pnvx1)+1.0)
    vxmat=vdist1(*,*,*,0)
    vymat=vdist1(*,*,*,0)
    vzmat=vdist1(*,*,*,0)
    izsize=((size(vdist1))(3))
    for iz= 0, izsize -1 do begin
        vxmat(*,*,iz)=vxmat2d
        vymat(*,*,iz)=vymat2d
        vzmat(*,*,iz)=vzs(iz)*(fltarr(pnvx1,pnvy1)+1.0)
    endfor
    
    vsize=size(vdist1)
    Tempx1=fltarr(vsize(4))
    Tempy1=fltarr(vsize(4))
    Tempz1=fltarr(vsize(4))
    time1=findgen(vsize(4))*dt*nout*iskip*vdist_out_subcycle1
    for i=0,vsize(4)-1 do begin
        avgvx=mean(vdist1(*,*,*,i)*vxmat)/n0d1
        Tempx1(i)=mean(vdist1(*,*,*,i)*(vxmat-avgvx)^2)/n0d1*md1/kb
        avgvy=mean(vdist1(*,*,*,i)*vymat)/n0d1
        Tempy1(i)=mean(vdist1(*,*,*,i)*(vymat-avgvy)^2)/n0d1*md1/kb
        avgvz=mean(vdist1(*,*,*,i)*vzmat)/n0d1
        Tempz1(i)=mean(vdist1(*,*,*,i)*(vzmat-avgvz)^2)/n0d1*md1/kb
    endfor

    ;Calculate the combined dist1 and dist2 temp:
endif

; Plot the temps:
maxtemp=max([Tempx0,Tempy0,Tempz0])
if (ndist ge 2) then maxtemp=max([maxtemp,max([Tempx1,Tempy1,Tempz1])])
;if (ndist ge 3) then maxtemp=max([maxtemp,max([Tempx2,Tempy2])])
if (n_elements(Tempx12) ge 2) then maxtemp=max([maxtemp,max([Tempx12,Tempy12])])

mintemp=min([Tempx0,Tempy0,Tempz0])
if (ndist ge 2) then mintemp=min([mintemp,min([Tempx1,Tempy1,Tempz1])])
;if (ndist ge 3) then mintemp=min([mintemp,min([Tempx2,Tempy2])])
if (n_elements(Tempx12) ge 2) then mintemp=min([mintemp,min([Tempx12,Tempy12])])
plot, time0, Tempx0, ystyle=1, yrange=[mintemp,maxtemp],  $
  title='Temperatures: ' + title, ytitle='kT (J)', xtitle='time (s)'
xyouts,time0(vsize(4)-1),Tempx0(vsize(4)-1),' x0'
oplot,time0,Tempy0
xyouts,time0(vsize(4)-1),Tempy0(vsize(4)-1),' y0'
if vsize(3) gt 1 then begin
    oplot,time0,Tempz0
    xyouts,time0(vsize(4)-1),Tempz0(vsize(4)-1),' z0'
endif

if (ndist ge 2) then begin
    oplot, time1, Tempx1
    xyouts,time1(vsize(4)-1),Tempx1(vsize(4)-1),' x1'
    oplot,time1,Tempy1
    xyouts,time1(vsize(4)-1),Tempy1(vsize(4)-1),' y1'
    if vsize(3) gt 1 then begin
        oplot,time1,Tempz1
        xyouts,time1(vsize(4)-1),Tempz1(vsize(4)-1),' z1'
    endif
    if (n_elements(Tempx12) ge 2) then begin
        oplot, time1, Tempx12
        xyouts,time1(vsize(4)-1),Tempx12(vsize(4)-1),' x12'
        oplot,time1,Tempy12
        xyouts,time1(vsize(4)-1),Tempy12(vsize(4)-1),' y12'
    endif
endif
;if (ndist ge 3) then begin
;    oplot, time, Tempx2
;    xyouts,time(vsize(4)-1),Tempx2(vsize(4)-1),' x2'
;    oplot,time,Tempy2
;    xyouts,time(vsize(4)-1),Tempy2(vsize(4)-1),' y2'
;endif
date_plot

end
