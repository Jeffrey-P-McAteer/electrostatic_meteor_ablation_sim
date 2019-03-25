; input the J field in 1 - 3 -D .  Also set the x,y and z scales
@params_in.pro
if (file_info('Jx.bin')).size gt long(4)*nx2 then $
  Jx=readarray('Jx.bin',[nx2,ny2,nz2],order=[2,1,0],/binary,skip=iskip)

if (file_info('Jy.bin')).size gt long(4)*nx2 then $
  Jy=readarray('Jy.bin',[nx2,ny2,nz2],order=[2,1,0],/binary,skip=iskip) $
  else Jy=0

if (file_info('Jz.bin')).size gt long(4)*nx2 then $
  Jz=readarray('Jz.bin',[nx2,ny2,nz2],order=[2,1,0],/binary,skip=iskip) $
  else Jz=0

Jsize=size(Jx)

xv=findgen(Jsize(1))*dx*nout_avg
yv=findgen(Jsize(2))*dy*nout_avg
zv=findgen(Jsize(3))*dz*nout_avg
;Time coordinate
tv=findgen(Jsize(4))*dt*nout*iskip

J2avg=fltarr(Jsize[4])
J2max=fltarr(Jsize[4])

for i=0,Jsize[4]-1 do begin
    J2avg(i)=mean(Jx(*,*,*,i)^2)
    if (ny2 gt 1) then J2avg(i) += mean(Jy(*,*,*,i)^2)
    if (nz2 gt 1) then J2avg(i) += mean(Jz(*,*,*,i)^2)
    
    if (ny2 gt 1 and nz2 gt 1) then $
      J2max(i)=max(Jx(*,*,*,i)^2+Jy(*,*,*,i)^2+Jz(*,*,*,i)^2) $
    else if (ny2 gt 1) then $
      J2max(i)=max(Jx(*,*,*,i)^2+Jy(*,*,*,i)^2) $
    else $
      J2max(i)=max(Jx(*,*,*,i)^2)
    
endfor

plot,tv,J2max,/ylog,linestyle=1, yrange=[min(J2avg[1:*]),max(J2max)], $
  title='J!E2!N Max and Avg', $
  xtitle='time (s)', $
  ytitle='J!E2!N (SI units)'
oplot,tv,J2avg
date_plot,title
end

