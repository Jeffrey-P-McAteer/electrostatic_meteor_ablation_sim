; input the E field in 1 - 3 -D .  Also set the x,y and z scales
@params_in.pro
if (file_info('Ex.bin')).size gt long(4)*nx2 then $
  Ex=readarray('Ex.bin',[nx2,ny2,nz2],order=[2,1,0],/binary,skip=iskip)

if (file_info('Ey.bin')).size gt long(4)*nx2 then $
  Ey=readarray('Ey.bin',[nx2,ny2,nz2],order=[2,1,0],/binary,skip=iskip) $
  else Ey=0

if (file_info('Ez.bin')).size gt long(4)*nx2 then $
  Ez=readarray('Ez.bin',[nx2,ny2,nz2],order=[2,1,0],/binary,skip=iskip) $
  else Ez=0

Esize=size(Ex)

xv=findgen(Esize(1))*dx*nout_avg
yv=findgen(Esize(2))*dy*nout_avg
zv=findgen(Esize(3))*dz*nout_avg
;Time coordinate
tv=findgen(Esize(4))*dt*nout*iskip

E2avg=fltarr(Esize[4])
E2max=fltarr(Esize[4])

for i=0,Esize[4]-1 do begin
    E2avg(i)=mean(Ex(*,*,*,i)^2)
    if (ny2 gt 1) then E2avg(i) += mean(Ey(*,*,*,i)^2)
    if (nz2 gt 1) then E2avg(i) += mean(Ez(*,*,*,i)^2)
    
    if (ny2 gt 1 and nz2 gt 1) then $
      E2max(i)=max(Ex(*,*,*,i)^2+Ey(*,*,*,i)^2+Ez(*,*,*,i)^2) $
    else if (ny2 gt 1) then $
      E2max(i)=max(Ex(*,*,*,i)^2+Ey(*,*,*,i)^2) $
    else $
      E2max(i)=max(Ex(*,*,*,i)^2)
    
endfor

plot,tv,E2max,/ylog,linestyle=1, yrange=[min(E2avg[1:*]),max(E2max)], $
  title='E!E2!N Max and Avg', $
  xtitle='time (s)', $
  ytitle='E!E2!N (SI units)'
oplot,tv,E2avg
date_plot,title

if n_elements(want_E_plots) eq 0 then want_E_plots=0
if (want_E_plots) then begin

    kxv=shift(findgen(Esize(1))*2*!PI/(nx*dx),nx2/2)
    kxv(0:nx2/2-1)=kxv(0:nx2/2-1)-2*!PI/(dx*nout_avg)

    kyv=shift(findgen(Esize(2))*2*!PI/(ny*dy),ny2/2)
    if (ny2 gt 1) then kyv(0:ny2/2-1)=kyv(0:ny2/2-1)-2*!PI/(dy*nout_avg)
    
    kzv=shift(findgen(Esize(3))*2*!PI/(nz*dz),nz2/2)
    if (nz2 gt 1) then kzv(0:nz2/2-1)=kzv(0:nz2/2-1)-2*!PI/(dz*nout_avg)

    wv=findgen(Esize(4)/2)*2*!PI/(nt*dt*nout*iskip)
    
    Efft=(shift(abs(fft(Ey[*,0,0,*])),nx2/2))[*,*,*,0:Esize(4)/2-1]
    image_plot,Efft,kxv,wv,/zlog,/legend,xtitle='k',ytitle='omega'
    date_plot, title

    Efft=(shift(abs(fft(Ez[*,0,0,*])),nx2/2))[*,*,*,0:Esize(4)/2-1]
    image_plot,Efft,kxv,wv,/zlog,/legend,xtitle='k',ytitle='omega'
    date_plot, title

    Efft=(shift(abs(fft(Ez[*,0,0,*])),nx2/2))[*,*,*,0:Esize(4)/2-1]
    image_plot,Efft,kxv,wv,/zlog,/legend,xtitle='k',ytitle='omega'
    date_plot, title


endif

end

