; input the B field in 1 - 3 -D .  Also set the x,y and z scales
@params_in.pro
if (file_info('Bx.bin')).size gt long(4)*nx2 then $
  Bx=readarray('Bx.bin',[nx2,ny2,nz2],order=[2,1,0],/binary,skip=iskip)

if (file_info('By.bin')).size gt long(4)*nx2 then $
  By=readarray('By.bin',[nx2,ny2,nz2],order=[2,1,0],/binary,skip=iskip) $
  else By=0

if (file_info('Bz.bin')).size gt long(4)*nx2 then $
  Bz=readarray('Bz.bin',[nx2,ny2,nz2],order=[2,1,0],/binary,skip=iskip) $
  else Bz=0

Bsize=size(Bx)

xv=findgen(Bsize(1))*dx*nout_avg
yv=findgen(Bsize(2))*dy*nout_avg
zv=findgen(Bsize(3))*dz*nout_avg
;Time coordinate
tv=findgen(Bsize(4))*dt*nout*iskip

B2avg=fltarr(Bsize[4])
B2max=fltarr(Bsize[4])

for i=0,Bsize[4]-1 do begin
    B2avg(i)=mean(Bx(*,*,*,i)^2)
    if (ny2 gt 1) then B2avg(i) += mean(By(*,*,*,i)^2)
    if (nz2 gt 1) then B2avg(i) += mean(Bz(*,*,*,i)^2)
    
    if (ny2 gt 1 and nz2 gt 1) then $
      B2max(i)=max(Bx(*,*,*,i)^2+By(*,*,*,i)^2+Bz(*,*,*,i)^2) $
    else if (ny2 gt 1) then $
      B2max(i)=max(Bx(*,*,*,i)^2+By(*,*,*,i)^2) $
    else $
      B2max(i)=max(Bx(*,*,*,i)^2)
    
endfor

plot,tv,B2max,/ylog,linestyle=1, yrange=[min(B2avg[1:*]),max(B2max)], $
  title='B!E2!N Max and Avg', $
  xtitle='time (s)', $
  ytitle='B!E2!N (SI units)'
oplot,tv,B2avg
date_plot,title

if n_elements(want_B_plots) eq 0 then want_B_plots=0
if (want_B_plots) then begin

    kxv=shift(findgen(Bsize(1))*2*!PI/(nx*dx),nx2/2)
    kxv(0:nx2/2-1)=kxv(0:nx2/2-1)-2*!PI/(dx*nout_avg)

    kyv=shift(findgen(Bsize(2))*2*!PI/(ny*dy),ny2/2)
    if (ny2 gt 1) then kyv(0:ny2/2-1)=kyv(0:ny2/2-1)-2*!PI/(dy*nout_avg)
    
    kzv=shift(findgen(Bsize(3))*2*!PI/(nz*dz),nz2/2)
    if (nz2 gt 1) then kzv(0:nz2/2-1)=kzv(0:nz2/2-1)-2*!PI/(dz*nout_avg)

    nt=Bsize(4)
    wv=findgen(Bsize(4))*2*!PI/(nt*dt*nout*iskip)
    
    Bfft=(shift(abs(fft(Bx[*,0,0,*])),nx2/2))[*,*,*,0:Bsize(4)/2-1]
    image_plot,Bfft,kxv,wv[0:Bsize(4)/2-1],/zlog,/legend,xtitle='k',ytitle='omega'
    date_plot, title

    Bfft=(shift(abs(fft(By[*,0,0,*])),nx2/2))[*,*,*,0:Bsize(4)/2-1]
    image_plot,Bfft,kxv,wv[0:Bsize(4)/2-1],/zlog,/legend,xtitle='k',ytitle='omega'
    date_plot, title

;    Bfft=(shift(abs(fft(Bz[*,0,0,*])),nx2/2))[*,*,*,0:Bsize(4)/2-1]
;    image_plot,Bfft,kxv,wv[0:Bsize(4)/2-1],/zlog,/legend,xtitle='k',ytitle='omega'
;    date_plot, title


endif

end
