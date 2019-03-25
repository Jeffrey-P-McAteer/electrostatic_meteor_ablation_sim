; This routine averages iavg images of den1 together and outputs them
; as an array

; Set zlog_on=1 for logrithmic scales (default) and zlog_on=0 for linear
if n_elements(zlog_on) eq 0 then zlog_on=1
if n_elements(decades) eq 0 then decades=3
legend_on=1

nt=densize(4)
den1ktavg=fltarr(densize[1],densize[2],densize[3],(densize[4]-1)/iavg+1)

for i=0,nt-1 do begin
    den1kt=shift(reform(abs(fft(den1(*,*,*,i))),nx2,ny2,nz2),nx2/2,ny2/2,nz2/2)
    den1ktavg(*,*,*,i/iavg) += den1kt
endfor
den1ktavg /= iavg
if (i MOD iavg) ne 0 then den1ktavg(*,*,*,(i-1)/iavg) *= iavg/(i MOD iavg)



;Set up the k vector for the axes
;and construct k^2 array
kxv=shift(findgen(densize(1))*2*!PI/(nx*dx),nx2/2)
kxv(0:nx2/2-1)=kxv(0:nx2/2-1)-2*!PI/(dx*nout_avg)
k2=kxv^2
if densize(0) gt 2 and ny2 gt 1 then begin
    kyv=shift(findgen(densize(2))*2*!PI/(ny*dy),ny2/2)
    kyv(0:ny2/2-1)=kyv(0:ny2/2-1)-2*!PI/(dy*nout_avg)
    k2=((fltarr(ny2)+1.)##kxv^2+kyv^2##(fltarr(nx2)+1.))
    if Densize(0) gt 3 and nz2 gt 1 then begin
        kzv=shift(findgen(densize(3))*2*!PI/(nz*dz),nz2/2)
        kzv(0:nz2/2-1)=kzv(0:nz2/2-1)-2*!PI/(dz*nout_avg)
        k2store=k2
        k2=fltarr(nx2,ny2,nz2)
        for iz=0, nz2-1 do k2(*,*,iz)= k2store + kzv(iz)^2
    endif
endif

for ip=0,(densize[4]-1)/iavg do begin
    t_start=nout*dt*ip*iavg
    t_stop=min([nout*dt*(ip+1)*iavg,nout*dt*densize[4]])

                                ;Find the range over which there is interesting content:
    range=image_range(den1ktavg[*,*,0,ip] gt max(den1ktavg[*,*,0,ip])/100.)
    ixmin=nx2/2
    ixmax=range[0,1]
    iymin=range[1,0]
    iymax=range[1,1]


    label='n!De!N(k): t='+ strcompress(string(t_start,format='(G8.4)')) $
      + ' - ' + strcompress(string(t_stop,format='(G8.4)')) $
      + ' max=' + strcompress(string(max(den1ktavg),format='(G8.2)'))

    image_plot,den1ktavg(ixmin:ixmax,iymin:iymax,0,ip), $
      kxv(ixmin:ixmax), kyv(iymin:iymax), /aspect, $
      zlog=zlog_on, nlabels=decades+1, legend=legend_on, $
      zrange = [0,max(den1ktavg[*,*,0,ip])], $
      title=label, $
      xtitle='k!Dx!N / !4k!3!D0!N', $
      ytitle='ky / !4k!3!D0!N'

    if zlog_on then date_plot,title + ' log ' + strcompress(string(decades)) + ' decades '$
    else date_plot,title + ' linear'
    
    endfor

end
