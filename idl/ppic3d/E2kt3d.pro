; Routine to plot k vs. t power spectra in 1, 2 and 3 D.
; Starts with the potential, phi and generates E^2(k) at various t

; Set zlog_on=1 for logrithmic scales (default) and zlog_on=0 for linear
if n_elements(zlog_on) eq 0 then zlog_on=1

; Set legend_on=0 for no legend (default) or legend_on=1 for a legend
if n_elements(legend_on) eq 0 then legend_on=0

Esize=size(phi)

if n_elements(ifirst) eq 0 then ifirst=0
if n_elements(ilast) eq 0 then ilast=Esize(4)-1
if n_elements(iskip) eq 0 then iskip=1
if n_elements(avg_3d) eq 0 then avg_3d=0
if n_elements(decades) eq 0 then decades=3
if (avg_3d eq 0) then avg_string='' else avg_string='avg'


;Set up the k vector for the axes
;and construct k^2 array
kxv=shift(findgen(Esize(1))*2*!PI/(nx*dx),nx2/2)
kxv(0:nx2/2-1)=kxv(0:nx2/2-1)-2*!PI/(dx*nout_avg)
k2=kxv^2
if Esize(0) gt 2 and ny2 gt 1 then begin
    kyv=shift(findgen(Esize(2))*2*!PI/(ny*dy),ny2/2)
    kyv(0:ny2/2-1)=kyv(0:ny2/2-1)-2*!PI/(dy*nout_avg)
    k2=((fltarr(ny2)+1.)##kxv^2+kyv^2##(fltarr(nx2)+1.))
    if Esize(0) gt 3 and nz2 gt 1 then begin
        kzv=shift(findgen(Esize(3))*2*!PI/(nz*dz),nz2/2)
        kzv(0:nz2/2-1)=kzv(0:nz2/2-1)-2*!PI/(dz*nout_avg)
        k2store=k2
        k2=fltarr(nx2,ny2,nz2)
        for iz=0, nz2-1 do k2(*,*,iz)= k2store + kzv(iz)^2
    endif
endif


;1D plots:
if (Esize(2) eq 1 and Esize(3) eq 1) then begin
    E2=reform(phi)
    for i = 0, Esize(4)-1 do begin
        E=gradient(phi(*,0,0,i),dx=dx*nout_avg)
        E=shift(fft(E),nx2/2)
        E2(*,i)=abs(E*conj(E))
    endfor
    label='E!E2!N(k,t):'
    image_plot, reform(E2), kxv, tv, zlog=zlog_on, nlabels=decades+1, $
      title=label, xtitle='k!Dx!N / !4k!3!D0!N', ytitle='time (s)', $
      legend=legend_on,/interp
    if zlog_on then date_plot,title + ' log ' + strcompress(string(decades)) + ' decades '$
      else date_plot,title + ' linear'

; 2D electric field figures (Assume z is neglected)
endif else if (Esize(3) eq 1) then begin
                                ;Only plot the inner component of k space
    phifft=shift(reform(fft(phi(*,*,*,phisize(4)-1)),nx2,ny2,nz2),nx2/2,ny2/2,nz2/2)
    E2=float(phifft*conj(phifft))*k2
    range=image_range(E2 gt max(E2)/100.)
                                ; Extent to image
    ixmin=range[0,0]
    ixmax=range[0,1]
    if ixmax eq ixmin then ixmax=ixmin+1
    iymin=Esize(2)/2
    iymax=range[1,1]
    if iymax eq iymin then iymax=iymin+1
    if size(range,/n_dimensions) ge 3 then begin
        izmin=range[2,0]
        izmax=range[2,1]
    endif


    for i=ifirst,ilast,iskip do begin
        t=nout*dt*i
        
        phifft=shift(fft(reform(phi(*,*,*,i))),nx2/2,ny2/2)
        E2=float(phifft*conj(phifft))*k2

        label='E!E2!N(k,t): t='+ strcompress(string(t,format='(G8.4)')) $
          + ' max=' + strcompress(string(max(E2),format='(G8.2)'))
        image_plot,E2(ixmin:ixmax,iymin:iymax), $
          kxv(ixmin:ixmax), kyv(iymin:iymax), /aspect, $
          zlog=zlog_on, nlabels=decades+1, legend=legend_on, $
          zrange = [0,max(E2)], $
          title=label, $
          xtitle='k!Dx!N / !4k!3!D0!N', $
          ytitle='ky / !4k!3!D0!N' ,/interp
        if zlog_on then date_plot,title + ' log ' + strcompress(string(decades)) + ' decades '$
        else date_plot,title + ' linear'

endfor

; 3D electric field figures (Sum around k_perp)
endif else begin
    for i=ifirst,ilast,iskip do begin
        t=nout*dt*i

        phifft=shift(fft(reform(phi(*,*,*,i))),nx2/2,ny2/2,nz2/2)
        E2a=float(phifft*conj(phifft))*k2

                                ; X-Y cross-sections:	

        if avg_3d then E2=total(E2a,3)/Esize(3) $
        else E2=reform(E2a(*,*,Esize(3)/2))
        

                                ; Extent to image
        range=image_range(E2 gt max(E2)/100.)
        ixmin=range[0,0]
        ixmax=range[0,1]
        if ixmax eq ixmin then ixmax=ixmin+1
        iymin=Esize(2)/2
        iymax=range[1,1]
        if iymax eq iymin then iymax=iymin+1

        label='E!E2!N(k,t): 3D '+ avg_string + ' x-y t=' $ 
          + strcompress(string(t,format='(G8.4)')) $
          + ' max=' + strcompress(string(max(E2),format='(G8.2)'))
        image_plot, E2(ixmin:ixmax,iymin:iymax), $
          kxv(ixmin:ixmax), kyv(iymin:iymax), $
          /aspect, zlog=zlog_on, nlabels=decades+1, legend=legend_on, $
          zrange = [0,max(E2)], $
          title=label, $
          xtitle='k!Dx!N / !4k!3!D0!N', $
          ytitle='ky / !4k!3!D0!N' ,/interp

        if zlog_on then date_plot,title + ' log ' + strcompress(string(decades)) + ' decades '$
        else date_plot,title + ' linear'


                                ; X-Z cross-section

        if avg_3d then E2=total(E2a,2)/Esize(2) $
        else E2=reform(E2a(*,Esize(2)/2,*))

                                ; Extent to image
        range=image_range(E2 gt max(E2)/100.)
        ixmin=range[0,0]
        ixmax=range[0,1]
        if ixmax eq ixmin then ixmax=ixmin+1
        izmin=Esize(2)/2
        izmax=range[1,1]
        if izmax eq izmin then izmax=izmin+1

        label='E!E2!N: 3D '+ avg_string + ' x-z t=' $ 
          + strcompress(string(t,format='(G8.4)')) $
          + ' max=' + strcompress(string(max(E2),format='(G8.2)'))
        image_plot, E2(ixmin:ixmax,izmin:izmax), $
          kxv(ixmin:ixmax), kzv(izmin:izmax), $
          /aspect, zlog=zlog_on, nlabels=decades+1, legend=legend_on, $
          zrange = [0,max(E2)], $
          title=label, $
          xtitle='k!Dx!N / !4k!3!D0!N', $
          ytitle='k!Dz!N / !4k!3!D0!N' ,/interp
        if zlog_on then date_plot,title + ' log ' + strcompress(string(decades)) + ' decades '$
        else date_plot,title + ' linear'

                                ; Y-Z cross-section

        if avg_3d then E2=total(E2a,1)/Esize(1) $
        else E2=reform(E2a(Esize(1)/2,*,*))

                                ; Extent to image
        range=image_range(E2 gt max(E2)/100.)
        iymin=range[0,0]
        iymax=range[0,1]
        if iymax eq iymin then iymax=iymin+1
        izmin=Esize(2)/2
        izmax=range[1,1]
        if izmax eq izmin then izmax=izmin+1


        label='E!E2!N: 3D '+ avg_string + ' y-z t=' $ 
          + strcompress(string(t,format='(G8.4)')) $
          + ' max=' + strcompress(string(max(E2),format='(G8.2)'))
        image_plot, E2(iymin:iymax,izmin:izmax), $
          kyv(iymin:iymax), kzv(izmin:izmax), $
          /aspect, zlog=zlog_on, nlabels=decades+1, legend=legend_on,$
          zrange = [0,max(E2)], $
          title=label, $
          xtitle='k!Dy!N / !4k!3!D0!N', $
          ytitle='k!Dz!N / !4k!3!D0!N' ,/interp
        if zlog_on then date_plot,title + ' log ' + strcompress(string(decades)) + ' decades '$
        else date_plot,title + ' linear'


     endfor


endelse
end
