; Routine to plot k vs. t power spectra in 1, 2 and 3 D.
; Starts with the potential, phi and generates E^2(k) at various t

@phi3d

Esize=size(phi)

if n_elements(ifirst) eq 0 then ifirst=0
if n_elements(ilast) eq 0 then ilast=Esize(4)-1
if n_elements(iskip) eq 0 then iskip=1
if n_elements(avg_3d) eq 0 then avg_3d=0
if n_elements(decades) eq 0 then decades=3
if (avg_3d eq 0) then avg_string='' else avg_string='avg'

; Extent to image
ixmin=Esize(1)/2-(Esize(1)*6/16)
ixmax=Esize(1)/2+(Esize(1)*6/16)-1
iymin=Esize(2)/2
iymax=Esize(2)/2+(Esize(2)*7/16)-1
izmin=Esize(3)/2-(Esize(3)/2)
izmax=Esize(3)/2+(Esize(3)/2)-1


;Set up the k vector for the axes
;and construct k^2 array
kxv=shift(findgen(Esize(1))*2*!PI/(nx*dx),nx2/2)
kxv(0:nx2/2-1)=kxv(0:nx2/2-1)-2*!PI/(dx*nout_avg)
k2=kxv^2
if Esize(0) gt 2 then begin
    kyv=shift(findgen(Esize(2))*2*!PI/(ny*dy),ny2/2)
    kyv(0:ny2/2-1)=kyv(0:ny2/2-1)-2*!PI/(dy*nout_avg)
    k2=((fltarr(ny2)+1.)##kxv^2+kyv^2##(fltarr(nx2)+1.))
endif
if Esize(0) gt 3 and nz2 gt 1 then begin
    kzv=shift(findgen(Esize(3))*2*!PI/(nz*dz),nz2/2)
    kzv(0:nz2/2-1)=kzv(0:nz2/2-1)-2*!PI/(dz*nout_avg)
    k2store=k2
    k2=fltarr(nx2,ny2,nz2)
    for iz=0, nz2-1 do k2(*,*,iz)= k2store + kzv(iz)^2
endif


;1D plots:
if (Esize(2) eq 1 and Esize(3) eq 1) then begin
    E2=reform(phi)
    for i =0, Esize(4) do begin
        E=gradient(E(*,0,0,i),dx=dx*nout_avg)
        E=shift(fft(E),nx2/2)
        E2(i)=abs(E*conj(E))
    endfor
    label='E!E2!N(k,t):'
    image_plot, reform(E2), kxv, tv, /zlog, nlabels=decades+1, $
      title=label, xtitle='k!Dx!N / !4k!3!D0!N', ytitle='time (s)'
    date_plot,title + ' log ' + strcompress(string(decades)) + ' decades '

; 2D electric field figures (Assume z is neglected)
endif else if (Esize(3) eq 1) then begin
    for i=ifirst,ilast,iskip do begin
        t=nout*dt*i
        
        phifft=shift(fft(reform(phi(*,*,*,i))),nx2/2,ny2/2)
        E2=float(phifft*conj(phifft))*k2

        label='E!E2!N(k,t): t='+ strcompress(string(t,format='(G8.4)')) $
          + ' max=' + strcompress(string(max(E2),format='(G8.2)'))
        image_plot,E2(ixmin:ixmax,iymin:iymax), kxv(ixmin:ixmax), kyv(iymin:iymax), /aspect, $
          /zlog, nlabels=decades+1, $
          zrange = [0,max(E2)], $
          title=label, $
          xtitle='k!Dx!N / !4k!3!D0!N', $
          ytitle='ky / !4k!3!D0!N'
        date_plot, title + ' log ' + strcompress(string(decades)) + ' decades '
endfor

; 3D electric field figures (Sum around k_perp)
endif else begin

    for i=ifirst,ilast,iskip do begin
        t=nout*dt*i

        phifft=shift(fft(reform(phi(*,*,*,i))),nx2/2,ny2/2,nz2/2)
        E2a=float(phifft*conj(phifft))*k2

;
;       Eliminate the strictly perpendicular modes to evaluate the
;       amplitude of the non perpendicular modes
        E2a(nx2/2,*,*)=0

                                ; X-Y cross-sections:	

        if avg_3d then E2=total(E2a,3)/Esize(3) $
        else E2=reform(E2a(*,*,Esize(3)/2))
        
        label='E!E2!N(k,t): 3D '+ avg_string + ' x-y t=' $ 
          + strcompress(string(t,format='(G8.4)')) $
          + ' max=' + strcompress(string(max(E2),format='(G8.2)'))
        image_plot, E2(ixmin:ixmax,iymin:iymax), kxv(ixmin:ixmax), kyv(iymin:iymax), $
          /aspect, /zlog, nlabels=decades+1, $
          zrange = [0,max(E2)], $
          title=label, $
          xtitle='k!Dx!N / !4k!3!D0!N', $
          ytitle='ky / !4k!3!D0!N'
        date_plot, title + ' log ' + strcompress(string(decades)) + ' decades '

                                ; X-Z cross-section

        if avg_3d then E2=total(E2a,2)/Esize(2) $
        else E2=reform(E2a(*,Esize(2)/2,*))

        label='E!E2!N: 3D '+ avg_string + ' x-z t=' $ 
          + strcompress(string(t,format='(G8.4)')) $
          + ' max=' + strcompress(string(max(E2),format='(G8.2)'))
        image_plot, E2(ixmin:ixmax,izmin:izmax), kxv(ixmin:ixmax), kzv(izmin:izmax), $
          /aspect, /zlog, nlabels=decades+1, $
          zrange = [0,max(E2)], $
          title=label, $
          xtitle='k!Dx!N / !4k!3!D0!N', $
          ytitle='k!Dz!N / !4k!3!D0!N'
        date_plot, title + ' log ' + strcompress(string(decades)) + ' decades '

                                ; Y-Z cross-section

        if avg_3d then E2=total(E2a,1)/Esize(1) $
        else E2=reform(E2a(Esize(1)/2,*,*))

        label='E!E2!N: 3D '+ avg_string + ' y-z t=' $ 
          + strcompress(string(t,format='(G8.4)')) $
          + ' max=' + strcompress(string(max(E2),format='(G8.2)'))
        image_plot, E2(iymin:iymax,izmin:izmax), kyv(iymin:iymax), kzv(izmin:izmax), $
          /aspect, /zlog, nlabels=decades+1, $
          zrange = [0,max(E2)], $
          title=label, $
          xtitle='k!Dy!N / !4k!3!D0!N', $
          ytitle='k!Dz!N / !4k!3!D0!N'
        date_plot, title + ' log ' + strcompress(string(decades)) + ' decades '

     endfor


endelse
end
