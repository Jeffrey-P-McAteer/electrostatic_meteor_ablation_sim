; Routine to plot k vs. t power spectra in 1, 2 and 3 D.
; In 3-D average over the dimension not shown.
; Starts with the potential, phi and generates E^2(k) at various t


@params_in.pro

;== open and assoc the data file ======================================
Openr, unit, 'phi.bin', /GET_LUN
phi_a = Assoc(unit, Fltarr(nz2,ny2,nx2))     
nphi = (Fstat(phi_a)).size / (4L * nz2*ny2*nx2 )            
    
if n_elements(ifirst) eq 0 then ifirst=0
if n_elements(ilast) eq 0 then ilast=nphi-1
if n_elements(iskip) eq 0 then iskip=1

phisize=size(Transpose(convert_endian(phi_a[0]),[2,1,0]))

; Extent to image
ixmin=phisize(1)/2-(phisize(1)*6/16)
ixmax=phisize(1)/2+(phisize(1)*6/16-1)
iymin=phisize(2)/2-(phisize(2)*2/16)
iymax=phisize(2)/2+(phisize(2)*2/16-1)
izmin=phisize(3)/2-(phisize(3)*2/16)
izmax=phisize(3)/2+(phisize(3)*2/16-1)

;Set up the k vector for the axes
kxv=shift(findgen(phisize(1))*2*!PI/(nx*dx),nx2/2)
kxv(0:nx2/2-1)=kxv(0:nx2/2-1)-2*!PI/(dx*nout_avg)
if phisize(0) ge 2 then begin
    kyv=shift(findgen(phisize(2))*2*!PI/(ny*dy),ny2/2)
    kyv(0:ny2/2-1)=kyv(0:ny2/2-1)-2*!PI/(dy*nout_avg)
endif
if phisize(0) ge 3 then begin
    kzv=shift(findgen(phisize(3))*2*!PI/(nz*dz),nz2/2)
    kzv(0:nz2/2-1)=kzv(0:nz2/2-1)-2*!PI/(dz*nout_avg)
endif
;Calculate the simulation volume
vol=nx2*dx
if ny gt 1 then vol=vol*ny2*dy
if nz gt 1 then vol=vol*nz2*dz

;1D plots:
if (phisize(2) eq 1 and phisize(3) eq 1) then begin
    E2=fltarr(nx2,nphi)
    for i =0, nphi,iskip do begin
        E=reform(convert_endian(phi_a[i]))
        E=gradient(E,dx=dx*nout_avg)
        E=shift(fft(E),nx2/2)/vol
        E2(i)=abs(E*conj(E))
    endfor
    label='E!E2!N(k,t):'
    image_plot,reform(E2), kxv, tv, /zlog, /legend, nlabels=4, $
      title=label, xtitle='kx (m)', ytitle='time (s)'
    date_plot,title + ' log 3 decades '

; 2D electric field figures (Assume z is neglected)
endif else if (phisize(3) eq 1) then begin
    for i=ifirst,ilast,iskip do begin
        t=nout*dt*i
        E2=Reform(Transpose(convert_endian(phi_a[i]),[2,1,0]))
        E2=gradient(E2,dx=dx*nout_avg,dy=dy*nout_avg)
        E2=shift(fft(E2),nx2/2,ny2/2)/vol
        E2=total(abs(E2*conj(E2)),3)

        label='E!E2!N(k,t): t='+ strcompress(string(t,format='(G8.4)')) $
          + ' max=' + strcompress(string(max(E2),format='(G8.2)'))
        image_plot,E2(ixmin:ixmax,iymin:iymax), $
          kxv(ixmin:ixmax), kyv(iymin:iymax), $
          /zlog, nlabels=4, /aspect,$
          title=label, $
          xtitle='kx / !4k!3!D0!N', $
          ytitle='ky / !4k!3!D0!N',/legend
        date_plot, title + ' log 3 decades '
endfor


; 3D electric field figures (Sum around k_perp)
endif else begin
; Set up theta_integration:
    pi2dx = !PI / (dx * NOUT_AVG)
    pi2dy = !PI / (dy * NOUT_AVG)
    pi2dz = !PI / (dz * NOUT_AVG)

    kx = FIndgen(nx2) / nx2 * 2.*pi2dx
    greater_pi = Where(kx GT pi2dx, count)
    IF count GT 0 THEN kx[greater_pi] = kx[greater_pi] - 2.*pi2dx

    ky = FIndgen(ny2) / ny2 * 2.*pi2dy
    greater_pi = Where(ky GT pi2dy, count)
    IF count GT 0 THEN ky[greater_pi] = ky[greater_pi] - 2.*pi2dy

    kz = FIndgen(nz2) / nz2 * 2.*pi2dz
    greater_pi = Where(kz GT pi2dz, count)
    IF count GT 0 THEN kz[greater_pi] = kz[greater_pi] - 2.*pi2dz

    bounds = Fltarr(4)
    nperp = ny2 > nz2
    dummy = int_theta(reform(ky # kz, ny2, nz2) ,ky, kz, /grid, $
                     points=nperp, bounds=bounds)

    k_perp_max = bounds[2]
    dkperp = k_perp_max / nperp
    kpv=findgen(nperp)*dkperp

    ; X-Y cross-sections:	
    for i=ifirst,ilast,iskip do begin
        t=nout*dt*i
        E2=Transpose(convert_endian(phi_a[i]),[2,1,0])
        E2=gradient(E2, dx=dx*nout_avg, dy=dy*nout_avg, dz=dz*nout_avg)
        E2=fft(E2,/overwrite)/vol
        E2=total(abs(E2*conj(E2)),4)
        E2int=fltarr(nx2,nperp)
        for ix=0,nx2-1 do $
          E2int(ix,*)=int_theta(reform(E2(ix,*,*)), ky, kz, /grid, /recall)
        E2=shift(E2int,nx2/2-1,0)
        E2=E2(ixmin:ixmax,0:iymax>izmax)

        label='E!E2!N(k,t): 3D avg kx-ky t=' $ 
          + strcompress(string(t,format='(G8.4)')) $
          + ' max=' + strcompress(string(max(E2),format='(G8.2)'))
        image_plot,E2, kxv(ixmin:ixmax), kpv(0:(size(E2))(2)-1), /aspect,  $
          /zlog, nlabels=4, $
          title=label, $
          xtitle='kx / !4k!3!D0!N', $
          ytitle='ky / !4k!3!D0!N'
        date_plot, title + ' log 3 decades '
     endfor

    ; Y-Z cross-section
;xxx    for i=ifirst,ilast,iskip do begin
;xxx        t=nout*dt*i
;xxx        E2=Transpose(convert_endian(phi_a[i]),[2,1,0])
;xxx        E2=gradient(E2, dx=dx*nout_avg, dy=dy*nout_avg, dz=dz*nout_avg)
;xxx        E2=shift(fft(E2,/overwrite),nx2/2,ny2/2,nz2/2,0)/vol
;xxx        E2=total(abs(E2*conj(E2)),4)
;xxx        E2=total(E2,1)/phisize(1)
;xxx
;xxx        E2=E2(iymin:iymax,izmin:izmax)
;xxx
;xxx        label='E!E2!N: 3D avg ky-kz t=' $ 
;xxx          + strcompress(string(t,format='(G8.4)')) $
;xxx          + ' max=' + strcompress(string(max(E2),format='(G8.2)'))
;xxx        image_plot,E2, kyv(iymin:iymax), kzv(izmin:izmax), /aspect,  $
;xxx          /zlog, nlabels=4, $
;xxx          zrange = [0,max(E2)], $
;xxx          title=label, $
;xxx          xtitle='ky / !4k!3!D0!N', $
;xxx          ytitle='kz / !4k!3!D0!N'
;xxx        date_plot, title + ' log 3 decades '
;xxx     endfor


endelse
end
