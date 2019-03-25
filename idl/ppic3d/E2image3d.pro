; Generate images of E^2
; Use the variables ifirst,ilast, iskip to determine how many images

;Calculate E2 in real space:
if n_elements(phi) le 2 then begin
@phi3d
endif

Esize=size(phi)

if n_elements(ifirst) eq 0 then ifirst=0
if n_elements(ilast) eq 0 then ilast=0
if ilast eq 0 then ilast=Esize(4)-1
if n_elements(iskip) eq 0 then iskip=1


;1D plots:
if (Esize(2) eq 1 and Esize(3) eq 1) then begin
    E2=reform(phi)
    for i =0, Esize(4)-1 do E2(i)=gradient(phi(*,0,0,i),dx=dx*nout_avg)^2
    label='E!E2!N:'
    image_plot,reform(E2), xv, tv, /zlog, /legend, nlabels=4, $
      title=label, xtitle='x (m)', ytitle='time (s)'
    date_plot,title + ' log 3 decades '

; 2D electric field figures (Assume z is neglected)
endif else if (Esize(3) eq 1) then begin
    for i=ifirst,ilast,iskip do begin
        t=nout*dt*i
        E2=total(gradient(phi(*,*,0,i),dx=dx*nout_avg,dy=dy*nout_avg)^2,3)
        label='E!E2!N: t='+ strcompress(string(t,format='(G8.4)')) $
          + ' max=' + strcompress(string(max(E2),format='(G8.2)'))
        image_plot,E2, xv, yv, /aspect, $
          /zlog, nlabels=3, $
          zrange = [0,max(E2)], $
          title=label, $
          xtitle='x / !4k!3!D0!N', $
          ytitle='y / !4k!3!D0!N',/legend
        date_plot, title + ' log 2 decades '
endfor

; 3D electric field figures (Cross-sections)
endif else begin
    ; X-Y cross-sections:	
    for i=ifirst,ilast,iskip do begin
        t=nout*dt*i
        E2=total(gradient(phi(*,*,*,i), dx=dx*nout_avg, $
                          dy=dy*nout_avg, dz=dz*nout_avg)^2, 4)

        label='E!E2!N: 3D cross section x-y t=' $ 
          + strcompress(string(t,format='(G8.4)')) $
          + ' max=' + strcompress(string(max(E2),format='(G8.2)'))
        image_plot,reform(E2(*,*,Esize(3)/2)), xv, yv, /aspect, $
          /zlog, nlabels=3, $
          zrange = [0,max(E2)], $
          title=label, $
          xtitle='x / !4k!3!D0!N', $
          ytitle='y / !4k!3!D0!N',/legend
        date_plot, title + ' log 2 decades '
     endfor

    ; X-Z cross-section
    for i=ifirst,ilast,iskip do begin
        t=nout*dt*i
        E2=total(gradient(phi(*,*,*,i), dx=dx*nout_avg, $
                          dy=dy*nout_avg, dz=dz*nout_avg)^2, 4)

        label='E!E2!N: 3D cross section x-z t=' $ 
          + strcompress(string(t,format='(G8.4)')) $
          + ' max=' + strcompress(string(max(E2),format='(G8.2)'))
        image_plot,reform(E2(*,Esize(2)/2,*)), xv, zv, /aspect, $
          /zlog, nlabels=3, $
          zrange = [0,max(E2)], $
          title=label, $
          xtitle='x / !4k!3!D0!N', $
          ytitle='y / !4k!3!D0!N',/legend
        date_plot, title + ' log 2 decades '
     endfor

     ; Y-Z cross-section
     for i=ifirst,ilast,iskip do begin
        t=nout*dt*i
        E2=total(gradient(phi(*,*,*,i), dx=dx*nout_avg, $
                          dy=dy*nout_avg, dz=dz*nout_avg)^2, 4)

        label='E!E2!N: 3D cross section y-z t=' $ 
          + strcompress(string(t,format='(G8.4)')) $
          + ' max=' + strcompress(string(max(E2),format='(G8.2)'))
        image_plot,reform(E2(Esize(1)/2,*,*)), yv, zv, /aspect, $
          /zlog, nlabels=3, $
          zrange = [0,max(E2)], $
          title=label, $
          xtitle='x / !4k!3!D0!N', $
          ytitle='y / !4k!3!D0!N',/legend
        date_plot, title + ' log 2 decades '
    endfor

endelse
end
