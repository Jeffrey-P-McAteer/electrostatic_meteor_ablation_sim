; Generate images of E^2
; Use the variables ifirst,ilast, iskip to determine how many images

if ( n_elements(A2) le 2)  then message, 'Please run envelope3d first'

Esize=size(phi)

if n_elements(ifirst) eq 0 then ifirst=0
if n_elements(ilast) eq 0 then ilast=Esize(4)-1
if n_elements(iskip) eq 0 then iskip=1


;1D plots:
if (Esize(2) eq 1 and Esize(3) eq 1) then begin
    label='A!E2!N:'
    image_plot,reform(A2), xv, tv, /zlog, /legend, nlabels=4, $
      title=label, xtitle='x (m)', ytitle='time (s)'
    date_plot,title + ' log 3 decades '

; 2D electric field figures (Assume z is neglected)
endif else if (Esize(3) eq 1) then begin
    for i=ifirst,ilast,iskip do begin
        t=nout*dt*i
        label='A!E2!N: t='+ strcompress(string(t,format='(G8.4)')) $
          + ' max=' + strcompress(string(max(A2),format='(G8.2)'))
        image_plot,refrom(A2(*,*,i)), xv, yv, /aspect, $
          /zlog, nlabels=4, $
          zrange = [0,max(A2(*,*,i))], $
          title=label, $
          xtitle='x / !4k!3!D0!N', $
          ytitle='y / !4k!3!D0!N',/legend
        date_plot, title + ' log 3 decades '
endfor

; 3D electric field figures (Cross-sections)
endif else begin
    if (n_elements(A2index) ne 3) then begin
        A2max=max(A2,A2i)
        A2index=calc_index(A2,A2i)
    endif
                                ; X-Y cross-sections:	
    for i=ifirst,ilast,iskip do begin
        t=nout*dt*i
        A2sec=sqrt(reform(A2(*,*,A2index(2),i)))
        label='A!E2!N: 3D cross section x-y, z=' $
          + strcompress(string(A2index(2)*dz*nout_avg)) + ' t=' $ 
          + strcompress(string(t,format='(G8.4)')) $
          + ' max=' + strcompress(string(max(A2sec),format='(G8.2)'))
        image_plot,A2sec, xv, yv, /aspect, $
;          /zlog, nlabels=4, $
          zrange = [0,max(A2sec)], $
          title=label, $
          xtitle='x / !4k!3!D0!N', $
          ytitle='y / !4k!3!D0!N',/legend
        date_plot, title; + ' log 3 decades '

                                ; X-Z cross-section
        A2sec=sqrt(reform(A2(*,A2index(1),*,i)))
        label='A!E2!N: 3D cross section x-z, y=' $
          + strcompress(string(A2index(1)*dy*nout_avg)) + ' t=' $ 
          + strcompress(string(t,format='(G8.4)')) $
          + ' max=' + strcompress(string(max(A2sec),format='(G8.2)'))
        image_plot,A2sec, xv, zv, /aspect, $
;          /zlog, nlabels=4, $
          zrange = [0,max(A2sec)], $
          title=label, $
          xtitle='x / !4k!3!D0!N', $
          ytitle='y / !4k!3!D0!N',/legend
        date_plot, title; + ' log 3 decades '
        
                                ; Y-Z cross-section
        A2sec=sqrt(reform(A2(A2index(0),*,*,i)))
        label='A!E2!N: 3D cross section y-z, x=' $
          + strcompress(string(A2index(0)*dx*nout_avg)) + ' t=' $ 
          + strcompress(string(t,format='(G8.4)')) $
          + ' max=' + strcompress(string(max(A2sec),format='(G8.2)'))
        image_plot, A2sec, yv, zv, /aspect, $
;          /zlog, nlabels=4, $
          zrange = [0,max(A2sec)], $
          title=label, $
          xtitle='x / !4k!3!D0!N', $
          ytitle='y / !4k!3!D0!N',/legend
        date_plot, title; + ' log 3 decades '

 endfor

endelse
end
