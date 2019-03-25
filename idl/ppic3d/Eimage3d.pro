; Generate images of Ex, Ey, Ez
; Use the variables ifirst,ilast, iskip to determine how many images

@phi3d

Esize=size(phi)

if n_elements(ifirst) eq 0 then ifirst=0
if n_elements(ilast) eq 0 then ilast=Esize(4)-1
if n_elements(iskip) eq 0 then iskip=1
if n_elements(iskip_phi) eq 0 then iskip_phi=1
if n_elements(aspect_flag) eq 0 then aspect_flag=1

;1D plots:
if (Esize(2) eq 1 and Esize(3) eq 1) then begin
    Ex=reform(phi)
    for i =0, Esize(4) do Ex(i)=gradient(E(*,0,0,i),dx=dx*nout_avg)
    label='E!Ix!N:'
    image_plot,reform(E), xv, tv*iskip_phi, /legend, $
      title=label, xtitle='x (m)', ytitle='time (s)'
    date_plot,title

; 2D electric field figures (Assume z is neglected)
endif else if (Esize(3) eq 1) then begin
    for i=ifirst,ilast,iskip do begin
        t=nout*dt*i*iskip_phi
        E=gradient(phi(*,*,0,i),dx=dx*nout_avg,dy=dy*nout_avg)
        label='E!Ix!N: t='+ strcompress(string(t,format='(G8.4)')) $
          + ' max=' + strcompress(string(max(E(*,*,0)),format='(G8.2)'))
        image_plot, reform(E(*,*,0)), xv, yv, aspect=aspect_flag, $
          title=label, $
          xtitle='x / !4k!3!D0!N', $
          ytitle='y / !4k!3!D0!N',/legend
        date_plot, title

        label='E!Iy!N: t='+ strcompress(string(t,format='(G8.4)')) $
          + ' max=' + strcompress(string(max(E(*,*,1)),format='(G8.2)'))
        image_plot, reform(E(*,*,1)), xv, yv, aspect=aspect_flag, $
          title=label, $
          xtitle='x / !4k!3!D0!N', $
          ytitle='y / !4k!3!D0!N',/legend
        date_plot, title
endfor

; 3D electric field figures (Cross-sections)
endif else begin
    ; X-Y cross-sections:	
    for i=ifirst,ilast,iskip do begin
        t=nout*dt*i*iskip_phi
        E=gradient(phi(*,*,*,i), dx=dx*nout_avg, $
          	   dy=dy*nout_avg, dz=dz*nout_avg)

        label='E!Ix!N: 3D cross section x-y t=' $ 
          + strcompress(string(t,format='(G8.4)')) + ' max=' $
          + strcompress(string(max(E(*,*,Esize(3)/2,0)),format='(G8.2)'))
        image_plot,reform(E(*,*,Esize(3)/2,0)), xv, yv, aspect=aspect_flag, $
          title=label, $
          xtitle='x / !4k!3!D0!N', $
          ytitle='y / !4k!3!D0!N',/legend
        date_plot, title

        label='E!Iy!N: 3D cross section x-y t=' $ 
          + strcompress(string(t,format='(G8.4)')) + ' max=' $
          + strcompress(string(max(E(*,*,Esize(3)/2,1)),format='(G8.2)'))
        image_plot,reform(E(*,*,Esize(3)/2,1)), xv, yv, aspect=aspect_flag, $
          title=label, $
          xtitle='x / !4k!3!D0!N', $
          ytitle='y / !4k!3!D0!N',/legend
        date_plot, title

        label='E!Iz!N: 3D cross section x-y t=' $ 
          + strcompress(string(t,format='(G8.4)')) + ' max=' $
          + strcompress(string(max(E(*,*,Esize(3)/2,2)),format='(G8.2)'))
        image_plot,reform(E(*,*,Esize(3)/2,2)), xv, yv, aspect=aspect_flag, $
          title=label, $
          xtitle='x / !4k!3!D0!N', $
          ytitle='y / !4k!3!D0!N',/legend
        date_plot, title

     endfor

    ; X-Z cross-section
    for i=ifirst,ilast,iskip do begin
        t=nout*dt*i*iskip_phi
        E=gradient(phi(*,*,*,i), dx=dx*nout_avg, $
                          dy=dy*nout_avg, dz=dz*nout_avg)

        label='E!Ix!N: 3D cross section x-z t=' $ 
          + strcompress(string(t,format='(G8.4)')) + ' max=' + $
            strcompress(string(max(E(*,Esize(2)/2,*,0)),format='(G8.2)'))
        image_plot,reform(E(*,Esize(2)/2,*,0)), xv, zv, aspect=aspect_flag, $
          title=label, $
          xtitle='x / !4k!3!D0!N', $
          ytitle='y / !4k!3!D0!N',/legend
        date_plot, title + ' log 2 decades '

        label='E!Iy!N: 3D cross section x-z t=' $ 
          + strcompress(string(t,format='(G8.4)')) + ' max=' + $
            strcompress(string(max(E(*,Esize(2)/2,*,1)),format='(G8.2)'))
        image_plot,reform(E(*,Esize(2)/2,*,1)), xv, zv, aspect=aspect_flag, $
          title=label, $
          xtitle='x / !4k!3!D0!N', $
          ytitle='y / !4k!3!D0!N',/legend
        date_plot, title + ' log 2 decades '

        label='E!Iz!N: 3D cross section x-z t=' $ 
          + strcompress(string(t,format='(G8.4)')) + ' max=' + $
            strcompress(string(max(E(*,Esize(2)/2,*,2)),format='(G8.2)'))
        image_plot,reform(E(*,Esize(2)/2,*,2)), xv, zv, aspect=aspect_flag, $
          title=label, $
          xtitle='x / !4k!3!D0!N', $
          ytitle='y / !4k!3!D0!N',/legend
        date_plot, title + ' log 2 decades '
     endfor

     ; Y-Z cross-section
     for i=ifirst,ilast,iskip do begin
        t=nout*dt*i*iskip_phi
        E=gradient(phi(*,*,*,i), dx=dx*nout_avg, $
                          dy=dy*nout_avg, dz=dz*nout_avg)

        label='E!Ix!N: 3D cross section y-z t=' $ 
          + strcompress(string(t,format='(G8.4)')) + ' max=' $
          + strcompress(string(max(E(Esize(1)/2,*,*,0)),format='(G8.2)'))
        image_plot,reform(E(Esize(1)/2,*,*,0)), yv, zv, aspect=aspect_flag, $
          title=label, $
          xtitle='x / !4k!3!D0!N', $
          ytitle='y / !4k!3!D0!N',/legend
        date_plot, title + ' log 2 decades '

        label='E!Iy!N: 3D cross section y-z t=' $ 
          + strcompress(string(t,format='(G8.4)')) + ' max=' $
          + strcompress(string(max(E(Esize(1)/2,*,*,1)),format='(G8.2)'))
        image_plot,reform(E(Esize(1)/2,*,*,1)), yv, zv, aspect=aspect_flag, $
          title=label, $
          xtitle='x / !4k!3!D0!N', $
          ytitle='y / !4k!3!D0!N',/legend
        date_plot, title + ' log 2 decades '

        label='E!Iz!N: 3D cross section y-z t=' $ 
          + strcompress(string(t,format='(G8.4)')) + ' max=' $
          + strcompress(string(max(E(Esize(1)/2,*,*,2)),format='(G8.2)'))
        image_plot,reform(E(Esize(1)/2,*,*,2)), yv, zv, aspect=aspect_flag, $
          title=label, $
          xtitle='x / !4k!3!D0!N', $
          ytitle='y / !4k!3!D0!N',/legend
        date_plot, title + ' log 2 decades '
    endfor

endelse
end
