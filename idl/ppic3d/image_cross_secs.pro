pro, image_cross_secs, phi, xv, yv, zv, tv, _EXTRA_e= = e

; Generate images of phi
; Use the variables ifirst,ilast, iskip to determine how many images

;Calculate E2 in real space:
if n_elements(phi) le 2 then begin
    message, "Error: array passed to image_cross_secs is not an array"
endif
phisize=size(phi)

if n_elements(ifirst) eq 0 then ifirst=0
if n_elements(ilast) eq 0 then ilast=phisize(4)-1
if n_elements(iskip) eq 0 then iskip=1


;1D plots:
if (phisize(2) eq 1 and phisize(3) eq 1) then begin
    label='Potential:'
    image_plot,phi, xv, tv, /legend, nlabels=4, $
      title=label, xtitle='x (m)', ytitle='time (s)'
    date_plot,title 

; 2D electric field figures (Assume z is neglected)
endif else if (phisize(3) eq 1) then begin
    for i=ifirst,ilast,iskip do begin
        t=nout*dt*i
        label='Potential: t='+ strcompress(string(t,format='(G8.4)')) $
          + ' max=' + strcompress(string(max(phi(*,*,0,i)),format='(G8.2)'))
        image_plot, reform(phi(*,*,0,i)), xv, yv, /aspect, $
          nlabels=4, $
          title=label, $
          xtitle='x / !4k!3!D0!N', $
          ytitle='y / !4k!3!D0!N',/legend
        date_plot, title 
endfor

; 3D electric field figures (Cross-sections)
endif else begin
    ; X-Y cross-sections:	
    if (dx*nx/(dy*ny)>4. or dx*nx/(dy*ny)<0.25) then aspect=0 else aspect=-1
    for i=ifirst,ilast,iskip do begin
        t=nout*dt*i
        label='Potential: 3D cross section x-y t=' $ 
          + strcompress(string(t,format='(G8.4)')) $
          + ' max=' + strcompress(string(max(phi(*,*,0,i)),format='(G8.2)'))
        image_plot,reform(phi(*,*,0,i)), xv, yv, aspect=aspect, $
          nlabels=4, $
          title=label, $
          xtitle='x / !4k!3!D0!N', $
          ytitle='y / !4k!3!D0!N',/legend
        date_plot, title 
     endfor

    ; X-Z cross-section
    if (dx*nx/(dz*nz)>4. or dx*nx/(dz*nz)<0.25) then aspect=0 else aspect=-1
    for i=ifirst,ilast,iskip do begin
        t=nout*dt*i
        label='Potential: 3D cross section x-z t=' $ 
          + strcompress(string(t,format='(G8.4)')) $
          + ' max=' + strcompress(string(max(phi(*,0,*,i)),format='(G8.2)'))
        image_plot,reform(phi(*,0,*,i)), xv, zv, aspect=aspect, $
          nlabels=4, $
          title=label, $
          xtitle='x / !4k!3!D0!N', $
          ytitle='z / !4k!3!D0!N',/legend
        date_plot, title 
     endfor

     ; Y-Z cross-section
    if (dy*ny/(dz*nz)>4. or dy*ny/(dz*nz)<0.25) then aspect=0 else aspect=-1
     for i=ifirst,ilast,iskip do begin
        t=nout*dt*i
        label='Potential: 3D x-section y-z t=' $ 
          + strcompress(string(t,format='(G8.4)')) $
          + ' max=' + strcompress(string(max(phi(0,*,*,i)),format='(G8.2)'))
        image_plot,reform(phi(0,*,*,i)), yv, zv, aspect=aspect, $
          nlabels=4, $
          title=label, $
          xtitle='y / !4k!3!D0!N', $
          ytitle='z / !4k!3!D0!N',/legend
        date_plot, title 
    endfor

endelse
end
