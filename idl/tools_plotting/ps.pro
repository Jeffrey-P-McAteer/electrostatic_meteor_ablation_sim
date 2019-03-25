;******************************************************************
  pro ps, fname, color=color,  no_psfont = no_psfont, _EXTRA = e
;******************************************************************

  if (n_elements(fname) eq 0) then fname='idl.ps'
  print, 'Creating PostScript file: ',fname
  set_plot,'ps'
;
; Black and White
;
  if not keyword_set(color) then begin
    print, 'Black and White'
      device, filename=fname, bits_per_pixel=8,_EXTRA = e, color=0
    loadct,0
;
; Color
;
  endif else begin
      ;print, 'Color - Warning you must now load your colormap - default ct=5'
      device, filename=fname,  bits_per_pixel=24, /color, _EXTRA = e 
  endelse
  loadct,5
;
    if (keyword_set(no_psfont)) then !p.font=-1 else !p.font=0
;
  end
;-)

