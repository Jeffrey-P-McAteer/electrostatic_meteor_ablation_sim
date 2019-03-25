;******************************************************************
  pro start_ps, fname, color=color, land=land, font_size=font_size, $
                no_psfont = no_psfont
;******************************************************************

  if (n_elements(fname) eq 0) then fname='idl.ps'
  print, 'Creating PostScript file: ',fname
  set_plot,'ps'
  if not keyword_set(no_psfont) then !p.font = 0
  if not keyword_set(font_size) then font_size = 12
;
; Black and White
;
  if not keyword_set(color) then begin
    print, 'Black and White'
    if keyword_set(land) then             $
      device, /landscape, filename=fname, $
             xoffset=1.27,                $ ;; 0.5 inches (in cm) from bottom
             xsize=25.4,                  $ ;;10.0 inches (in cm) width
             yoffset=26.67,               $ ;;10.5 inches (in cm) from right
             ysize=19.05,                 $ ;; 7.5 inches (in cm) height
             font_size = font_size,       $
             /times,color=0               $
    else                                  $
      device, /portrait, filename=fname,  $
             xoffset=1.27,                $ ;; 0.5 inches (in cm) from left
             xsize=19.05,                 $ ;; 7.5 inches (in cm) width
             yoffset=1.27,                $ ;; 0.5 inches (in cm) from bottom
             ysize=25.4,                  $ ;;10.0 inches (in cm) height
             font_size = font_size,       $
             /times
;
; Color
;
  endif else begin
    print, 'Color'
    if keyword_set(land) then             $
      device, /landscape, filename=fname, $
      xoffset=1.27,                $ ;; 0.5 inches (in cm) from bottom
      xsize=25.4,                  $ ;;10.0 inches (in cm) width
      yoffset=26.67,               $ ;;10.5 inches (in cm) from right
      ysize=19.05,                 $ ;; 7.5 inches (in cm) height
      font_size = font_size,       $
      /times, /color,              $
      BITS_PER_PIXEL=8             $

    else                                  $
      device, /portrait, filename=fname,  $
             xoffset=1.27,                $ ;; 0.5 inches (in cm) from left
             xsize=19.05,                 $ ;; 7.5 inches (in cm) width
             yoffset=1.27,                $ ;; 0.5 inches (in cm) from bottom
             ysize=25.4,                  $ ;;10.0 inches (in cm) height
             font_size = font_size,       $
             /times, /color,              $
      BITS_PER_PIXEL=8
  endelse
;
  end
;-)
