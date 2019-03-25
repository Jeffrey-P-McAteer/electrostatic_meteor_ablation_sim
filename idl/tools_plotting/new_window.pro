pro new_window , name=fname, reset=reset, no_psfont = no_psfont,  _EXTRA = e
;This IDL proceedure opens a new window or resets to the first window
;for either X or PS mode 

common window_number, winnum

if n_elements(winnum) eq 0 or keyword_set(reset) then winnum=0 $
else winnum=winnum+1

if (!D.NAME eq 'X') then window,winnum, _EXTRA = e $

else if (!D.NAME eq 'PS') then begin

    ;Specify a default name:
    if not keyword_set(fname) then begin
        fname='idl'+strcompress(string(winnum),/remove_all)+'.ps'
        print, 'Creating PostScript file: ',fname
    endif

    if (keyword_set(no_psfont)) then !p.font=-1 else !p.font=0

    ;Open the postscript device:
    device,/close
    set_plot,'ps'
    device, filename=fname, $
      /portrait, $
      xoffset=1.905,$
      xsize=17.78, $
      yoffset=12.7, $
      ysize=12.7,$
      /times, $
      _EXTRA = e 
        
endif	

end

