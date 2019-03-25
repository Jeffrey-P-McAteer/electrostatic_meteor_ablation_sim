FUNCTION fortread, FileName, shortint=shortint, double=double, reorder_bytes=reorder_bytes, dims=dims
;Read a file of real*4 binary fortran
;
;dims is a vector giving the N-1 dimensions of the output vector, 
;     So if you want to read a 3D array of size 4x8x16, set dims=[4,8]
;     (and fortread will determine the size of the last dimension).

openr, 11, FileName

intlen=long(0)
if (n_elements(shortint) ne 0) then intlen=fix(0)
reallen=4
if (keyword_set(double)) then reallen=8

nt=0
while not EOF(11) do begin
	;Determine the number of elements in the fortran record:
	np=intlen
	readu, 11, np
	if (keyword_set(reorder_bytes)) then byteorder, np, /ftoxdr $
          else begin
                                ; If the file is not an integral
                                ; number times this byte size, then
                                ; it is probably endian reversed
            fsize_bytes=long((fstat(11)).size)
            if ( (fsize_bytes mod (np+reallen*2)) ne 0) then begin
                byteorder, np, /ftoxdr
                                ; test to see if the numbers are better
                if ( (fsize_bytes mod (np+reallen*2)) ne 0) then $
                  byteorder, np, /ftoxdr $
                else reorder_bytes=1
            endif
        endelse
        np=np/reallen

	on_ioerror, ioerrjmp
	;Read the elements in the fortran record
	if (keyword_set(double)) then vec=dblarr(np) $
	                         else vec=fltarr(np)
	readu, 11, vec

        ; Reverse the byte ordering if the keyword is set or if the data looks
	; scrambled and this is an alpha
	if (keyword_set(reorder_bytes)) then byteorder, vec, /ftoxdr $
	else if (!version.arch eq 'alpha') then $
	  if (max(vec) gt 1e30 or min(vec) lt -1e30) then begin
	    byteorder, vec, /ftoxdr
	    ; test to see if the numbers are better
	    if (max(vec) gt 1e30 or min(vec) lt -1e30) then $
		byteorder, vec, /ftoxdr $
	    else reorder_bytes=1
   	  endif


	;Read the dumb int at the end of the record
	idum=intlen
	readu, 11, idum

	;Copy array into xp 
	if nt eq 0 then begin
	  nsize=long(2)^24/(np*reallen)
          if (keyword_set(double)) then xp=dblarr(np,nsize) $
                                     else xp=fltarr(np,nsize)
	endif

        if nt eq nsize then begin 
            xp=[[xp],[fltarr(np,nsize)]]
            nsize=nsize*2
        endif

	xp(*,nt)=vec

	nt=nt+1
	if nt gt nsize then begin
	  print, 'Error: Array too Large'
	  stop
	endif

endwhile

ioerrjmp:
close, 11

xp=xp(*,0:nt-1)

; Reformat array to match the size dims
if (keyword_set(dims)) then begin
    xp=reform(xp,[[dims],nt])
endif

return, xp

END
