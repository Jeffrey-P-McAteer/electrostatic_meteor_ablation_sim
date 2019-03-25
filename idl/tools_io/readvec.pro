FUNCTION readvec, FileName, nelem, ntimes, DOUBLE=double, binary=binary, skip= skip, maxsize=maxsize, reorder_bytes=reorder_bytes
;This function reads in the contents of a compressed string of numbers

if (n_elements(nelem) eq 0) then nelem=long(1)
if (nelem eq 0) then nelem=long(1)
	
;Read ntimes
if (n_elements(maxsize) eq 0) then begin
	if (!version.arch eq 'alpha') then maxsize=long(long(2)^27) $
	else maxsize=long(long(2)^27)
endif	
	
maxsize=maxsize/nelem
if (n_elements(ntimes) eq 0) then ntimes=maxsize
if (ntimes eq 0) then ntimes=maxsize
if keyword_set(double) then begin
    data=dblarr(nelem,ntimes)
    vec=dblarr(nelem)
endif else begin
    data=fltarr(nelem,ntimes)
    vec=fltarr(nelem)
endelse

if not keyword_set(skip) then skip=1

nt=long(0)

on_ioerror, openerr
openr, 11,  FileName

normend = 0
on_ioerror, errjmp

errjmp : if (not EOF(11) and nt ne 0) then $
  print, ' Warning: Found problem in file ... Continuing'

i=long(0)
while (not EOF(11)) and (nt lt ntimes) do begin
    if (not keyword_set(binary)) then readf, 11, vec $
    else begin
	readu, 11, vec
	; Reverse the byte ordering if the keyword is set or if the data looks
	; scrambled and this is an alpha
	if (keyword_set(reorder_bytes)) then byteorder, vec, /ftoxdr $
	else if (!version.arch eq 'alpha') then $
	  if (not (max(vec) lt 1e18 and min(vec) gt -1e18)) then begin
	    byteorder, vec, /ftoxdr
	    ; test to see if the numbers are better
	    if (max(vec) gt 1e30 or min(vec) lt -1e30) then $
		byteorder, vec, /ftoxdr $
	    else reorder_bytes=1
   	  endif
    endelse
    if (i mod skip eq 0 ) then begin
        data(*,nt)=vec
        nt=nt+1
    endif
    i = i + 1
endwhile

if (not EOF(11)) then print, ' Warning: Complete file not read (file too big)'
close,11

if (nt gt 0) then begin
  return, data(*,0:nt-1)
endif else begin
  print, ' Error: Insufficient data found in file: ', FileName, ' in dir:' 
  pwd
  print, ' Returning found data (Is this a binary or ASCII file)'
  return, vec
endelse

openerr: 
print," Error opening file: ",FileName, " Returning ..."
close, 11
return, 0

end
