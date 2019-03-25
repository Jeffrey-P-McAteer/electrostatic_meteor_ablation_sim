FUNCTION readcomp, FileName, nelem, ntimes
;This function reads in the contents of a compressed string of numbers


if n_elements(nelem) eq 0 then begin
	print, 'Error: you must specify the length of each array'
	return, 0
end
	

; Uncompress the file
spawn,'/bin/rm -f tmp.readcomp.tmp'
fstr='zcat ' + FileName + ' > tmp.readcomp.tmp'
spawn, fstr

;Read ntimes
if n_elements(ntimes) eq 0 then ntimes=2.^22/nelem;
data=fltarr(nelem,ntimes)
vec=fltarr(nelem)
nt=0

openr, 11,  'tmp.readcomp.tmp'
normend = 0
on_ioerror, errjmp

errjmp : if (not EOF(11)) then $
  print, ' Warning: Found problem in file at line',nt,' ... Continuing'

while (not EOF(11)) and (nt lt ntimes) do begin
	readf,11, vec
        data(*,nt)=vec
	nt=nt+1
endwhile


close,11
spawn,'/bin/rm -f tmp.readcomp.tmp &'

if (nt gt 0) then begin
  d=data(*,0:nt-1) 
  return,d 
endif else begin
  print, ' Error: No data found in file: ', FileName, ' in dir:' 
  pwd
endelse

end
