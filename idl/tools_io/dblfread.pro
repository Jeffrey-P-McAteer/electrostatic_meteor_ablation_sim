FUNCTION dblfread, FileName

openr, 11, FileName

nt=0
while not EOF(11) do begin
	;Determine the number of elements in the fortran record:
	np=long(0)
	readu, 11, np
	np=np/8

	on_ioerror, ioerrjmp
	;Read the elements in the fortran record
	vec=dblarr(np)
	readu, 11, vec

	;Read the dumb int at the end of the record
	idum=long(0)
	readu, 11, idum

	;Copy array into xp 
	if nt eq 0 then begin
	  nsize=2.^24/(np*8.)
	  xp=dblarr(np,nsize)
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

return, xp

END
