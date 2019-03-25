FUNCTION fread, FileUnit
; Read a single record of binary fortran real*4 data from FileUnit
; The file should be opened with a command like: "openr, 11, FileName"

on_ioerror, errjmp

;Determine the number of elements in the fortran record:
np=long(0)
readu, FileUnit, np
np=np/4

;Read the elements in the fortran record
vec=fltarr(np)
readu, FileUnit, vec

;Read the dumb int at the end of the record
idum=long(0)
readu, FileUnit, idum

return, vec

errjmp:

return, 0

END
