PRO writearray, FileName, arr, order=orderarray, $
                    binary=binary, $
                    first=first, last=last, skip= skip, lineskip=lineskip, $
                    reorder_bytes=reorder_bytes, $
                    no_reorder_bytes=no_reorder_bytes, infileunit=infileunit
;+
; 
; NAME: 
;
; NAME:
;	WRITEARRAY
;
; PURPOSE:
;       To write an array as a string of numbers, either in asc or binary
;
; CATEGORY:
;	General I/O.
;
; CALLING SEQUENCE:
;	WRITEARRAY, 'FileName', arr
;
; INPUTS:
;	FileName:   The name of the file containing data
;       arr:        The array that will be written
;       orderarray: Array indicating order of data:
;                   default order=[0,1,2,...] but commonly reversed
;                   order=[...,2,1,0]  
;       binary:     Switch to indicate the data is in binary format
;       double:     Switch to indicate the data is double precision
;       first:      Integer indicating the element number (not byte) at
;                   which to begin storing data starting with 0
;       last:       Integer indicating the element number (not byte) at
;                   which to stop storing data.
;       skip:       The number of arrays (not elements or bytes) to 
;                   alternately skip
;	lineskip:   The number of lines to skip at the beginning
;       reorder_bytes:  Convert from big endian to little or visa
;                       versa - routine attempts to do this automatically
;       no_reorder_bytes:  Prevent Convertion from big endian to little or visa
;                          versa
;       fileunit    Use fileunit instead of filename (useful for
;                   reading a segment). User must close file.
;
; OUTPUTS:
;        Array: An array of dimension [nx,ny,nz,...,nt]
;-

;Read ntimes

;on_error,2                      ;Return to caller if an error occurs

; Reverse the byte ordering if the keyword is set or if the data looks
; scrambled

array = arr

on_ioerror, file_error
;Open the file!
if (n_elements(infileunit) eq 0 ) then infileunit=0

if infileunit le 0 then begin
    fileunit=11
    openw, fileunit,  FileName
endif else fileunit=infileunit

; Calculate dimensions and sizes
ndim_dim=(size(array))[0]
ndim=(size(array))[1:ndim_dim]
ndim_product=long64(ndim(0))
for i = 1, ndim_dim-1 do ndim_product=ndim_product*ndim(i)

if  not keyword_set(no_reorder_bytes) then begin
    if (keyword_set(reorder_bytes)) then byteorder, array, /ftoxdr $
    else if  (max(array) gt 1e18  or  min(array) lt -1e18) then begin
        byteorder, array, /ftoxdr
                                ; test to see if the numbers are better, if not revert
        if (max(array) gt 1e18 or min(array) lt -1e18) then $
          byteorder, array, /ftoxdr $
        else reorder_bytes=1
    endif
endif


if n_elements(orderarray) ne 0 then begin
; Array needs reordering
    ndim_reorder=ndim
    for i=0,ndim_dim-1 do ndim_reorder(i)=ndim(orderarray(i))
    array=transpose(array,[orderarray])
    array=reform(array,[ndim_reorder],/overwrite)
endif

; Determine the file size in bytes
fsize_bytes=long64((fstat(fileunit)).size)
word_size=4
if (keyword_set(double)) then word_size=8

if n_elements(first) eq 0 then first=0

if keyword_set(binary) then tlast=fsize_bytes/word_size else tlast=fsize_bytes
if n_elements(last) eq 0 then last=tlast
if last lt 0 or last gt tlast then last=tlast

if n_elements(skip) eq 0 then skip=1

req_size=(last-first)/skip


;Calculate the size of the added last dimension:
nlast=req_size/ndim_product


; Skip to first element of the file:
if (first ne 0) then begin
   bytes=first*word_size
   array=bytarr(bytes)
   for i=0, first*word_size/bytes-1 do $
      if (not keyword_set(binary)) then readf, fileunit, array $
      else readu, fileunit, array
endif

dummy=''
if (not keyword_set(binary) and n_elements(lineskip) gt 0) then $
    for i=0,lineskip-1 do readf, fileunit, dummy

if keyword_set(binary) and skip eq 1 then begin
    writeu, fileunit, array
endif else begin
   ;on_ioerror, file_end
   if (not keyword_set(binary)) then  printf, fileunit, array else writeu, fileunit, array
endelse

if (infileunit le 0) then close,fileunit
return
file_error: 
print," Error opening or writing to file: ",FileName, " on unit", $
  fileunit, " Returning ..."
if (infileunit eq 0) then close, fileunit

end
