FUNCTION readarray, FileName, ndim, order=orderarray, $
                    DOUBLE=double, binary=binary, $
                    first=first, last=last, skip= skip, lineskip=lineskip, $
                    maxsize=maxsize, reorder_bytes=reorder_bytes, $
                    no_reorder_bytes=no_reorder_bytes, infileunit=infileunit
; 
; NAME: 
;
; NAME:
;	READARRAY
;
; PURPOSE:
;       To read in a string of numbers, either in asc or binary, and
;       to output a formatted array.
;
; CATEGORY:
;	General I/O.
;
; CALLING SEQUENCE:
;	READARRAY, 'FileName', [nx,ny,nz,...]
;
; INPUTS:
;	FileName:   The name of the file containing data
;       ndim:       The dimensions of the array less one - the last
;                   dimension size may be determined from the file size
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
;       maxsize:    Limit in bytes on the array size
;       reorder_bytes:  Convert from big endian to little or visa
;                       versa - routine attempts to do this automatically
;       no_reorder_bytes:  Prevent Convertion from big endian to little or visa
;                          versa
;       fileunit    Use fileunit instead of filename (useful for
;                   reading a segment). User must close file.
;
; OUTPUTS:
;        Array: An array of dimension [nx,ny,nz,...,nt]

;Read ntimes
if (n_elements(maxsize) eq 0) then begin
    maxsize=long64(long64(2)^34)
endif

;on_error,2                      ;Return to caller if an error occurs

on_ioerror, file_error
;Open the file!
if (n_elements(infileunit) eq 0 ) then infileunit=0

if infileunit le 0 then begin
    fileunit=11
    openr, fileunit,  FileName
endif else fileunit=infileunit

; Calculate dimensions and sizes
ndim_dim=size(ndim,/n_elements)
if n_elements(ndim) eq 0 then ndim=0
ndim_product=long64(ndim(0))
for i = 1, ndim_dim-1 do ndim_product=ndim_product*ndim(i)

; Determine the file size in bytes
fsize_bytes=long64((fstat(fileunit)).size)
if (fsize_bytes eq 0) then $
  message, string('File, ', FileName, ' empty or nonexistant'),/ioerror
word_size=4
if (keyword_set(double)) then word_size=8

if long64(ndim_product*long64(word_size)) gt maxsize then $
  message, string('Size of requested array ', ndim_product*word_size, $
                  ' exceeds maxsize memory limit ', maxsize)


if n_elements(first) eq 0 then first=0

if keyword_set(binary) then tlast=fsize_bytes/word_size else tlast=fsize_bytes

if n_elements(last) eq 0 then last=tlast
if last lt 0 or last gt tlast then last=tlast

if n_elements(skip) eq 0 then skip=1

req_size=(last-first)/skip

if req_size*word_size gt maxsize then begin
    last=maxsize/word_size*skip+first
    print, 'Size of array ', req_size*word_size, $
  ' exceeds maxsize limit ', maxsize, ' truncating at element', last
    req_size=(last-first)/skip
endif

;Calculate the size of the added last dimension:
nlast=req_size/ndim_product > 1

; Skip to first element of the file:
if (first ne 0) then begin
    if (not keyword_set(binary)) then begin
        if (first < ndim_product*nlast or first*word_size < maxsize/2) then $
          bytes=first*word_size else bytes=maxsize/2
        array=bytarr(bytes)
        for i=0, first*word_size/bytes-1 do $
          readf, fileunit, array
        if first*word_size mod bytes ne 0 then begin
            array = bytarr(first*word_size mod bytes)
            readf, fileunit, array
        endif
    endif else $
      point_lun,fileunit,long64(first*1.0*word_size)
endif

dummy=''
if (not keyword_set(binary) and n_elements(lineskip) gt 0) then $
    for i=0,lineskip-1 do readf, fileunit, dummy

array=make_array(dimension=[ndim,nlast],double=double)


if keyword_set(binary) and skip eq 1 then begin
    readu, fileunit, array
endif else begin
                                ; We will have to read the data in sections
    on_ioerror, file_end
    tmp=make_array(dimension=[ndim],double=double) 
    for count=0L,nlast-1 do begin
        if (not keyword_set(binary)) then  readf, fileunit, tmp else readu, fileunit, tmp
        array(count*ndim_product:(count+1)*ndim_product-1)=tmp
        if (not keyword_set(binary)) then $
          for i=1,skip-1 do readf, fileunit, tmp $
        else $
          point_lun,fileunit,long64((ndim_product*1L*skip*count+first)*word_size)
    endfor
endelse

file_end: 
;Shorten the array to the actual size:
if (not keyword_set(binary)) then begin
    array2=make_array(dimension=[ndim,count],double=double)
    for i=0L,count-1 do array2(i*ndim_product:(i+1)*ndim_product-1)= $
      array(i*ndim_product:(i+1)*ndim_product-1)
    array=array2
    array2=0
endif

; Reverse the byte ordering if the keyword is set or if the data looks
; scrambled

!except=0 ; turn off math error reporting
if  not keyword_set(no_reorder_bytes) then begin
    if (keyword_set(reorder_bytes)) then byteorder, array, /ftoxdr $
    else begin
        notfinite = where(finite(array)eq 0,count)
        if count gt 0 then $
          byteorder, array, /ftoxdr $
        else $
          if  (max(array,/nan) gt 1e18  or  $
               min(array,/nan) lt -1e18) then $
          byteorder, array, /ftoxdr 
                                ; test to see if the numbers are better, if not revert
         if (max(array,/nan) gt 1e18 or min(array,/nan) lt -1e18) then $
           byteorder, array, /ftoxdr $
         else reorder_bytes=1

    endelse
endif
;; clear math errors and turn reporting back on 
accum_errs = check_math()
!except=1

if (infileunit le 0) then close,fileunit

if n_elements(orderarray) ne 0 then begin
; Array needs reordering
    ndim_reorder=ndim
    for i=0,ndim_dim-1 do ndim_reorder(i)=ndim(orderarray(i))
    array=reform(array,[ndim_reorder,nlast],/overwrite)
    array=transpose(array,[orderarray,ndim_dim])
endif

return, array

file_error: 
print," Error opening or reading file: ",FileName, " on unit", $
  fileunit, " Returning ..."
if (infileunit eq 0) then close, fileunit
return, 0

end
