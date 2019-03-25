function read_domains, FileName, ndim,  ndomains=ndomains, fileunits=fileunits, _EXTRA = e
; 
; NAME:
;	READ_DOMAINS
;
; PURPOSE:
;       To read in binary data from a series of subdirectories and
;       combines it into a single output array.
;
; CATEGORY:
;	Eppic I/O
;
; CALLING SEQUENCE:
;	READ_DOMAIN, 'FileName', [nx,ny,nz,...]
;
; INPUTS:
;	FileName:   The name of the file containing data
;       ndim:       The dimensions of the array less one - the last
;                   dimension size may be determined from the file size
;       segment     Read a segment: do not close the file after each read
;       OTHER:      Passed to READARRAY routine
;
; OUTPUTS:
;        Array: An array of dimension [nx*nsubdomains,ny,nz,...,nt]

; Determine the number of subdomains with this file:  Must be labeled
; domain000/filename, domain001/filename, domain002/filename, ...


filename_search_strg = "domain*/" + filename
filename_list=file_search(filename_search_strg)

if n_elements(ndomains) eq 0 then ndomains=n_elements(filename_list)

if (n_elements(fileunits) eq 0) then fileunits=intarr(ndomains)-1 ; Make fileunits an array of -1

if fileunits[0] eq 0 then begin
    fileunits=indgen(ndomains)+13
    for ifile=0, ndomains-1 do openr, fileunits[ifile],  filename_list[ifile]
endif

array_tmp=readarray(filename_list[0],ndim, infileunit=fileunits[0], _EXTRA = e )
if (size(array_tmp))[0] gt 3 then size_last=(size(array_tmp))[4] else size_last=1
size_array=n_elements(array_tmp)/size_last
if size_array le 1 then return, 0

; Do all the files contain the same amount of data:
file_characteristics=file_info(filename_list)
file_sizes=file_characteristics.size
if size(array_tmp,/type) eq 5 then byte_per_element = 8 else byte_per_element = 4
if max(file_sizes) lt n_elements(array_tmp)*byte_per_element then $
    print, 'Warning (read_domains): Arrays not all same size, using minimum'

if min(file_sizes) lt n_elements(array_tmp)*byte_per_element then begin
    print, 'Warning (read_domains): Arrays not all same size, using minimum'
    size_last=min(file_sizes)/byte_per_element/(ndim[0]*ndim[1]*ndim[2])
endif

;array=reform(array[*,*,*,0:size_last-1],[ndim,size_last])
ndim2 = ndim
ndim2[0] *= ndomains
array=fltarr([ndim2,size_last])
array[0:ndim[0]-1,*,*,*]=array_tmp[*,*,*,0:size_last-1]

for ifile=1, ndomains-1 do begin
   print, 'reading domain ',ifile
    array_tmp=readarray(filename_list[ifile],ndim,infileunit=fileunits[ifile],_EXTRA = e)
    array_tmp=reform(array_tmp[*,*,*,0:size_last-1],[ndim,size_last])
    size_tmp=n_elements(array_tmp)/((size(array_tmp))[4])
; Append it to array:
    if size_tmp eq size_array and ((size(array_tmp))[4]) ge size_last then begin
        array[ifile*ndim[0]:(ifile+1)*ndim[0]-1,*,*,*]=array_tmp[*,*,*,0:size_last-1]
    endif  else begin
        print, "Warning: Array ", filename_list[ifile], " has size ", $
          n_elements(array_tmp)/((size(array_tmp))[4])
        print, "  This does not match size of initial array ", $
        n_elements(array)/((size(array))[4])
    endelse
endfor

if (NOT finite(total(array))) then stop,"Error in read_domains: array not finite"
return, array

end
