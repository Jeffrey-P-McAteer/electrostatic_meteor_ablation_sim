function read_domains, FileName, ndim,  ndomains=ndomains, $
                       fileunits=fileunits, basepath=basepath, $
                       non_periodic_x=non_periodic_x, _EXTRA = e

;+
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
;                   dimension size may be determined from the file
;                   size
;       ndomains:   The number of domains, defaults to number of files
;                   found 
;       fileunits:  a list of fileunits, defaults to a list matching
;                   the number of files found
;       basepath:   Where to look for the files, defaults to 'domain*'
;       OTHER:      Passed to READARRAY routine
;
; OUTPUTS:
;        Array: An array of dimension [nx*nsubdomains,ny,nz,...,nt]
;
; Determine the number of subdomains with this file:  Must be labeled
; domain000/filename, domain001/filename, domain002/filename, ...
;-

if n_elements(basepath) eq 0 then basepath="domain*/"

filename_search_strg = basepath + filename
filename_list=file_search(filename_search_strg)

if n_elements(filename_list) eq 0 then begin
    print,"No files found; seach string =  '",$
      filename_search_strg,"'; returning..."
    return, 1
endif

if (n_elements(filename_list) eq 1) && (filename_list[0] eq '') then begin
    print,"No files found; seach string = '",$
      filename_search_strg,"'; returning..."
    return, 1
endif

print,"Reading domain information from ..."
print,"... '",filename_list[0],"'"
;for i=0,n_elements(filename_list)-1 do print,"... '",filename_list[i],"'"
if n_elements(ndomains) eq 0 then ndomains=n_elements(filename_list)

; Make fileunits an array of -1
if (n_elements(fileunits) eq 0) then fileunits=intarr(ndomains)-1 

if fileunits[0] eq 0 then begin
    fileunits=indgen(ndomains)+13
    for ifile=0, ndomains-1 do openr, fileunits[ifile],  filename_list[ifile]
endif

ndim_tmp=ndim
shiftl=0
shiftr=0
if (n_elements(non_periodic_x) ne 0) then shiftl = non_periodic_x[0]
if n_elements(non_periodic_x) ne 0 then shiftr=non_periodic_x[1]

ndim_tmp[0]=ndim_tmp[0]+shiftl+shiftr


array=readarray(filename_list[0],ndim_tmp, infileunit=fileunits[0],_EXTRA = e)
if (size(array))[0] gt 3 then size_last=(size(array))[4] else size_last=1
size_array=n_elements(array)/size_last
;print,filename_list[0],(size(array))[1],(size(array))[2]
if size_array le 1 then return, 0

if (ndomains gt 1) then $
    array=reform(array[0:ndim_tmp[0]-1-shiftr,*,*,0:size_last-1],$
                 [ndim_tmp[0]-shiftr,ndim_tmp[1],ndim_tmp[2],size_last])




if (ndomains gt 1) then begin
    ;; loop over all but last domain
    for ifile=1, ndomains-2 do begin
        print,"... '",filename_list[ifile],"'"
        ndim_tmp = ndim
        ndim_tmp[0]=ndim[0]+shiftl+shiftr

        ;; read in array with left and right padding
        array_tmp=readarray(filename_list[ifile],$
                            ndim_tmp,infileunit=fileunits[ifile],_EXTRA = e)

        size_last_tmp=size_last
        if (size(array_tmp))[0] gt 3 then $
          size_last_tmp=(size(array_tmp))[4] else size_last_tmp=1
        if size_last_tmp lt size_last then size_last = size_last_tmp
        
        ;; crop left and right padding
        array_tmp=$
          reform(array_tmp[shiftl:ndim_tmp[0]-1-shiftr,*,*,0:size_last_tmp-1],$
                 [ndim,size_last_tmp]) 
          
        size_tmp=n_elements(array_tmp)/((size(array_tmp))[4])
        
        ;; Append it to array:
        if size_tmp eq size_array then begin 
            if size_last_tmp ge size_last then begin
                array=[array[*,*,*,0:size_last-1],$
                       array_tmp(*,*,*,0:size_last-1)]
            endif  else begin
                if n_elements(non_periodic_x) ne 0 then begin
                    array=[array,array_tmp(*,*,*,0:size_last-1)]
                endif else begin
                    print, "Warning: Array ", filename_list[ifile], $
                      " has time ",size_last_tmp
                    print, "  This does not match time of initial array ", $
                      size_last
                    print, " Truncating origonal array so sizes match..."
                    if ((size(array_tmp))[4] lt size_last) then $
                      size_last=(size(array_tmp))[4]
                    array=[array[*,*,*,0:size_last-1],$
                           array_tmp[*,*,*,0:size_last-1]]
                endelse
            endelse
        endif else begin
            if n_elements(non_periodic_x) ne 0 then begin
                array=[array(*,*,*,0:size_last-1),$
                       array_tmp(*,*,*,0:size_last-1)]
            endif else begin
                print,"Warning: Arrays not the same size. Current domain (",$
                  filename_list[ifile],"): ",$
                  size_tmp,"; Previous domain(s):",$
                  size_array," ... returning ..."
                return, 0
            endelse
        endelse
    endfor ;;; done with all but last domain
    ;; now last domain
    print,"... '",filename_list[ifile],"'"    
    ndimRead = ndim
    ndim_tmp=ndim

    ndim_tmp[0]=ndim_tmp[0]+shiftr
    ndimRead[0]=ndim[0]+shiftl+shiftr


    ;; read all including padding
    array_tmp=readarray(filename_list[ifile],ndimRead,$
                        infileunit=fileunits[ifile],_EXTRA = e)
    size_last_tmp=size_last
    if (size(array_tmp))[0] gt 3 then $
      size_last_tmp=(size(array_tmp))[4] else size_last_tmp=1

    ;; crop left only, and extra time 
    array_tmp=reform(array_tmp[shiftl:ndimRead[0]-1,*,*,0:size_last_tmp-1],$
                     [ndim_tmp,size_last_tmp])

    size_tmp=n_elements(array_tmp)/((size(array_tmp))[4])

; Append it to array:
    if size_tmp eq size_array then begin
        if ((size(array_tmp))[4]) ge size_last then begin
            array=[array,array_tmp(*,*,*,0:size_last-1)]
        endif  else begin
            if n_elements(non_periodic_x) ne 0 then begin
                array=[array,array_tmp(*,*,*,0:size_last-1)]
            endif else begin
                print, "Warning: Array ", filename_list[ifile], $
                  " has time ",size_last_tmp
                print, "  This does not match time of initial array ", $
                  size_last
                print, " Truncating origonal array so sizes match..."
                if ((size(array_tmp))[4] lt size_last) then $
                  size_last=(size(array_tmp))[4]
                array=[array[*,*,*,0:size_last-1],$
                       array_tmp[*,*,*,0:size_last-1]]
            endelse
        endelse
    endif else begin
        if n_elements(non_periodic_x) ne 0 then begin
            array=[array,array_tmp(*,*,*,0:size_last-1)]
        endif else begin
            print,"Warning: Arrays not the same size. Current domain: ",$
              size_tmp,"; Previous domain(s):",size_array,"... returning ..."
            return, 0
        endelse
    endelse

endif
return, array
    
end
