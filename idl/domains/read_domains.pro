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

;; get some info about the file
shiftl=0
shiftr=0
if (n_elements(non_periodic_x) ne 0) then shiftl = non_periodic_x[0]
if n_elements(non_periodic_x) ne 0 then shiftr=non_periodic_x[1]

ncells = long(ndim[0]+shiftr+shiftl)*long(ndim[1])*long(ndim[2])
totaltime = timesteps(Filename,ncells,ndomains,double=True)

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

ndim_tmp[0]=ndim_tmp[0]+shiftl+shiftr

;; get timing info from extra

extraTags = TAG_NAMES(e)

first = 0
whereFirst = WHERE(STRCMP(extraTags,'FIRST') EQ 1)
if (whereFirst ne -1) then first = e.(whereFirst)

last = -1*ncells
whereLast = WHERE(STRCMP(extraTags,'LAST') EQ 1)
if (whereLast ne -1) then last = e.(whereLast)
if last/ncells gt totaltime then last = totaltime*ncells
if last/ncells lt 0 then last = totaltime*ncells

skip = 1
whereSkip = WHERE(STRCMP(extraTags,'SKIP') EQ 1)
if (whereSkip ne -1) then skip = e.(whereSkip)

tsize = (last-first)/ncells/skip

;; build full array

array = findgen([shiftl+shiftr+ndim[0]*ndomains,ndim[1],ndim[2],tsize])

;; gather each domains array, and stitch together

array_tmp=readarray(filename_list[0],ndim_tmp, infileunit=fileunits[0],_EXTRA = e)

if n_elements(array) le 1 then return, 0


if (ndomains gt 1) then begin
    array[0:shiftl+ndim[0]-1,*,*,*] = array_tmp[0:shiftl+ndim[0]-1,*,*,0:tsize-1]
endif else begin
    array = array_tmp
endelse

if (ndomains gt 1) then begin
    ;; loop over all but last domain
    for ifile=1, ndomains-2 do begin
        print,"... '",filename_list[ifile],"'"
        ndim_tmp = ndim
        ndim_tmp[0]=ndim[0]+shiftl+shiftr

        ;; read in array with left and right padding
        array_tmp=readarray(filename_list[ifile],$
                            ndim_tmp,infileunit=fileunits[ifile],_EXTRA = e)

        array[shiftl+(ifile)*ndim[0]:shiftl+(ifile+1)*ndim[0]-1,*,*,*] = $
          array_tmp[shiftl:shiftl+ndim[0]-1,*,*,0:tsize-1]

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
; Append it to array:
    array[shiftl+(ifile)*ndim[0]:shiftl+(ifile+1)*ndim[0]+shiftr-1,*,*,*] = $
      array_tmp[shiftl:shiftl+ndim[0]+shiftr-1,*,*,0:tsize-1]

    
endif
return, array
    
end
