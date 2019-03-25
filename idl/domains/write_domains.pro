pro write_domains, FileName, array, ndomains,$
                       fileunits=fileunits, basepath=basepath, _EXTRA = e
;+
; NAME:
;	WRITE_DOMAINS
;
; PURPOSE:
;       To write in binary data to a series of subdirectories by
;       spliting a single array into multiple files.
;
; CATEGORY:
;	Eppic I/O
;
; CALLING SEQUENCE:
;	WRITE_DOMAIN, 'FileName', [nx,ny,nz,...], ndomains
;
; INPUTS:
;	FileName:   The name of the file containing data
;       array:      The array to be split up and written
;       ndomains:   The number of domains
;       fileunits:  a list of fileunits, defaults to a list matching
;                   the number of domains
;       basepath:   Where to put the files, defaults to 'domain*'
;       OTHER:      Passed to WRITEARRAY routine
;
; OUTPUTS:
;        Array: An array of dimension [nx*nsubdomains,ny,nz,...,nt]
;
; Determine the number of subdomains with this file:  Must be labeled
; domain000/filename, domain001/filename, domain002/filename, ...
;-


if n_elements(basepath) eq 0 then basepath="domain*/"

if (n_elements(fileunits) eq 0) then fileunits=intarr(ndomains)-1 ; Make fileunits an array of -1

nsize = size(array)

ndim = nsize[0]
otherdims = nsize[2:ndim]
step = nsize[1]/ndomains
full_unfltrd = [strsplit(basepath+filename,'\*',/EXTRACT)]
nfull_unfltrd = (size(full_unfltrd))[1]
for ifile=0, ndomains-1 do begin
    if (ifile lt 10) then num_shift='00'
    if ((ifile lt 100) && (ifile gt 9)) then num_shift='0'
    if (ifile gt 99) then num_shift=''
    num_shift = num_shift + String(ifile)
    full_filename=full_unfltrd[0]
    print,"Making ",strcompress(full_filename+num_shift,/remove)
    file_mkdir,strcompress(full_filename+num_shift,/remove)
    for itxt=1,nfull_unfltrd-1 do begin
        full_filename = full_filename+num_shift+full_unfltrd[itxt]
    endfor
    full_filename=strcompress(full_filename,/remove)
    startpt = ifile*step
    endpt = (ifile+1)*step-1
    if (ndim eq 2) then begin 
        tmparray = array[startpt:endpt,*]
        ;print,"max tmparray: ",max(tmparray)
        ;print,"Tmparray ",ifile
        ;print,tmparray
        writearray,full_filename,tmparray,_EXTRA = e 
    endif
    if (ndim eq 3) then begin
        tmparray = array[startpt:endpt,*,*]
        ;print,"Tmparray ",ifile
        ;print,tmparray
        writearray,full_filename,tmparray,_EXTRA = e
    endif
endfor

end
