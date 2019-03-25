PRO combine_many_domains, FileName, outdir, combine=ncombine
; 
; NAME:
;	COMBINE_DOMAINS
;
; PURPOSE:
;       This routine combines ncombine domains into 1. 
;
; CATEGORY:
;	Eppic I/O
;
; CALLING SEQUENCE:
;	COMBINE_DOMAIN, 'FileName' 
;
; INPUTS:
;	FileName:   The name of the file containing data
;       Outdir:     The name of the output directory
;       ncombine:   The number of domains to combine
;
; OUTPUTS:

; Determine the number of subdomains with this file:  Must be labeled
; domain000/filename, domain001/filename, domain002/filename, ...


@eppic.i
@params_in.pro

if n_elements(ncombine) eq 0 then ncombine=2
array_size=nx2/nsubdomains*ny2*nz2 ; Size of a single array (on 1 domain)

filename_search_strg = "domain*/" + filename
filename_list=file_search(filename_search_strg)

if n_elements(ndomains) eq 0 then ndomains=n_elements(filename_list)

ndim = [nz2,ny2,nx2/nsubdomains]

file_mkdir,outdir

for ifile=0, ndomains-1, ncombine do begin

    print, "Starting Array:", filename_list[ifile]
    array0=readarray(filename_list[ifile],ndim,/binary)
    if (size(array0))[0] gt 3 then size_last=(size(array0))[4] else size_last=1
    size_array=n_elements(array0)/size_last
    if size_array le 1 then message, 'Error: array sizes or 1'
    array0=reform(array0,[ndim,size_last])
    array_combine=array0

    for ifile2=1,ncombine-1 do begin
        print, "Adding:", filename_list[ifile+ifile2]
        array0=readarray(filename_list[ifile+ifile2],ndim,/binary)
        if (size(array0))[0] gt 3 then size_last=(size(array0))[4] else size_last=1
        size_array0=n_elements(array0)/size_last
        
        if size_array0 ne size_array then message, 'Error: array sizes differ in different domains'

        array0=reform(array0,[ndim,size_last])

        array_combine=float([[[array_combine]],[[array0]]])
    endfor
 
    cur_outdir=outdir+'/domain'+string(ifile/ncombine,format='(I3.3)')
    file_mkdir, cur_outdir

    if file_test(cur_outdir+'/'+filename) then $
      message, 'Error: File already exists, must delete manually'
    print, "Outing ", cur_outdir+'/'+filename
    openw,15, cur_outdir+'/'+filename
    writeu,15,array_combine
    close,15

endfor

;Make a modified eppic.i
new_eppic=outdir + '/' + 'eppic.i'
file_copy, 'eppic.i', new_eppic,/overwrite
openu, 44, new_eppic,  /append
printf, 44, ';Lines added to allow eppic idl processing'
printf, 44, 'nx = ' + string(nx*ncombine)
printf, 44, 'nsubdomains = '+ string(nsubdomains/ncombine)
print, 'NOTE: ADDING lines to ' + new_eppic
close, 44

;return,ifile
end

