pro compress_domain_data, InFile, OutFile, compress=compress, $
                               _EXTRA = e
; 
; NAME:
;	COMPRESS_DOMAIN_DATA
;
; PURPOSE:
;       To read in binary data from a file, FFT the data, eliminate
;       high frequency modes and write it out again in a smaller format
;
; CATEGORY:
;	Eppic I/O
;
; CALLING SEQUENCE:
;	data_compress, InFile, OutFile, ndim,  compres=compress, $
;                       _EXTRA = e
;
; INPUTS:
;	InFile:  The name of the file containing data
;       OutFile  The name of the file used for output
;       ndim:       The dimensions of the array less one - the last
;                   dimension size may be determined from the file size
;
;       OTHER:      Passed to READARRAY routine
;
; OUTPUTS:
;       OutFile will contain the compressed data

; START CODE:
;

; Determine the number of subdomains with this file:  Must be labeled
; domain000/filename, domain001/filename, domain002/filename, ...

@eppic.i
@params_in.pro

array_size=nx2/nsubdomains*ny2*nz2 ; Size of a single array (on 1 domain)
if n_elements(istart) eq 0 then istart=0
if n_elements(istop) eq 0 then istop=-1

;Open output file
openw, 12, OutFile

nt=(file_info("domain000/"+Infile)).size/4L/array_size

;fileunits=intarr(nsubdomains)

for i=0, nt-1 do begin 

    array=read_domains(InFile,[nx2/nsubdomains,ny2,nz2], ndomains=nsubdomains, $
                       order=[2,1,0],/binary,skip=iskip,$
                       first=long64(i*1.0*array_size), $
                       last=long64(array_size*1.0*(i+1)))

    if (n_elements(array) le 1) then break
    if (NOT finite(total(array))) then stop,"Error: array not finite"

                                ;FFT Array
    af=shift(reform(fft(array),nx2,ny2,nz2),nx2/2,ny2/2,nz2/2)
                                ;Eliminate high frequency components
    ix0=nx2/2
    iy0=ny2/2
    iz0=nz2/2
    af=reform(af(ix0-ix0/compress:ix0+ix0/compress-1, $
                 iy0-iy0/compress:max([iy0+iy0/compress-1,0]), $
                 iz0-iz0/compress:max([iz0+iz0/compress-1,0]) ), $
              nx2/compress, max([ny2/compress,1]), max([nz2/compress,1]))
                                ;INVFFT Array:

    af=shift(af,[nx2,ny2,nz2]/(-2*compress))
    af=fft(af,/overwrite,/inverse)
    if (NOT finite(total(af))) then stop,"Error: array not finite"
    writeu, 12, float(transpose(af))
    print, i
endfor

close, /all

end
