function read_ft, FileName_str, ndim, _EXTRA = e
; 
; NAME:
;	READ_FT
;
; PURPOSE:
;       To read in binary data from a series of subdirectories and
;       combines it into a single output array.  The data will be
;       organized in a specific fasion:
;         ikx[0], iky[0], ..., data(ikx[0], iky[0], ...)
;         ikx[1], iky[1], ..., data(ikx[1], iky[1], ...)
;       where there are ndim indicies per data component
;
; CATEGORY:
;	Eppic I/O
;
; CALLING SEQUENCE:
;	READ_FT, 'FileName', 2
;
; INPUTS:
;	FileName_str:   The name of the file(s) containing data. Can
;                       contain wild cards. If multiple files, read_FT
;                       will append data onto bottom of output array
;       ndim:       The dimension of the data and the number of
;                   indicies per data component
;
; OUTPUTS:
;        Array: An array of structures containing many ik[ndim],data


; START:
; Determine the number of subdomains with this file:  Must be labeled
; domain000/filename*, domain001/filename*, domain002/filename*, ...
filename_list=file_search(filename_str)

ndomains=n_elements(filename_list)

if ndim eq 1 then FT={ikx:0, fk:Complex(0.0,0.0)}
if ndim eq 2 then FT={ikx:0, iky:0, fk:Complex(0.0,0.0)}
if ndim eq 3 then FT={ikx:0, iky:0, ikz:0, fk:Complex(0.0,0.0)}

; Assume our indicies are 2 Byte integers and our data is a 4 byte float
int_size=2
complex_size=8

for ifile=0, ndomains-1 do begin
; Determine how big the data set is from the file size:
    file_characteristics=file_info(filename_list[ifile])
    file_size=file_characteristics.size
    
    comp_length=ndim*int_size+complex_size
    n_comp=file_size/comp_length

    if (n_comp gt 0) then begin

        fstruct_tmp=replicate(FT,n_comp)

        openr,22,filename_list[ifile]
        readu, 22, fstruct_tmp

        ; Fix the stupid endian problem
        if (max(real_part(fstruct_tmp.fk)) gt 1E18 $
            or min(real_part(fstruct_tmp.fk)) lt -1E18 or  $
              Finite(max(real_part(fstruct_tmp.fk))) eq 0 $
            or Finite(min(real_part(fstruct_tmp.fk))) eq 0 ) then begin

            swap_endian_inplace, fstruct_tmp
            if (max(real_part(fstruct_tmp.fk)) gt 1E18 $
                or min(real_part(fstruct_tmp.fk)) lt -1E18 ) then $
              swap_endian_inplace, fstruct_tmp
        endif
              
        close, 22
    
; append the new array to the end of existing
        if n_elements(fstruct) eq 0 then fstruct=fstruct_tmp else $ 
          fstruct=[fstruct,fstruct_tmp]
    endif

endfor

if (n_elements(fstruct.fk) eq 0) then begin
    msg = 'No data found in files with name: ' + filename_str + ', returning'
    message, msg
endif

return, fstruct

end
