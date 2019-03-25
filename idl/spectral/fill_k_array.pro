function fill_k_array, fk_struct, ikrange=ikrange

; fill_k_array
;
; function: Fills an array from a structure list 
;
; Input:
;       fk_struct: array of structures arranged depending on dimension
;       as:
;                  1D:(ikx,x)
;                  2D:(ikx,iky,x)
;                  3D:(ikx,iky,ikz,x)
;       irange: Range of array to be inserted. Any indicies out of
;       range will be ignored. DEFAULT: max of list.
;               
ndim=n_tags(fk_struct)-1

if n_elements(ikRange) eq 0 then begin
    ikRange=intarr(ndim,2)
    for i=0,ndim-1 do ikRange[i,0]=min(fk_struct.(i))
    for i=0,ndim-1 do ikRange[i,1]=max(fk_struct.(i))
endif else begin
    ; Test ikrange to make sure it is the correct size
    message,"ikrange option not fully implemented"
    if size(ikRange) ne size(intarr(ndim,2)) then message, "ikRange in fill_k_array wrong dimensions"
endelse

fk=complexarr(ikrange(*,1)-ikrange(*,0)+1)

if ndim eq 3 then fk(fk_struct.ikx,fk_struct.iky,fk_struct.ikz)=fk_struct.fk $
else if ndim eq 2 then fk(fk_struct.ikx,fk_struct.iky)=fk_struct.fk $
else fk(fk_struct.ikx)=fk_struct.fk 

return, fk

end
