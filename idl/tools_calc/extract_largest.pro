FUNCTION extract_largest, array, threshold_fraction

asize=size(array)
alast=asize(asize(0))

amax=max(array)*threshold_fraction

icount=0
for i=0,alast-1 do if max(array(*,*,*,i)) ge amax then icount=icount+1

array2=make_array(dimension=[asize(1),asize(2),asize(3),icount], $
                 type=asize(asize(0)+1))

icount=0
for i=0,alast-1 do begin
    if max(array(*,*,*,i)) ge amax then begin
        array2(*,*,*,icount)=array(*,*,*,i)
        icount=icount+1
    endif
endfor

return,array2
end
