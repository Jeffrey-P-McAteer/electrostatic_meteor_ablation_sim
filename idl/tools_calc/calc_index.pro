; calculate the index of a variable given the 1d index as returned by
; max
function calc_index, matrix, index
; idl matrices evolve as do fortran matrices:
s=size(matrix)
iv=lonarr(s(0))
nelem=s(s(0)+2)
ind2=index
for i=s(0)-1,0,-1 do begin
    nelem=nelem/s(i+1)
    iv(i)=ind2/nelem
    ind2=ind2-iv(i)*nelem
endfor
return,iv
end

