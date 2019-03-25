; calculate the index of a matrix's maximum value
function max_index, matrix
; idl matrices evolve as do fortran matrices:
max=max(matrix,index)
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

