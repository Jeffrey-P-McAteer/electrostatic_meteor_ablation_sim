function shift_nowrap, matrix, index
;
; SHIFT_NOWRAP shifts a matrix but does not wrap the edges, 
;              instead placing zeroes on the shifted edges
;
s=shift(matrix,index)
; Eliminate the wrapped components
ssize=size(matrix,/structure)
for i=0,ssize.n_dimemsions-1 do begin
    mask=lint(s) mod 
