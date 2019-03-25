function gap_remove, xxx
; Removes any small segments in the last dim of xxx
;
small=1E-6
xsize=size(xxx)
yyy=xxx

xmaxdim=xsize(xsize(0))

; Get the sum square of the content of xxx
totalx=total(xxx^2,1)
for id=1,xsize(0)-2 do totalx=total(totalx,1)

max_total=max(totalx)
count=0
for i=0, xmaxdim - 1  do begin
    if (Totalx(i) < max_total*small^2) then begin
        if xsize(0) eq 2 then yyy(*,count)=xxx(*,i) else $
          if xsize(0) eq 3 then yyy(*,*,count)=xxx(*,*,i) else $
          if xsize(0) eq 4 then yyy(*,*,*,count)=xxx(*,*,*,i)
        count += 1
    endif
endfor

if xsize(0) eq 2 then yyy=yyy(*,0:count-1) else $
  if xsize(0) eq 3 then yyy=yyy(*,*,0:count-1) else $
  if xsize(0) eq 4 then yyy=yyy(*,*,*,0:count-1)

return, yyy
end
    
