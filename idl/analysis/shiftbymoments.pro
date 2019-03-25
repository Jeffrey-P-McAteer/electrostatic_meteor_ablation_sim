;; get density
@denep

;; get moments
moments0 = readarray('domain000/moments0.out',13,lineskip=1,skip=iskip)
moments1 = readarray('domain000/moments1.out',13,lineskip=1,skip=iskip)

;; calculate shift
totalt0=(size(moments0))[2]
totalt1=(size(moments1))[2]

shiftx0 = findgen(totalt0)
for it=1,totalt0-1 do shiftx0[it] = $
  floor(total(moments0[1,1:it])*moments0[0,1]*dt/dx/nout_avg)
shifty0 = findgen(totalt0)
for it=1,totalt0-1 do shifty0[it] = $
  floor(total(moments0[5,1:it])*moments0[0,1]*dt/dx/nout_avg)
shiftx1 = findgen(totalt1)
for it=1,totalt1-1 do shiftx1[it] = $
  floor(total(moments1[1,1:it])*moments1[0,1]*dt/dx/nout_avg)
shifty1 = findgen(totalt1)
for it=1,totalt1-1 do shifty1[it] = $
  floor(total(moments1[5,1:it])*moments1[0,1]*dt/dx/nout_avg)

;; shift
den0shift = den0
for it=1,totalt0-1 do $
  den0shift[*,*,0,it] = shift(den0[*,*,0,it],shiftx0[it],shifty0[it])
den1shift = den1
for it=1,totalt1-1 do $
  den1shift[*,*,0,it] = shift(den1[*,*,0,it],shiftx1[it],shifty1[it])

end
