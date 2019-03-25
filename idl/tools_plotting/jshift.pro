function jshift,img2,x,y

;allows for fractional shifting of an image!

img=img2

x1=fix(x)
x2=x1+1
if x lt 0 then x2=x1-1
dx=abs(x mod 1)

y1=fix(y)
y2=y1+1
if y lt 0 then y2=y1-1
dy=abs(y mod 1)

img=(1.-dx)*shift(img,x1,0)+dx*shift(img,x2,0)
img=(1.-dy)*shift(img,0,y1)+dy*shift(img,0,y2)

return,img

end
