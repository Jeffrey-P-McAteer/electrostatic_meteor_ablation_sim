
a=readarray('log',6)
nx=32
ny=2
nz=2
index=fix([a[0,*],a[1,*]*nx/4+a[2,*],a[3,*],a[4,*]])
b=fltarr(3,nx,ny,nz)
for i=0,long(3)*nx*ny*nz-1 do b[index[0,i],index[1,i],index[2,i],index[3,i]]=a[5,i]
plot,b[0,*,0,0]
oplot,b[0,*,1,0]
oplot,b[0,*,1,1]
oplot,b[0,*,0,1]
 
end
