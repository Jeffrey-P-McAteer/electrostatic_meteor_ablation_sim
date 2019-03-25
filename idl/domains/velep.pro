; input the potential in 1 - 3 -D .  Also set the x,y and z scales
@eppic.i
@params_in.pro
if (n_elements(velx0) lt 2) then $
  velx0=readarray('velx0.bin',[nx2,ny2,nz2],order=[2,1,0],/binary,skip=iskip)

if (ndist gt 1 and n_elements(velx1) lt 2) then $
    velx1=readarray('velx1.bin',[nx2,ny2,nz2],order=[2,1,0],/binary,skip=iskip)

if (ndist gt 2 and n_elements(velx2) lt 2) then $
    velx2=readarray('velx2.bin',[nx2,ny2,nz2],order=[2,1,0],/binary,skip=iskip)

if (ndist gt 3  and n_elements(velx3) lt 2) then $
    velx3=readarray('velx3.bin',[nx2,ny2,nz2],order=[2,1,0],/binary,skip=iskip)

if (n_elements(vely0) lt 2) then $
  vely0=readarray('vely0.bin',[nx2,ny2,nz2],order=[2,1,0],/binary,skip=iskip)

if (ndist gt 1 and n_elements(vely1) lt 2) then $
    vely1=readarray('vely1.bin',[nx2,ny2,nz2],order=[2,1,0],/binary,skip=iskip)

if (ndist gt 2 and n_elements(vely2) lt 2) then $
    vely2=readarray('vely2.bin',[nx2,ny2,nz2],order=[2,1,0],/binary,skip=iskip)

if (ndist gt 3  and n_elements(vely3) lt 2) then $
    vely3=readarray('vely3.bin',[nx2,ny2,nz2],order=[2,1,0],/binary,skip=iskip)

if (n_elements(velz0) lt 2) then $
  velz0=readarray('velz0.bin',[nx2,ny2,nz2],order=[2,1,0],/binary,skip=iskip)

if (ndist gt 1 and n_elements(velz1) lt 2) then $
    velz1=readarray('velz1.bin',[nx2,ny2,nz2],order=[2,1,0],/binary,skip=iskip)

if (ndist gt 2 and n_elements(velz2) lt 2) then $
    velz2=readarray('velz2.bin',[nx2,ny2,nz2],order=[2,1,0],/binary,skip=iskip)

if (ndist gt 3  and n_elements(velz3) lt 2) then $
    velz3=readarray('velz3.bin',[nx2,ny2,nz2],order=[2,1,0],/binary,skip=iskip)

end

