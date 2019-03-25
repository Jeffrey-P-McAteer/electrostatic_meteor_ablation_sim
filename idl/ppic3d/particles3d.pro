; This routine brings the particle data into memory.
@params_in.pro

np=lonarr(ndist)
np(0)=npd0/npout
if (ndist gt 1) then np(1)=npd1/npout
if (ndist gt 2) then np(2)=npd2/npout
if (ndist gt 3) then np(3)=npd3/npout
if (ndist gt 4) then np(4)=npd4/npout
if (ndist gt 5) then np(5)=npd5/npout



  ntf=0
  i=0
  vx0=read2d('vx0',np(i),ntf)
  ntf=n_elements(vx0)/np(i)
  x0=read2d('x0',np(i),ntf)
  vy0=read2d('vy0',np(i),ntf)
  y0=read2d('y0',np(i),ntf)
  vz0=read2d('vz0',np(i),ntf)
  z0=read2d('z0',np(i),ntf)
  tvp=findgen( (size(vx0))[2] )*dt*npout


  if (ndist gt 1) then begin
    ntf=0
    i=1
    vx1=read2d('vx1',np(i),ntf)
    ntf=n_elements(vx1)/np(i)
    x1=read2d('x1',np(i),ntf)
    vy1=read2d('vy1',np(i),ntf)
    y1=read2d('y1',np(i),ntf)
    vz1=read2d('vz1',np(i),ntf)
    z1=read2d('z1',np(i),ntf)
  endif

  if (ndist gt 2) then begin
    ntf=0
    i=2
    vx2=read2d('vx2',np(i),ntf)
    ntf=n_elements(vx2)/np(i)
    x2=read2d('x2',np(i),ntf)
    vy2=read2d('vy2',np(i),ntf)
    y2=read2d('y2',np(i),ntf)
    vz2=read2d('vz2',np(i),ntf)
    z2=read2d('z2',np(i),ntf)
  endif

  if (ndist gt 3) then begin
    ntf=0
    i=3
    vx3=read2d('vx3',np(i),ntf)
    ntf=n_elements(vx3)/np(i)
    x3=read2d('x3',np(i),ntf)
    vy3=read2d('vy3',np(i),ntf)
    y3=read2d('y3',np(i),ntf)
    vz3=read2d('vz3',np(i),ntf)
    z3=read2d('z3',np(i),ntf)
  endif

  if (ndist gt 4) then begin
    ntf=0
    i=4
    vx4=read2d('vx4',np(i),ntf)
    ntf=n_elements(vx4)/np(i)
    x4=read2d('x4',np(i),ntf)
    vy4=read2d('vy4',np(i),ntf)
    y4=read2d('y4',np(i),ntf)
    vz4=read2d('vz4',np(i),ntf)
    z4=read2d('z4',np(i),ntf)
  endif





end
