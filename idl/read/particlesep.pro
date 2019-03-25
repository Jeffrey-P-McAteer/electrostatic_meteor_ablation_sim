; This routine brings the particle data into memory.
@eppic.i
@params_in.pro


np=lonarr(ndist)
np(0)=npd0/npout
if (ndist gt 1) then np(1)=npd1/npout
if (ndist gt 2) then np(2)=npd2/npout
if (ndist gt 3) then np(3)=npd3/npout
if (ndist gt 4) then np(4)=npd4/npout
if (ndist gt 5) then np(5)=npd5/npout

filename='x0.bin'
filename_search_strg = "domain*/" + filename
filename_list=file_search(filename_search_strg)



x0=read_domains('x0.bin',[np(0),1,1],/binary,skip=iskip)
y0=read_domains('y0.bin',[np(0),1,1],/binary,skip=iskip)
z0=read_domains('z0.bin',[np(0),1,1],/binary,skip=iskip)

vx0=read_domains('vx0.bin',[np(0),1,1],/binary,skip=iskip)
vy0=read_domains('vy0.bin',[np(0),1,1],/binary,skip=iskip)
vz0=read_domains('vz0.bin',[np(0),1,1],/binary,skip=iskip)

if (ndist gt 1) then x1=read_domains('x1.bin',[np(1),1,1],/binary,skip=iskip)
if (ndist gt 1) then y1=read_domains('y1.bin',[np(1),1,1],/binary,skip=iskip)
if (ndist gt 1) then z1=read_domains('z1.bin',[np(1),1,1],/binary,skip=iskip)

if (ndist gt 1) then vx1=read_domains('vx1.bin',[np(1),1,1],/binary,skip=iskip)
if (ndist gt 1) then vy1=read_domains('vy1.bin',[np(1),1,1],/binary,skip=iskip)
if (ndist gt 1) then vz1=read_domains('vz1.bin',[np(1),1,1],/binary,skip=iskip)



