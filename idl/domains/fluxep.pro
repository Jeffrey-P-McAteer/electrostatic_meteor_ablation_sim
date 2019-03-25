; Read in fluxsities
@eppic.i
@params_in.pro
nx3=nx2/nsubdomains

if (n_elements(istart) eq 0) then istart=0
if (n_elements(iend) eq 0) then iend=-1

sizepertime=long64((nx3*1.0*ny2))
if(ndim_space eq 3) then sizepertime = long64(sizepertime*nz2)


if n_elements(reread) eq 0 then reread=1

if (reread or n_elements(fluxx0) lt 2) then $
 fluxx0=read_domains('fluxx0.bin',[nx2/nsubdomains,ny2,nz2], ndomains=nsubdomains, $
                   order=[2,1,0],/binary,skip=iskip, first=long64(istart*sizepertime),last=long64(iend*sizepertime))
if (reread or n_elements(fluxy0) lt 2) then $
 fluxy0=read_domains('fluxy0.bin',[nx2/nsubdomains,ny2,nz2], ndomains=nsubdomains, $
                   order=[2,1,0],/binary,skip=iskip, first=long64(istart*sizepertime),last=long64(iend*sizepertime))
if (reread or n_elements(fluxz0) lt 2) then $
 fluxz0=read_domains('fluxz0.bin',[nx2/nsubdomains,ny2,nz2], ndomains=nsubdomains, $
                   order=[2,1,0],/binary,skip=iskip, first=long64(istart*sizepertime),last=long64(iend*sizepertime))


if ((reread or n_elements(fluxx1) lt 2) and ndist gt 1) then $
 fluxx1=read_domains('fluxx1.bin',[nx2/nsubdomains,ny2,nz2], ndomains=nsubdomains, $
                   order=[2,1,0],/binary,skip=iskip, first=long64(istart*sizepertime),last=long64(iend*sizepertime))
if ((reread or n_elements(fluxy1) lt 2) and ndist gt 1) then $
 fluxy1=read_domains('fluxy1.bin',[nx2/nsubdomains,ny2,nz2], ndomains=nsubdomains, $
                   order=[2,1,0],/binary,skip=iskip, first=long64(istart*sizepertime),last=long64(iend*sizepertime))
if ((reread or n_elements(fluxz1) lt 2) and ndist gt 1) then $
 fluxz1=read_domains('fluxz1.bin',[nx2/nsubdomains,ny2,nz2], ndomains=nsubdomains, $
                   order=[2,1,0],/binary,skip=iskip, first=long64(istart*sizepertime),last=long64(iend*sizepertime))

if (ndist gt 2 and (reread or n_elements(fluxz2) lt 2)) then $
 fluxx2=read_domains('fluxx2.bin',[nx2/nsubdomains,ny2,nz2], ndomains=nsubdomains, $
                   order=[2,1,0],/binary,skip=iskip, first=long64(istart*sizepertime),last=long64(iend*sizepertime))
if (ndist gt 2 and (reread or n_elements(fluxy2) lt 2)) then $
 fluxy2=read_domains('fluxy2.bin',[nx2/nsubdomains,ny2,nz2], ndomains=nsubdomains, $
                   order=[2,1,0],/binary,skip=iskip, first=long64(istart*sizepertime),last=long64(iend*sizepertime))
if (ndist gt 2 and (reread or n_elements(fluxz2) lt 2)) then $
 fluxz2=read_domains('fluxz2.bin',[nx2/nsubdomains,ny2,nz2], ndomains=nsubdomains, $
                   order=[2,1,0],/binary,skip=iskip, first=long64(istart*sizepertime),last=long64(iend*sizepertime))


