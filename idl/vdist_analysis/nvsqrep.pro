; Read in nvsqrsities
@params_in.pro

if (n_elements(nvsqrx0) lt 2) then $
 nvsqrx0=read_domains('nvsqrx0.bin',[nx2/nsubdomains,ny2,nz2], ndomains=nsubdomains, $
                   order=[2,1,0],/binary,skip=iskip)
if (n_elements(nvsqry0) lt 2) then $
 nvsqry0=read_domains('nvsqry0.bin',[nx2/nsubdomains,ny2,nz2], ndomains=nsubdomains, $
                   order=[2,1,0],/binary,skip=iskip)
if (n_elements(nvsqrz0) lt 2) then $
 nvsqrz0=read_domains('nvsqrz0.bin',[nx2/nsubdomains,ny2,nz2], ndomains=nsubdomains, $
                   order=[2,1,0],/binary,skip=iskip)


if (ndist ge 1 and n_elements(nvsqrx1) lt 2 and ndist ge 2) then $
 nvsqrx1=read_domains('nvsqrx1.bin',[nx2/nsubdomains,ny2,nz2], ndomains=nsubdomains, $
                   order=[2,1,0],/binary,skip=iskip)
if (ndist ge 1 and n_elements(nvsqry1) lt 2 and ndist ge 2) then $
 nvsqry1=read_domains('nvsqry1.bin',[nx2/nsubdomains,ny2,nz2], ndomains=nsubdomains, $
                   order=[2,1,0],/binary,skip=iskip)
if (ndist ge 1 and n_elements(nvsqrz1) lt 2 and ndist ge 2) then $
 nvsqrz1=read_domains('nvsqrz1.bin',[nx2/nsubdomains,ny2,nz2], ndomains=nsubdomains, $
                   order=[2,1,0],/binary,skip=iskip)

if (ndist ge 2 and n_elements(nvsqrx2) lt 2 and ndist ge 3) then $
 nvsqrx2=read_domains('nvsqrx2.bin',[nx2/nsubdomains,ny2,nz2], ndomains=nsubdomains, $
                   order=[2,1,0],/binary,skip=iskip)
if (ndist ge 2 and n_elements(nvsqry2) lt 2 and ndist ge 3) then $
 nvsqry2=read_domains('nvsqry2.bin',[nx2/nsubdomains,ny2,nz2], ndomains=nsubdomains, $
                   order=[2,1,0],/binary,skip=iskip)
if (ndist ge 2 and n_elements(nvsqrz2) lt 2 and ndist ge 3) then $
 nvsqrz2=read_domains('nvsqrz2.bin',[nx2/nsubdomains,ny2,nz2], ndomains=nsubdomains, $
                   order=[2,1,0],/binary,skip=iskip)


