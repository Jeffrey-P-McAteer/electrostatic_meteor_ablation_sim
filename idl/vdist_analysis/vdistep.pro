; Read in vdist arrays
@params_in.pro

if n_elements(reread) eq 0 then reread=1
if (n_elements(iskip) eq 0) then iskip = 1

if (n_elements(pnvy0) eq 0 or ndim_space eq 1) then pnvy0=1
if (n_elements(pnvz0) eq 0  or ndim_space eq 1) then pnvz0=1

if (reread or n_elements(vdist0) lt 2) then $
 vdist_tmp=read_domains('vdist0.bin',[pnvx0,pnvy0,pnvz0], ndomains=nsubdomains, $
                   order=[2,1,0],/binary,skip=iskip)
; Each domain is now a piece of the matrix concatinated in the
; x-direction. Combine these:
vdist0=vdist_tmp[0:pnvx0-1,*,*,*]
for i=1,nsubdomains-1 do vdist0 += vdist_tmp[i*pnvx0:(i+1)*pnvx0-1,*,*,*]

if ((reread or n_elements(vdist1) lt 2) and ndist gt 1) then $
 vdist_tmp=read_domains('vdist1.bin',[pnvx1,pnvy1,pnvz1], ndomains=nsubdomains, $
                   order=[2,1,0],/binary,skip=iskip)
; Each domain is now a piece of the matrix concatinated in the
; x-direction. Combine these:
vdist1=vdist_tmp[0:pnvx1-1,*,*,*]
for i=1,nsubdomains-1 do vdist1 += vdist_tmp[i*pnvx1:(i+1)*pnvx1-1,*,*,*]


if ((reread or n_elements(vdist2) lt 2) and ndist gt 2) then $
  vdist_tmp=read_domains('vdist2.bin',[pnvx2,pnvy2,pnvz2], ndomains=nsubdomains, $
                         order=[2,1,0],/binary,skip=iskip)
; Each domain is now a piece of the matrix concatinated in the
; x-direction. Combine these:
if n_elements(pnvx2) ne 0 then $
   vdist2=vdist_tmp[0:pnvx2-1,*,*,*]
if n_elements(pnvx2) ne 0 then $
   for i=1,nsubdomains-1 do $
      vdist2 += vdist_tmp[i*pnvx2:(i+1)*pnvx2-1,*,*,*]

vdist_tmp=0



