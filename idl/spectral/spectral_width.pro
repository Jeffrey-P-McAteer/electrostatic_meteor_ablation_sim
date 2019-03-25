function spectral_width,data,ft=ft,bin=bin,h5=h5,plot=plot,ps=ps,indir=indir,outdir=outdir
;returns a 2d array containing [k,sigma] for a range of ks

@eppic.i
  cd, current=c
  if n_elements(outdir) eq 0 then outdir = c + '/'
  if n_elements(indir) eq 0 then indir = c + '/'

  if ndim_space lt 3 then nz=1
  if ndim_space lt 3 then dz=0.
  if ndim_space lt 2 then ny=1
  if ndim_space lt 2 then dy=0.

  nx2 = long(nx*nsubdomains/nout_avg) > 1
  ny2 = long(ny/nout_avg) > 1
  nz2 = long(nz/nout_avg) > 1
  nx3 = nx2/nsubdomains
  
  if n_elements(ft) gt 0 then begin ;ft h5
     den1 = read_h5(indir + data,0,/ft)
  endif else if n_elements(h5) gt 0 then begin ;regular h5
     den1 = read_h5(indir + data,1)
  endif else begin ;binary
     if n_elements(iskip) eq 0 then iskip = 1
     if (n_elements(istart) eq 0) then istart=0
     if (n_elements(iend) eq 0) then iend=-1
     sizepertime=long64((nx3*1.0*ny2))
     if(ndim_space eq 3) then sizepertime = sizepertime*nz2
     den1=read_domains(indir + data + '.bin',[nx3,ny2,nz2], ndomains=nsubdomains, $
                       order=[2,1,0],/binary,skip=iskip,$
                       first=long64(istart*sizepertime),last=long64(iend*sizepertime))
  endelse

  spec_plot=0
  @meteor_spectra
  
  krange = n_elements(kvec)
  width = findgen(2,krange)
  width[0,*] = kvec
  for i =  1, krange-1 do begin
     denf = denkw[i,*]
     gfit = gaussfit(wv/kvec[i],denf,coeff)
     width[1,i] = coeff[2]
  endfor

  if (n_elements(plot) ne 0) then begin ;plot sigma(k) 
     if (n_elements(ps) ne 0) then ps,outdir+'spectral_width.ps',/color
     !p.multi = [0,2,3]
     
;     for i =  1, krange-1 do begin
        ;plot,wv/kvec[i],denf,title='Spectra at ' + string(kvec[i]) + 'm',xtitle='Vphase (m/s)',ytitle='Power (dB)'
        ;oplot,wv/kvec[i],gfit,thick=3
;     endfor
     
     plot,width[0,*],width[1,*],title='Spectral Width',xtitle='k (1/m)',ytitle='sigma'
     
     if (n_elements(ps) ne 0) then stop_ps
  endif ;plot
  
  return,width
  
end
