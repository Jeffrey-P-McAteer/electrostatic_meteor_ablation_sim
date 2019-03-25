if (n_elements(den1) le 2) then begin
@den3d
endif
nt=densize(4)
if n_elements(den1f) lt 2 then begin
    den1f=den1
    for i=0,nt-1 do den1f(*,*,*,i)=shift(reform(abs(fft(den1(*,*,*,i))),nx2,ny2,nz2),nx2/2,ny2/2,nz2/2)
    for i=0,nt-1 do den1f(nx2/2,ny2/2,nz2/2,i)=0
endif

ix1=nx2/2
iy1=ny2/2

plot,tv,den1f(ix1+1,iy1,0,*),/ylog,xtitle='time (s)',ytitle='den(k,t) amplitude'
xyouts,tv[nt-1],den1f(ix1+1,iy1,0,nt-1),'Mode (1,0)'
oplot,tv,den1f(ix1+4,iy1,0,*),linestyle=2
xyouts,tv[nt-1],den1f(ix1+4,iy1,0,nt-1),'Mode (4,0)'
oplot,tv,den1f(ix1+20*4,iy1-6*4,0,*),linestyle=1
xyouts,tv[nt-1],den1f(ix1+20*4,iy1-6*4,0,nt-1),'Mode (20,-6)'

end
