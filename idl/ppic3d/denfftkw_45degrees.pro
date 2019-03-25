; Make an image of den(w,kx,-kx) by just slicing those k out.

if (n_elements(den0) le 2) then begin
@den3d
; Use only the second half of the simulation:
nden=densize[densize[0]]
den0=den1[*,*,*,nden/2:nden-1]
densize=size(den0)
nden=densize[densize[0]]
endif
;To reduce memory usage:
den1=0.

;Use the second half of the data only
dsize=size(reform(den0))
if dsize[0] eq 4 then begin
    denfft=reform(den0[*,*,0,*])
    dsize=size(reform(denfft))
endif else if dsize[0] eq 3 then denfft=reform(den0) $
else message, 'Input array not required 3-D'

;To reduce memory usage:
den0=0.

;Fourier transform the input matrix:
denfft=shift(abs(fft(denfft)), dsize[1]/2, dsize[2]/2, dsize[3]/2)

;Since the time transform should be a negative transform, swap the
;time dimension:
nt=dsize[dsize[0]]
for i = 0, nt/2-1 do begin
    temp=denfft(*,*,i)
    denfft(*,*,i)=denfft(*,*,nt-i-1)
    denfft(*,*,nt-i-1)=temp
endfor

;Zero the D.C. component
denfft(dsize[1]/2, dsize[2]/2, dsize[3]/2)=0


if dsize(1) ne dsize(2) then message, "Error: Only works on square array"
denfft45=reform(denfft(*,0,*))
for i=0,dsize(1)-1 do denfft45(i,*)=denfft(i,dsize(1)-1-i,*)

kv45=shift(findgen(densize(1))*2*!PI/(nx*dx),nx2/2)
kv45(0:nx2/2-1)=kv45(0:nx2/2-1)-2*!PI/(dx*nout_avg)
kv45 *= sqrt(2.)/2.

wv=shift(findgen(nt)*2*!PI/(nt*dt*nout*iskip),nt/2)
wv(0:nt/2-1)=wv(0:nt/2-1)-2*!PI/(dt*nout*iskip)

range=image_range(denfft45 gt max(denfft45)/100.)

image_plot,denfft45[range[0,0]:range[0,1],range[1,0]:range[1,1]], $
  kv45[range[0,0]:range[0,1]], wv[range[1,0]:range[1,1]], $
  /legend,/zlog, $
  title='den(0:*,*:0,w)',xtitle='k at -45 degrees (m)',ytitle='omega (1/s)'

image_plot,denfft45[range[0,0]:range[0,1],range[1,0]:range[1,1]], $
  kv45[range[0,0]:range[0,1]], wv[range[1,0]:range[1,1]], $
  /legend, $
  title='den(0:*,*:0,w)',xtitle='k at -45 degrees (1/m)',ytitle='omega (1/s)'


; Now make for +45 degrees
denfft45=reform(denfft(*,0,*))
for i=0,dsize(1)-1 do denfft45(i,*)=denfft(i,i,*)
range=image_range(denfft45 gt max(denfft45)/100.)

image_plot,denfft45[range[0,0]:range[0,1],range[1,0]:range[1,1]], $
  kv45[range[0,0]:range[0,1]], wv[range[1,0]:range[1,1]], $
  /legend,/zlog, $
  title='den(0:*,*:0,w)',xtitle='k at +45 degrees (m)',ytitle='omega (1/s)'

image_plot,denfft45[range[0,0]:range[0,1],range[1,0]:range[1,1]], $
  kv45[range[0,0]:range[0,1]], wv[range[1,0]:range[1,1]], $
  /legend, $
  title='den(0:*,*:0,w)',xtitle='k at +45 degrees (1/m)',ytitle='omega (1/s)'

end
