; Make an image of den(w,|k|) by averaging k in a circle.

if (n_elements(den0) le 2) then begin
@den3d
endif

;To reduce memory usage:
den1=0.

dsize=size(reform(den0))
if dsize[0] eq 4 then begin
    denfft=reform(den0[*,*,0,*])
    dsize=size(reform(denfft))
endif else if dsize[0] eq 3 then denfft=reform(den0) $
else message, 'Input array not required 3-D'

;To reduce memory usage:
den0=0.

nx=dsize[1]
ny=dsize[2]
kx_min=2*!PI/(nx*dx)
ky_min=2*!PI/(ny*dy)
k_min=sqrt(kx_min^2+ky_min^2)

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

denkw=arc_interpolate(denfft)

if nx gt ny then kvec= findgen(nx2/2) / nx2 * 2. * !PI / (dx * NOUT_AVG) $
else kvec=findgen(ny2/2) / ny2 * 2. * !PI / (dy * NOUT_AVG)
wv=shift(findgen(nt)*2*!PI/(nt*dt*nout*iskip),nt/2)
wv(0:nt/2-1)=wv(0:nt/2-1)-2*!PI/(dt*nout*iskip)

range=image_range(denkw gt max(denkw)/1000. )

image_plot, denkw[0:range[0,1],range[1,0]:range[1,1]], kvec[0:range[0,1]], wv[range[1,0]:range[1,1]], $
 /zlog, /legend, nlabels=4, $
  title='Total Spectral Energy: den(|k|,w)', xtitle = '|k| (1/m)', ytitle='omega (1/s)'
date_plot,title

; Make a second plot with only 2 orders of magnitude 
range=image_range(denkw gt max(denkw)/100. )
if range[0,1] eq 0 then range[0,1]=1
if range[1,0] eq range[1,1] then range[1,1]=range[1,1]+1
image_plot, denkw[0:range[0,1],range[1,0]:range[1,1]], kvec[0:range[0,1]], wv[range[1,0]:range[1,1]], $
 /zlog, /legend, nlabels=3, $
  title='Total Spectral Energy: den(|k|,w)', xtitle = '|k| (1/m)', ytitle='omega (1/s)'
date_plot,title


end
