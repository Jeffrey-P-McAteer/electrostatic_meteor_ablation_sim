;produces spectra averaged along the meteor
@eppic.i

if (n_elements(spec_plot)) then spec_plot=0

sz = size(den1)
nt=sz[sz[0]]

dentmp = findgen(sz[1],sz[2],sz[4])
print,'averaging in z'
for t = 0, sz[4]-1 do for i = 0, sz[1]-1 do for j = 0, sz[2]-1 do dentmp[i,j,t] = mean(den1[i,j,*,t])

print,'averaging in y'
denavg = findgen(sz[1],sz[4])
for t = 0, nt-1 do for i = 0, sz[1]-1 do denavg[i,t] = mean(dentmp[i,*,t])

H_time=Hanning(nt/2)
for i=0,nt/2-1 do denavg(*,i)=denavg(*,i)*H_time(i)

print,'starting fft'
denfft = shift(abs(fft(denavg-max(denavg))), sz[1]/2, sz[4]/2)

;Since the time transform should be a negative transform, swap the time dimension:
for i = 0, nt/2-1 do begin
    temp=denfft(*,i)
    denfft(*,i)=denfft(*,nt-i-1)
    denfft(*,nt-i-1)=temp
endfor

;Plot only half the spectra
denkw = denfft[sz[1]/2:sz[1]-1,*]

if nx gt ny then kvec= findgen(nx2/2) / nx2 * 2. * !PI / (dx * NOUT_AVG) $
else kvec=findgen(ny2/2) / ny2 * 2. * !PI / (dy * NOUT_AVG)
wv=shift(findgen(nt)*2*!PI/(nt*dt*nout*iskip),nt/2)
wv(0:nt/2-1)=wv(0:nt/2-1)-2*!PI/(dt*nout*iskip)

denkw[0,*]=0

if (spec_plot gt 0) then begin
   range=image_range((denkw))
   print, 'plotting'
   image_plot, denkw[0:range[0,1],range[1,0]:range[1,1]], kvec[0:range[0,1]], wv[range[1,0]:range[1,1]], $
               /zlog, /legend, nlabels=4, title='den(k,w) Spectra', xtitle = 'k (1/m)', ytitle='omega (1/s)'
   date_plot,title
   
   ;Make a second plot with only 2 orders of magnitude 
   range=image_range(denkw gt max(denkw)/100. )
   if range[0,1] eq 0 then range[0,1]=1
   if range[1,0] eq range[1,1] then range[1,1]=range[1,1]+1
   image_plot, denkw[0:range[0,1],range[1,0]:range[1,1]], kvec[0:range[0,1]], wv[range[1,0]:range[1,1]], $
               /zlog, /legend, nlabels=3, title='den(k,w) Spectra', xtitle = 'k (1/m)', ytitle='omega (1/s)' 
   date_plot,title
endif
