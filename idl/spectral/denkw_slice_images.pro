; Slices of a 3D array along defined as den(ikx,iky,wv) where we want
; to choose a particular angle along ikx,iky (usually the oblique angles).
; WARNING: only works for nx=ny systems!

pro denkw_slice_images, den, angles, Length, Duration, krange=krange, wrange=wrange, filename=filename, _extra=extra

;Make den(w,k) figures for 2D system
;If filename is specified than save a file with the data for the slice
;with a name starting with filename.

densize=size(den)
nhalf=densize[4]/2

;Put a Hanning window over the temporal component of the data
H_time=Hanning(nhalf)
for i=0,nhalf-1 do den(*,*,*,i)=den(*,*,*,i)*H_time(i)

fden=abs(shift(fft(reform(den[*,*,0,nhalf:nhalf*2-1])),densize[1]/2,densize[2]/2,nhalf/2))

kmin=2*!PI/Length
kxv=(findgen(densize[1])-densize[1]/2)*kmin[0]
kyv=(findgen(densize[2])-densize[2]/2)*kmin[1]
if n_elements(length) gt 2 and densize[3] gt 1 then kzv=(findgen(densize[3])-densize[3]/2)*k[2]

wmin=2*!PI/Duration
wv=(findgen(nhalf)-nhalf/2)*wmin*2 ; One only gets half of the last dimension

; Obtain the appropriate density slice

ikx=findgen(densize[1])
iky=findgen(densize[2])

;if wrange is set, limit the range of w to that given
if n_elements(wrange) eq 2 then begin
   wvr=value_locate(wv,wrange)
   if wvr[0] lt 0 then wvr[0] = 0
   if wvr[1]-wvr[0] lt 3 then print, 'denkw_slice_images: wrange gives to few points.  Ignoring...'
   endif else wvr=[0,nhalf-1]

   
;if krange is set, limit the range of k to that given
;if n_elements(krange) eq 2 then begin

for i=0,n_elements(angles)-1 do begin
   ikxa=(ikx-densize[1]/2)*cos(angles[i]*!PI/180.)+densize[1]/2
   ikya=(iky-densize[2]/2)*sin(angles[i]*!PI/180.)+densize[2]/2

   den_slice=reform(fden(*,0,*))
   
   for iw=0,nhalf-1 do den_slice[*,iw]=interpolate(fden[*,*,iw],ikxa,ikya,missing=0.)

   ks=sqrt( ( (ikxa-densize[1]/2)*kmin[0])^2 + ( (ikya-densize[2]/2)*kmin[1])^2 )
   ks[0:densize[1]/2-1] *= -1.

 ;if krange is set, limit the range of k to that given
   if n_elements(krange) eq 2 then begin
      ksr=value_locate(ks,krange)
      if ksr[0] lt 0 then ksr[0] = 0
      if ksr[1]-ksr[0] lt 3 then print, 'denkw_slice_images: krange gives to few points.  Ignoring...'
   endif else ksr=[0,densize[1]-1]

   if n_elements(filename) gt 0 then begin
      fname=filename+Strcompress(string(angles(i)),/remove_all)+'.sav'
      ds=den_slice[ksr[0]:ksr[1],wvr[0]:wvr[1]]
      kss=ks[ksr[0]:ksr[1]]
      wvs=wv[wvr[0]:wvr[1]]
      save,ds,kss,wvs,ksr,wvr,filename=fname
   endif

   image_plot,den_slice[ksr[0]:ksr[1],wvr[0]:wvr[1]],ks[ksr[0]:ksr[1]],wv[wvr[0]:wvr[1]],/zlog,/legend, $
              title='den(k,w) at angle = '+Strcompress(string(angles(i))), $
              xtitle='k, 1/m (0 -> PI / Debye length )', $
              ytitle = 'w 1/s (0-> PI / (dt*nout_freq)', _extra=extra

endfor

end
