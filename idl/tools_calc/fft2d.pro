Function fft2d, Ain, _extra=e
;
; NAME:
;          fft2d
; PURPOSE:
;          To fft in 2d an n-d array where the first 2 dimensions are
;          transformed and the remaining dimensions are iterated
;          through.
;
Asize=size(Ain)
if (Asize[0] lt 2)then message, 'FFT2D cannot fourier transform a 1-D array'

Aout=complex(Ain)

if Asize[0] gt 2 then i1max=Asize[3]-1 else i1max=0
if Asize[0] gt 3 then i2max=Asize[4]-1 else i2max=0
if Asize[0] gt 4 then i3max=Asize[5]-1 else i3max=0

for i1=0, i1max do begin
    for i2=0, i2max do begin
        for i3=0, i3max do begin
            Aout[*,*,i1,i2,i3]=shift(fft(Ain[*,*,i1,i2,i3], $
                                         _EXTRA = e),Asize[1]/2,Asize[2]/2)
        endfor
    endfor
endfor

return, Aout

end
