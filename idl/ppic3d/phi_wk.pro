; Generate the ky=0 dispersion plot

@phi3d
if (n_elements(iky) eq 0) then iky=0

nx2=nx/nout_avg > 1
ny2=ny/nout_avg > 1
nz2=nz/nout_avg > 1


Esize=size(phi)
nt=Esize(4)
    
if (n_elements(E_in_k_space) eq 0) then E_in_k_space='false'
if (E_in_k_space ne 'true') then begin
        
  ; Temporal damping:    
  H_time=Hanning(nt)
  for i=0,nt-1 do begin
    phi(*,*,*,i)=phi(*,*,*,i)*H_time(i)
  endfor
   
  ;Convert phi into k space
  phi=fft(phi,1,/overwrite)
  phi=shift(phi,nx2/2,ny2/2,nz/2,nt/2)
  E_in_k_space='true'
    
endif

kxv=shift(findgen(Esize(1))*2*!PI/(nx*dx),nx2/2)
kxv(0:nx2/2-1)=kxv(0:nx2/2-1)-2*!PI/(dx*nout_avg)

kyv=shift(findgen(Esize(2))*2*!PI/(ny*dy),ny2/2)
kyv(0:ny2/2-1)=kyv(0:ny2/2-1)-2*!PI/(dy*nout_avg)


if dz eq 0 then dz=dx
kzv=shift(findgen(Esize(3))*2*!PI/(nz*dz),nz2/2)
nzm=max([nz2/2-1,0])
kzv(0:nzm)=kzv(0:nzm)-2*!PI/(dz*nout_avg)

wv=shift(findgen(nt)*2*!PI/(nt*dt*nout*iskip),nt/2)
wv(0:nt/2-1)=wv(0:nt/2-1)-2*!PI/(dt*nout*iskip)

    ksqr=fltarr(nx2,ny2,nz2)
    for ikx=0, nx2-1 do $
    for iky=0, ny2-1 do $
    for ikz=0, nz2-1 do $
      ksqr(ikx,iky,ikz)=kxv(ikx)^2+kyv(iky)^2+kzv(ikz)^2
    
ikx_min=nx2/2
ikx_max=(nx2/2+1) mod nx
iky_min=ny2/2
iky_max=(ny2/2+1) mod ny
ikz_min=nz2/2
ikz_max=(nz2/2+1) mod nz


end
