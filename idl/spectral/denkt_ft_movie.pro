;program to read in the ft data and produce a movie of den1(kx,ky,t)
;pro, denkt_ft_movie
@eppic.i

; A factor to arbitrarily reduce k-space further (generally for memory savings)
kreduce=1

;Obtain a list of all files in denft1
cd, 'domain000/denft1'
filename_list=file_search('d*.bin')
cd, '../..'

; Gather and assemble data from all times
nts=n_elements(filename_list)

nx3=nx*nsubdomains

; Define the range desired (Throw out the top 3/4 of the spectum by default
if n_elements(denft_out_kmax1) eq 0 then denft_out_kmax1=2*!PI/dx *(1./8.)
; Fraction of total array desired
kfrac=round((2*!PI/dx )/denft_out_kmax1)*kreduce
ikxmax=nx3/2+round(nx3/2/kfrac)-1
ikxmin=nx3/2-round(nx3/2/kfrac)
ikymax=ny/2+round(ny/2/kfrac)-1
ikymin=ny/2-round(ny/2/kfrac)
ikzmax=nz/2+round(nz/2/kfrac)-1
ikzmin=nz/2-round(nz/2/kfrac)

if ndim_space eq 3 then $
  ikdim=[ikxmax-ikxmin+1,ikymax-ikymin+1,ikzmax-ikzmin+1] $
else ikdim=[ikxmax-ikxmin+1,ikymax-ikymin+1] $


for it=0,n_elements(filename_list)-1 do begin

  print, 'Reading: ' +  filename_list[it]
  d1ft=read_ft('domain*/denft1/'+filename_list[it],ndim_space)

  dktmp=fill_k_array(d1ft)
; Shift to center the image
  dktmp=shift(dktmp,nx3/2,ny/2,0)

  if (it eq 0) then begin
      dks=size(dktmp)

      dim=[ikdim[0:dks[0]-1],[nts]]
      dk=make_array(dim,type=dks[dks[0]+1])

      ; The ranges of ikx,iky,ikz are:
      ikx_r=[nx3/2-ikdim[0]/2,nx3/2+ikdim[0]/2-1]
      iky_r=[ny/2-ikdim[1]/2,ny/2+ikdim[1]/2-1]
      ikz_r=[0,ikdim[2]-1]

  endif

; Place the correct component of dktmp into dk

  dk[*,*,0:min([ikz_r[1],dks[3]-1]),it] = $
    dktmp[ikx_r[0]:ikx_r[1],iky_r[0]:iky_r[1],ikz_r[0]:min([ikz_r[1],dks[3]-1])]

endfor

dktmp=0

;Set up the k vector for the axes
dxdk=dx*nx3/dim(0)
kxv=shift(findgen(dim(0))*2*!PI/(nx3*dx),dim(0)/2)
kxv(0:dim(0)/2-1)=kxv(0:dim(0)/2-1)-2*!PI/(dxdk)
if phisize(0) ge 2 then begin
    kyv=shift(findgen(phisize(2))*2*!PI/(ny*dy),ny2/2)
    kyv(0:ny2/2-1)=kyv(0:ny2/2-1)-2*!PI/(dy*nout_avg)
endif
if phisize(0) ge 3 then begin
    kzv=shift(findgen(phisize(3))*2*!PI/(nz*dz),nz2/2)
    kzv(0:nz2/2-1)=kzv(0:nz2/2-1)-2*!PI/(dz*nout_avg)
endif


end
