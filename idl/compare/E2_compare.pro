; A routine to produce a comparison of electric field averages
; The second column of the conserved.out file contains the total
; electric field energy in the system

cd, curr=curr
dirroot=curr
dirroot='/project/eregion/meerso/eppic-archive/FBI/'
ndir=5
dir=strarr(ndir)
dir[0]='FBI3D2048E105nuIvme44V120728'
dir[1]='FBI3D2048E105nu.5nuIvme44V120722'
dir[2]='FBI2D512E105nuIvme44TeV120801'
dir[3]='FBI2D512E105nu.5nuIvme44TeV120731'
dir[4]='FBI2D2048E0105nuSqrtme44'

;dir[0]='FBI2D512x512E105nu.5'
;dir[1]='FBI2D512x512E105nu.5np4x'
;dir[2]='FBI2D512x512psi.2dx.2E140nofilter'
;dir[3]='FBI2D512x512psi.2dx.1n04E140nofilter'
;dir[4]='FBI2D2048x2048E70nu.5'
;dir[5]='FBI3D2048E70'
;dir[6]='FBI3D2048E140nu.5V111210'
;dir[7]='FBI3D2048E105nu.5V101219'
;dir[8]='FBI3D2048E105nu.25V101222'
;dir[9]='FBI2D2048x2048nu.25'
;dir[10]='FBI3D2048E105nu.5V101219'
;dir[11]='FBI2D2048x2048E140nu.5'
;dir[12]='FBI3D2048E070nu.5V141210'
;dir[13]='FBI3D2048x512x512E140nu0.5'
;dir[14]='FBI3D2048E105V101210'

cd,dirroot
max_nt=5000
E2all=fltarr(ndir,max_nt)
tvall=fltarr(ndir,max_nt)
nE2all=intarr(ndir)


nd=0
print,'processing: '+dirroot+dir[nd]+'/domain000' 
cd,dirroot+dir[nd]+'/domain000' 
spawn,'cp -fp ../eppic.i '+ dirroot+dir[nd]+'_eppic.i'
@../eppic.i 
if (ndim_space lt 2) then begin ny=1 & dy=1 
if (ndim_space lt 3) then begin nz=1 & dz=1
print,nx*nsubdomains,dx,ny,dy,nz,dz 
con=readarray('conserved.out',2,lineskip=1) 
ncon=(size(con))[2] 
nE2all[nd]=ncon 
E2all[nd,0:ncon-1]=con[1,*]/(eps/2.)/(nx*nsubdomains*dx*ny*dy*nz*dz) 
tvall[nd,0:ncon-1]=findgen(ncon)*dt*nout 

nd=1
print,'processing: '+dirroot+dir[nd]+'/domain000' 
cd,dirroot+dir[nd]+'/domain000' 
spawn,'cp -fp ../eppic.i '+ dirroot+dir[nd]+'_eppic.i'
@../eppic.i  
if (ndim_space lt 2) then begin ny=1 & dy=1
if (ndim_space lt 3) then begin nz=1 & dz=1
print,nx*nsubdomains,dx,ny,dy,nz,dz 
con=readarray('conserved.out',2,lineskip=1) 
ncon=(size(con))[2] 
nE2all[nd]=ncon 
E2all[nd,0:ncon-1]=con[1,*]/(eps/2.)/(nx*nsubdomains*dx*ny*dy*nz*dz) 
tvall[nd,0:ncon-1]=findgen(ncon)*dt*nout 

nd=2
print,'processing: '+dirroot+dir[nd]+'/domain000' 
cd,dirroot+dir[nd]+'/domain000' 
spawn,'cp -fp ../eppic.i '+ dirroot+dir[nd]+'_eppic.i'
@../eppic.i  
if (ndim_space lt 2) then begin ny=1 & dy=1
if (ndim_space lt 3) then begin nz=1 & dz=1
print,nx*nsubdomains,dx,ny,dy,nz,dz 
con=readarray('conserved.out',2,lineskip=1) 
ncon=(size(con))[2] 
nE2all[nd]=ncon 
E2all[nd,0:ncon-1]=con[1,*]/(eps/2.)/(nx*nsubdomains*dx*ny*dy*nz*dz) 
tvall[nd,0:ncon-1]=findgen(ncon)*dt*nout 

nd=3
print,'processing: '+dirroot+dir[nd]+'/domain000' 
cd,dirroot+dir[nd]+'/domain000' 
spawn,'cp -fp ../eppic.i '+ dirroot+dir[nd]+'_eppic.i'
@../eppic.i  
if (ndim_space lt 2) then begin ny=1 & dy=1
if (ndim_space lt 3) then begin nz=1 & dz=1
print,nx*nsubdomains,dx,ny,dy,nz,dz 
con=readarray('conserved.out',2,lineskip=1) 
ncon=(size(con))[2] 
nE2all[nd]=ncon 
E2all[nd,0:ncon-1]=con[1,*]/(eps/2.)/(nx*nsubdomains*dx*ny*dy*nz*dz) 
tvall[nd,0:ncon-1]=findgen(ncon)*dt*nout 

nd=4
print,'processing: '+dirroot+dir[nd]+'/domain000' 
cd,dirroot+dir[nd]+'/domain000' 
spawn,'cp -fp ../eppic.i '+ dirroot+dir[nd]+'_eppic.i'
@../eppic.i  
if (ndim_space lt 2) then begin ny=1 & dy=1
if (ndim_space lt 3) then begin nz=1 & dz=1 
print,nx*nsubdomains,dx,ny,dy,nz,dz 
con=readarray('conserved.out',2,lineskip=1) 
ncon=(size(con))[2] 
nE2all[nd]=ncon 
E2all[nd,0:ncon-1]=con[1,*]/(eps/2.)/(nx*nsubdomains*dx*ny*dy*nz*dz) 
tvall[nd,0:ncon-1]=findgen(ncon)*dt*nout 

; Plotting
nd=0
plot, tvall(nd,0:nE2all[nd]-1), sqrt(E2all[nd,0:nE2all[nd]-1]) $
  ,yrange=sqrt(max(E2all))*[0.01,1.], xrange=[0,max(tvall)] $
, title='RMS E field Comparing many runs', xtitle='Time (s)', ytitle='|E| (mV/m)'

xyouts,tvall(nd,nE2all[nd]-1),sqrt(E2all[nd,nE2all[nd]-1]),dir[nd]

for nd=1,ndir-1 do  $
  oplot, tvall(nd,0:nE2all[nd]-1), sqrt(E2all[nd,0:nE2all[nd]-1]),linestyle= (nd mod 6)
for nd=1,ndir-1 do  $
  xyouts,tvall(nd,nE2all[nd]-1),sqrt(E2all[nd,nE2all[nd]-1]),dir[nd] 

; Put a reference line at the driving field value
oplot,[0,max(Tvall)*1.2],[Ey0_external,Ey0_external],color=128
