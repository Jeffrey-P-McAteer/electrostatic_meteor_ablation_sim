; A routine to read in and scale velocity distributions from ppic3d data
@params_in.pro

@translate_vars

if (n_elements(iskip) eq 0) then iskip = 1


if (n_elements(pnvy0) eq 0 or ndim_space eq 1) then pnvy0=1
if (n_elements(pnvz0) eq 0  or ndim_space eq 1) then pnvz0=1
vsize=long(pnvx0)*pnvy0*pnvz0
vdist0=readarray('vdist0.bin',[pnvx0,pnvy0,pnvz0],order=[2,1,0],/binary,skip=iskip)
;Define the coordinates
vdistvx0=findgen(pnvx0)*(pvxmax0-pvxmin0)/pnvx0+pvxmin0
if pnvy0 gt 1 then vdistvy0=findgen(pnvy0)*(pvymax0-pvymin0)/pnvy0+pvymin0
if pnvz0 gt 1 then vdistvz0=findgen(pnvz0)*(pvzmax0-pvzmin0)/pnvz0+pvzmin0
;Define the time coordinate
nvdist=(size(vdist0))(4)
if n_elements(vdist_out_subcycle0) eq 0 then vdist_out_subcycle0=1
vdistt0=findgen(nvdist)*dt*nout*vdist_out_subcycle0*iskip

if (n_elements(pnvy1) eq 0 or ndim_space eq 1) then pnvy1=1
if (n_elements(pnvz1) eq 0  or ndim_space eq 1) then pnvz1=1
vsize=long(pnvx1)*pnvy1*pnvz1
if (ndist ge 2) then $
vdist1=readarray('vdist1.bin',[pnvx1,pnvy1,pnvz1],order=[2,1,0],/binary,skip=iskip)
;Define the coordinates
vdistvx1=findgen(pnvx1)*(pvxmax1-pvxmin1)/pnvx1+pvxmin1
if pnvy1 gt 1 then vdistvy1=findgen(pnvy1)*(pvymax1-pvymin1)/pnvy1+pvymin1
if pnvz1 gt 1 then vdistvz1=findgen(pnvz1)*(pvzmax1-pvzmin1)/pnvz1+pvzmin1
;Define the time coordinate
nvdist=(size(vdist1))(4)
if n_elements(vdist_out_subcycle1) eq 0 then vdist_out_subcycle1=1
vdistt1=findgen(nvdist)/nvdist*dt*nout*vdist_out_subcycle1*iskip


if (n_elements(pnvx2) eq 0) then pnvx2=1
if (n_elements(pnvy2) eq 0 or ndim_space eq 1) then pnvy2=1
if (n_elements(pnvz2) eq 0 or ndim_space eq 1) then pnvz2=1
vsize=long(pnvx2)*pnvy2*pnvz2
if (ndist ge 3) then $
  vdist2=readarray('vdist2.bin',[pnvx2,pnvy2,pnvz2],order=[2,1,0],/binary,skip=iskip)
;Define the coordinates
if pnvx2 gt 1 then vdistvx2=findgen(pnvx2)*(pvxmax2-pvxmin2)/pnvx2+pvxmin2
if pnvy2 gt 1 then vdistvy2=findgen(pnvy2)*(pvymax2-pvymin2)/pnvy2+pvymin2
if pnvz2 gt 1 then vdistvz2=findgen(pnvz2)*(pvzmax2-pvzmin2)/pnvz2+pvzmin2
;Define the time coordinate
if n_elements(vdist2) gt 2 then  nvdist=(size(vdist2))(4) else nvdist=1
if n_elements(vdist_out_subcycle2) eq 0 then vdist_out_subcycle2=1
vdistt2=findgen(nvdist)/nvdist*dt*nout*vdist_out_subcycle2*iskip

