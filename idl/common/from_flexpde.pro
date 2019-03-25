
PRO from_FlexPDE,filename,arrout,ndim=ndim,nx=l_nx,ny=l_ny,nz=l_nz

;+
; NAME:
;	from_FlexPDE.pro
;
; PURPOSE:
;	Reads in data written produced by FlexPDE.
;
; CATEGORY:
;	I/O
;
; CALLING SEQUENCE:
;	from_FlexPDE, filename, array, ndim=2
;
; INPUTS:
;	filename:	String name of the file produced by FlexPDE
;
;       arrout:          A varible that will hold the data in the file
;
;
; KEYWORD PARAMETERS:
;	ndim:	Assumed to be 2d, because that is the max that flexpde
;	        can handle (with our current license).
;
; SIDE EFFECTS:
;	'arrout' will be redefined to have the correct dimensions.
;
; RESTRICTIONS:
;	The following syntax must be used to produce the file with 
;       the assumed format:
;
;       contour(n) export (nx+1,ny+1,nz+1) file "filename" export format "#x#b#y#b#1" noheader
;
;       n can be any expression accepted by FlexPDE and nx, ny, nz
;       should have the same meaning as in EPPIC.
;
; MODIFICATION HISTORY:
; 	Written by:	Yann Tambouret, July 2009
;-

  if not keyword_set(ndim) then ndim=2

  nlines=file_lines(filename)

  if (ndim eq 2) then begin
     if not keyword_set(l_nx) then begin
        l_nx=sqrt(nlines)
        l_ny = l_nx
     endif
     if not keyword_set(l_ny) then begin
        l_ny = nlines/l_nx
     endif
     l_nz=1
  endif else begin
     if not keyword_set(l_nx) then begin
        l_nx=nlines^(1/3)
        l_ny = l_nx
        l_nz = l_nx
     endif
     if not keyword_set(l_ny) then begin
        l_ny = sqrt(nlines/l_nx)
        l_nz = l_ny
     endif
     if not keyword_set(l_nz) then begin
        l_ny = nlines/l_nx/l_ny
     endif
  endelse

  ;print,"File name is ",filename," and it has ",nlines," lines."
  openr,99,filename
  var1=0.0
  var2=0.0
  var3=0.0
  var4=0.0
; Read first line

  arrout = findgen(l_nx,l_ny,l_nz)
  for iz=0,l_nz-1 do begin
     for iy=0,l_ny-1 do begin
        for ix=0,l_nx-1 do begin
           if (ndim eq 2) then readf,99,xval,yval,value
           if (ndim eq 3) then readf,99,xval,yval,zval,value
           arrout[ix,iy,iz] = value
        endfor
     endfor
  endfor

  close,99
END
