;Routine to interpolate around two dimensions of an array in a circular
;arc to collapse these 2 dimensions into 1 dimension.  The center of
;the arc is assumed to be the center of the array.
;
; INPUT
;      array    a 2-D or 3-D array of data.
;               the first 2 dimensions will be interpolated
; Output


function arc_interpolate, array


;Set number per cell
n_per_cell=4

dims=size(array,/dimensions)
npts=min([dims[0],dims[1]])/2
out_dims=[npts]
gt2d=0.0
if (n_elements(dims) gt 2) then $
  gt2d =1.0
if (gt2d) then $
  out_dims=[npts,dims(2:n_elements(dims)-1)]


out=make_array(out_dims,type=size(array,/type))


out[0,*]=array[dims[0]/2,dims[1]/2,*]

if (gt2d) then $
  tarray=transpose(array)
last_dim=size(tarray,/n_dimensions)
for i=1,npts-1 do begin
    theta=findgen(i*n_per_cell*4)/(i*n_per_cell*4)*2*!PI
    ix=i*cos(theta)+dims[0]/2-1
    iy=i*sin(theta)+dims[1]/2-1

    if (gt2d) then $
      out[i,*]=transpose(total(interpolate(tarray,ix,iy),last_dim-1)/(i*n_per_cell*4)) $ 
    else $
      out[i] = (total(interpolate(array,ix,iy),1) /(i*n_per_cell*4))

endfor

return,out
end
