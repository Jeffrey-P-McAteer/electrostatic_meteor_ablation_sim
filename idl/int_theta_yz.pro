; This routine integrates a 4D array assuming that dimensions 2 and 3
; are y, z and that they should be converted to r, theta and then
; integrated around theta.
;
; Output: A 3D array
;
Function int_theta_y, a, y, z


; Do the first integration just to get the size of the output
a1rz=int_theta(reform(a[0,*,*,0]),y,z,/grid)
nth=nelements(a1rz)

na4=size(a4D)

a3D=make_array(na4[1],
; Make the output array

;Begin the loop

return, a3D
end

