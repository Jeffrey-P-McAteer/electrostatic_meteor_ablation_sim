
Pro vtk_output,array,prefix,spacing

;+
; 
; NAME: 
;
; NAME:
;	VTK_OUTPUT
;
; PURPOSE:
;       To write data to vtk legacy format files for use with Paraview.
;
; CATEGORY:
;	General I/O.
;
; CALLING SEQUENCE:
;	VTK_OUTPUT,data,'FileName',[dx,dy,dz,dt]
;
; INPUTS:
;       data:       A 3D or 4D array. 4th demension assumed to be time.
;	FileName:   The name of the file containing data
;       spacing:    a 4D array characterizing the spacing in each
;                   direction; this is optional. 
; OUTPUTS:
;       No outputs...
; 
; TODO:
;       - make opening file safe, test if fileunit is in use.
;-

array_size=size(array)
nx = array_size[1]
ny = array_size[2]
nz = array_size[3]
nt = array_size[4]

if (n_elements(spacing) eq 0) then spacing=[1,1,1,1]

fileunit=11
for t=0,nt-1 do begin
    if (t lt 10) then num_shift='00'
    if ((t lt 100) && (t gt 9)) then num_shift='0'
    if (t gt 99) then num_shift=''
    num_shift = num_shift + string(t,FORMAT="(I3.3)")
    

    openw, fileunit, prefix+num_shift+".vtk"
    

    printf, fileunit, "# vtk DataFile Version 1.4.2"
    printf,fileunit,"Data created by eppic" 
    printf,fileunit,"ASCII" 
    printf,fileunit,"DATASET STRUCTURED_POINTS" 
    printf,fileunit,"DIMENSIONS ",nx,ny,nz
    printf,fileunit,"ORIGIN 0 0 0" 
    printf,fileunit,"SPACING ",$
      spacing[0]," ",$
      spacing[1]," ",$
      spacing[2]
    printf,fileunit,"POINT_DATA ",nx*ny*nz 
    printf,fileunit,"SCALARS data double 1" 
    printf,fileunit,"LOOKUP_TABLE default" 
    
    for z=0,nz-1 do begin
        for y=0,ny-1 do begin
            for x=0,nx-1 do begin
                printf,fileunit,array[x,y,z,t]
            endfor
        endfor
    endfor
    
    close,fileunit
endfor

end
