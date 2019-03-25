
PRO toGnuplot,array,filePrefix

; Open file
  openw,99,filePrefix+".gnp"

; Copy array to local array, of correct resize
  local_array = reform(array)

  local_size = size(local_array)
  local_ndim = local_size[0]
  print,"size= ",local_size," ndim = ",local_ndim


; Write header
  printf,99,format='(%"# Data generated in IDL for use in Gnuplot")'
  

; Deal with array size
  CASE local_ndim OF


; if 1D
     1:BEGIN
        print,"You have a 1 D array"
        FOR i=0,local_size[1]-1 DO BEGIN
           printf,99,format='(%"%10d %20.10g")',i,local_array[i]
        END
     END
; if 2D
     2:BEGIN
        print,"You have a 2 D array"
        FOR i=0,local_size[1]-1 DO BEGIN
           FOR j=0,local_size[2]-1 DO BEGIN
              printf,99,format='(%"%10d %10d %20.10g")',i,j,local_array[i,j]
           END
           printf,99,format='(%"")'           
        END
     END
; if 3D
     3:BEGIN
        print,"You have a 3 D array"
     END
; if 4D
     4:BEGIN
        print,"You have a 4 D array"
     END
  ENDCASE
  printf,99,format='(%"\n")'
  close,99
  return
END
