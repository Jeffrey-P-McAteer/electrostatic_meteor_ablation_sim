;restores fourier transformed data that's been cut off, in 2 or 3D
;assumes the space is a square(2D) or cube(3D)
;if the space is not a square/cube, use mirror_fft_eppic, which
;uses the dimensions provided in the eppic.i file

FUNCTION mirror_fft, fft_matrix

  size = size(fft_matrix)

  IF size[0] eq 2 THEN BEGIN
  
    ;y, x
    array = cindgen(size[1],size[1])*0
    array[*,0:size[2]-1] = fft_matrix
    
    mirror = reverse(reverse(fft_matrix,2))
    mirror = conj(mirror)
    array[1:size[1]-1,size[1]-size[2]+1:size[1]-1] = mirror[0:size[1]-2,0:size[2]-2]  ;shifted up & over by 1
    array[0,size[1]-size[2]+1:size[1]-1] = mirror[size[1]-1,0:size[2]-2]              ;wrap last column to first, up by 1
    
  ENDIF ELSE BEGIN
    
    ;z, y, x
    array = cindgen(size[1],size[1],size[1])*0
    IF size[2] eq size[1] THEN array[*,*,0:size[3]-1] = fft_matrix ELSE BEGIN

      array[*,0:size[2]-1,0:size[3]-1] = fft_matrix
      mirror1 = reverse(reverse(reverse(fft_matrix,3),2))
      array[1:size[1]-1,size[1]-size[2]+1:size[1]-1,1:size[3]-1] = mirror1[0:size[1]-2,0:size[2]-2,0:size[3]-2]
      array[0,size[1]-size[2]+1:size[1]-1,1:size[3]-1] = mirror1[size[1]-2,0:size[2]-2,0:size[3]-2]
      
    ENDELSE
      
    mirror = reverse(reverse(reverse(array[*,*,0:size[3]-1],3),2))
    mirror = conj(mirror)
    array[1:size[1]-1,1:size[1]-1,size[1]-size[3]+1:size[1]-1] = mirror[0:size[1]-2,0:size[1]-2,0:size[3]-2]  ;shifted up & over by 1
    array[1:size[1]-1,0,size[1]-size[3]+1:size[1]-1] = mirror[0:size[1]-2,size[1]-1,0:size[3]-2]              ;wrap top to bottom, over by 1
    array[0,1:size[1]-1,size[1]-size[3]+1:size[1]-1] = mirror[size[1]-1,0:size[1]-2,0:size[3]-2]              ;wrap last column to first, up by 1
        
  ENDELSE
  
  return,array
  
END