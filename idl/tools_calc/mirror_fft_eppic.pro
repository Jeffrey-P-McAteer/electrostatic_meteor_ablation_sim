;restores fourier transformed data that's been cut off, in 2 of 3D
;uses dimensions from the eppic.i file
;unlike mirror_fft, works on non-square spaces

FUNCTION mirror_fft_eppic, matrix

  @eppic.i
  ndim = ndim_space
  out_size = [ndim,nx*nsubdomains,ny,nz]
  in_size = size(matrix)
  if (in_size[1] ne out_size[1]) && (ndim eq 2) then begin
     tmp = matrix
     matrix = cindgen(out_size[1],in_size[2])*0
     matrix[out_size[1]-in_size[1]:out_size[1]-1,0:in_size[2]-1] = tmp
     in_size=size(matrix)
  endif

  IF ndim eq 2 THEN BEGIN
     if (in_size[1] ne out_size[1]) then begin
        tmp = matrix
        matrix = cindgen(out_size[1],in_size[2])*0
        matrix[out_size[1]-in_size[1]:out_size[1]-1,0:in_size[2]-1] = tmp
        in_size=size(matrix)
     endif

    array = cindgen(out_size[1],out_size[2])*0
    array[*,0:in_size[2]-1] = matrix
    
    mirror = reverse(reverse(matrix,2))
    mirror = conj(mirror)
    array[1:out_size[1]-1,out_size[2]-in_size[2]+1:out_size[2]-1] = mirror[0:in_size[1]-2,0:in_size[2]-2]  ;shifted up & over by 1
    array[0,out_size[2]-in_size[2]+1:out_size[2]-1] = mirror[in_size[1]-1,0:in_size[2]-2]              ;wrap last column to first, up by 1
    
  ENDIF ELSE BEGIN
  
     if (in_size[1] ne out_size[1]) then begin
        tmp = matrix
        matrix = cindgen(out_size[1],in_size[2], in_size[3])*0
        matrix[out_size[1]-in_size[1]:out_size[1]-1,0:in_size[2]-1,0:in_size[3]-1] = tmp
        in_size=size(matrix)
        tmp = 0
     endif

    array = cindgen(out_size[1],out_size[2],out_size[3])*0
    
    IF in_size[2] eq out_size[2] THEN array[*,*,0:in_size[3]-1] = matrix ELSE BEGIN
    
      array[*,0:in_size[2]-1,0:in_size[3]-1] = matrix
      mirror1 = reverse(reverse(reverse(matrix,3),2))
      array[1:out_size[1]-1,out_size[2]-in_size[2]+1:out_size[1]-1,1:in_size[3]-1] = mirror1[0:in_size[1]-2,0:in_size[2]-2,0:in_size[3]-2]
      array[0,out_size[2]-in_size[2]+1:out_size[1]-1,1:in_size[3]-1] = mirror1[in_size[1]-2,0:in_size[2]-2,0:in_size[3]-2]
      matrix = 0
      mirror1 = 0
    ENDELSE
    
    mirror = reverse(reverse(reverse(array[*,*,0:out_size[3]-1],3),2))
    mirror = conj(mirror)
    array[1:out_size[1]-1,1:out_size[2]-1,out_size[3]-in_size[3]+1:out_size[3]-1] = mirror[0:in_size[1]-2,0:out_size[2]-2,0:in_size[3]-2]  ;shifted up & over by 1
    array[1:out_size[1]-1,0,out_size[3]-in_size[3]+1:out_size[3]-1] = mirror[0:in_size[1]-2,in_size[2]-1,0:in_size[3]-2]              ;wrap top to bottom, over by 1
    array[0,1:out_size[2]-1,out_size[3]-in_size[3]+1:out_size[3]-1] = mirror[in_size[1]-1,0:out_size[2]-2,0:in_size[3]-2]              ;wrap back segment to front, up by 1
    mirror = 0
  ENDELSE
  
  return,array
END
