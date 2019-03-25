function inverse_fft_array, data_struct
  
  ;create array
  ft_array = fill_k_array(data_struct)
  
  ;mirror data
  ft_full = mirror_fft_eppic(ft_array)
  
  ;transfrom the data back to real space
  array = fft(ft_full,/INVERSE)

  return, array

end
