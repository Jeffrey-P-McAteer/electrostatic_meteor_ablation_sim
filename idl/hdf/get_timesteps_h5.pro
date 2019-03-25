function get_timesteps_h5,dataname_str, dataType, get_times=get_times
  
  ;this is an inefficient method
  ;revising the system used to determine subcycling in eppic 
  ;to include all subcycling variables in eppic.i would be worthwhile

  ;if using parallel output & get_timesteps is returning giving an error,
  ;check that no parallel*.h5 files are left in the working directory 
  ;from a different run

  ;returns ntimes_array, containing the number of good timesteps, subcycle, & their locations
  
  ;if run has not been completed, data has unusual subcycling, or files from
  ;old runs exist in the same directory throw get_times flag to manually check which 
  ;timesteps have good data (otherwise assuems all files & timesteps have relevent data)

  CD, Current=pathbase
  t=0
  
  if (get_times ne 0) then get_times=1

  case get_times of
     0: begin
        
        if (dataType eq 1) or (dataType eq 3) then begin
           file = H5F_OPEN(pathbase + '/domain000.h5')
           ntimes = H5G_GET_NUM_OBJS(file)
        endif

        if (dataType eq 2) or (dataType eq 4) then begin
           files = file_search(pathbase + '/parallel/parallel*.h5')
           ntimes = n_elements(files)
        endif
        
        ntimes = (ntimes-1) + 1
        idx = indgen(ntimes)
        ntimes_array = [ntimes,idx]
     end

     1: begin

        if (dataType eq 2) or (dataType eq 4) then begin
           files = file_search(pathbase + '/parallel/parallel*.h5')
           times = n_elements(files)
           timestrs = strarr(times)
           ntimes = intarr(times)
           
           ;get filenames
           for i = 0, times-1 do begin
              f = files[i]
              ;identifies number of timesteps desired data is output in
              file = H5F_OPEN(f)
              dsets = H5G_GET_NUM_OBJS(file)
              for ndata = 0, dsets-1 do begin
                 name = H5G_GET_OBJ_NAME_BY_IDX(file,ndata)
                 if (name eq dataname_str) then begin
                    t++
                    ntimes[i] = 1
                 endif
              endfor
           endfor
           
           H5F_CLOSE, file
        endif else begin
           file = H5F_OPEN(pathbase + '/domain000.h5')
           times = H5G_GET_NUM_OBJS(file)
           ntimes = intarr(times)
           
           ;gets all group names and cuts off the restart timestamp
           timestrs = strarr(times)
           for i = 0, times-1 do begin
              group = H5G_GET_OBJ_NAME_BY_IDX(file,i)
              tbase = strmid(group, 0, 10)
              timestrs[i] = tbase
           endfor
           
           ;identifies repeat data and sets one copy to '0'
           for i = 0, times - 2 do begin
              for j = i + 1, times - 1 do begin
                 if timestrs[i] eq timestrs[j] then timestrs[i] = '0'
              endfor
           endfor
           
           ;identifies number of timesteps desired data is output in
           for i = 0, times - 1 do begin
              if timestrs[i] eq '0' then ntimes[i] = 0 else begin
                 group = H5G_GET_OBJ_NAME_BY_IDX(file,i)
                 group = H5G_OPEN(file,'/' + group)
                 dsets = H5G_GET_NUM_OBJS(group)
                 for ndata = 0, dsets-1 do begin
                    name = H5G_GET_OBJ_NAME_BY_IDX(group,ndata)
                    if (name eq dataname_str) then begin
                       t++
                       ntimes[i] = 1
                    endif
                 endfor
                 H5G_CLOSE, group
              endelse
           endfor
           H5F_CLOSE, file
        endelse
        
        good_data = where(ntimes eq '1')
        ndata = n_elements(good_data)
        ntimes_array=intarr(1+ndata)
        ntimes_array[0]=t
        ntimes_array[1:ndata] = good_data
     end  
  endcase

  tsize = size(ntimes_array)
  return, ntimes_array

end
