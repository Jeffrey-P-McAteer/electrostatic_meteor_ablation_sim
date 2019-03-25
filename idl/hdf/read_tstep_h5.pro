function read_tstep_h5, dataname_str, timestamp_str, dataType, ndim=ndim, ndomains=ndomains, $
                        data_size=data_size, ftout=ftout, x2d=x2d, y2d=y2d,z2d=z2d, $
                        BASEPATH=basepath

; NAME:
; READ_TSTEP_H5
;
; PURPOSE:
;       To read in HDF5 data from a series of datasets and
;       combine it into a single array  
;
; CATEGORY:
; Eppic I/O
;
; CALLING SEQUENCE:
;     result = READ_TSTEP_H5('dataname_str','timestamp_str', dataType[, group_idx=group_idx, 
;                   ndomains=ndomains, data_size=data_size, basepath='basepath'])
;
; INPUTS:
;  filename_str:   The name of the desired dataset
; 
; timestamp_str:   String of the desired timestamp
;
;      dataType:   1=h5 non-FT data
;                  2=ph5 non-FT data
;                  3=h5 FT data
;                  4=ph5 FT data
;
;      basepath:   The location of the HDF file, if not the current directory
;
;   REQUIRED FOR NON-PARALLEL READ-IN:
;
;          group_idx:   Index of the time group within the H5 file
; 
;           ndomains:   Number of subdomains.
;         
;          data_size:   Array of the form [*,nx*nsubdomains,ny,nz,...]
;
; OUTPUTS: returns darray, an array of the form [nx*nsubdomains,ny,nz] containing 
;          the data at a single timestep
  
  if not keyword_set(basepath) then CD, Current=basepath
  aug_z = 5 ;augment time string with preceding 0s
  for a = strlen(timestamp_str),aug_z do timestamp_str = '0' + timestamp_str
  print,'Reading timestep ... ' + timestamp_str
  
  if (n_elements(x2d) eq 0) then x2d=0
  if (n_elements(y2d) eq 0) then y2d=0
  if (n_elements(z2d) eq 0) then z2d=0

  if (dataType eq 3) or (dataType eq 4) then begin
     case ndim of
        1: FT={ikx:0, fk:Complex(0.0,0.0)}
        2: FT={ikx:0, iky:0, fk:Complex(0.0,0.0)}
        3: FT={ikx:0, iky:0, ikz:0, fk:Complex(0.0,0.0)}
     endcase
  endif else if (dataType eq 1) then begin
     if ndim eq 3 then darray = findgen(data_size[3],data_size[2],data_size[1]) $
     else darray = findgen(data_size[2],data_size[1])
  endif

  ;
  ;
  ; NON-PARALLEL HDF
  ;
  ;
  if (dataType eq 1 or dataType eq 3) then begin
     for i = 0L, ndomains-1 do begin   ;cycle through domains
        ;determine current domain # and open h5 file
        domain_n = strtrim(string(i),1)
        for k = strlen(domain_n),2 do domain_n = strcompress('0' + domain_n)
        file_path = basepath + '/domain' + domain_n + '.h5'
        file = H5F_OPEN(file_path)
        
        ;determine current timestep and read data & index datasets
        timename = 'time_' + timestamp_str
        dataset_path = timename + '/' + dataname_str
        dataset = H5D_OPEN(file,dataset_path)
        data = H5D_READ(dataset)
        if dataType eq 3 then begin     
        ;
        ; FT READ
        ;
           index_path = dataset_path + '_index'
           index_dataset = H5D_OPEN(file,index_path)
           index = H5D_READ(index_dataset)
        
           ;check for empty datasets
           ft_size = size(data)
           if (ft_size[0] eq 1 and data[0] eq 0) then continue
        
           ;convert data array into complex form
           if (ft_size[0] eq 1) then n_comp = 1 else n_comp = long(ft_size[2])
           data_comp = complex(data[0,0:n_comp-1],data[1,0:n_comp-1])
        
           fstruct_tmp=replicate(FT,n_comp)
        
           ;place data into structure
           for nstr = 0L, n_comp-1 do begin
              fstruct_tmp[nstr].ikx = index[0,nstr]
              if ndim gt 1 then fstruct_tmp[nstr].iky = index[1,nstr]
              if ndim eq 3 then fstruct_tmp[nstr].ikz = index[2,nstr]
              fstruct_tmp[nstr].fk = data_comp[nstr]
           endfor
        
           ;append the new array to the end of existing
           if n_elements(fstruct) eq 0 then fstruct=fstruct_tmp else fstruct=[fstruct,fstruct_tmp]
          
        endif else begin
        ;
        ; NON-FT READ
        ;
           case ndim OF
              2:  begin
                 ; data = reform(data,data_size[2],data_size[1]/ndomains)
                 s_location = (data_size[1]/ndomains)*i
                 e_location = (data_size[1]/ndomains)*(i+1)-1
                 darray[*,s_location:e_location,0] = data
                 data=0
              end
              3:  begin
                 ; data = reform(data,data_size[3],data_size[2],data_size[1]/ndomains)
                 s_location = (data_size[1]/ndomains)*i
                 e_location = (data_size[1]/ndomains)*(i+1)-1
                 darray[*,*,s_location:e_location] = data
                 data=0
              end
           endcase
        endelse
      
        ;release h5 datasets & file
        H5D_CLOSE, dataset
        if dataType eq 0 then H5D_CLOSE, index_dataset
        H5F_CLOSE, file
     endfor  ;end domain run-through  

     if (dataType eq 3) then begin
        if (n_elements(ftout) eq 0) then darray = inverse_fft_array(fstruct) else begin
           ft_array = fill_k_array(fstruct)
           darray = mirror_fft_eppic(ft_array)
        endelse
     endif
  ;
  ;
  ;PARALLEL HDF
  ;
  ;   
  endif else begin
     filename = 'parallel/parallel' + timestamp_str + '.h5'     
     file = H5F_OPEN(basepath + '/' + filename)

     dset = H5D_OPEN(file,'/' + dataname_str)
     data = H5D_READ(dset)
     ;
     ; NON-FT READ
     ;
     if dataType eq 2 then begin

        data_size = size(data)
        ndim = data_size[0]
        
        if ndim eq 2 then begin
           darray=reform(data,data_size[1],data_size[2],1)
        endif else darray = data
        data=0

     endif else begin
     ;
     ; FT-READ
     ;
        index_path = '/' + dataname_str + '_index'
        index_dset = H5D_OPEN(file,index_path)
        index = H5D_READ(index_dset)

        ft_size = size(data)

        ;convert data array into complex form
        if (ft_size[0] eq 1) then n_comp = 1 else n_comp = long(ft_size[2])
        data_comp = complex(data[0,0:n_comp-1],data[1,0:n_comp-1])

        fstruct_tmp=replicate(FT,n_comp)

        ;place data into structure
        for nstr = 0L, n_comp-1 do begin
           fstruct_tmp[nstr].ikx = index[0,nstr]
           if ndim gt 1 then fstruct_tmp[nstr].iky = index[1,nstr]
           if ndim eq 3 then fstruct_tmp[nstr].ikz = index[2,nstr]
           fstruct_tmp[nstr].fk = data_comp[nstr]
        endfor
        data_comp=0
        index=0
        
        ;append the new array to the end of existing
        if n_elements(fstruct) eq 0 then fstruct=fstruct_tmp else fstruct=[fstruct,fstruct_tmp]
        
        if (n_elements(ftout) eq 0) then darray = inverse_fft_array(fstruct) else begin
           ft_array = fill_k_array(fstruct)
           darray = mirror_fft_eppic(ft_array)
           ft_array=0
        endelse
        fstruct=0

        H5D_CLOSE, index_dset
     endelse
     H5D_CLOSE, dset
     H5F_CLOSE, file
  endelse

  ;  non-FT read in puts array in [[z,]y,x] order. correction:
  if (dataType ne 3) && (dataType ne 4) then darray = transpose(darray)

  dsize = size(darray)
  if (x2d ne 0) then return, darray[dsize[1]/2,*,*,*] $
  else if (y2d ne 0) then return, darray[*,dsize[2]/2,*,*] $
  else if (z2d ne 0) then return, darray[*,*,dsize[3]/2,*] $
  else return, darray

end
