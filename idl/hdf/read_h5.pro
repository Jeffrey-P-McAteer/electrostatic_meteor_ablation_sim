function read_h5, dataname_str, dataType=dataType, startt=startt, endt=endt, start_frac=start_frac,$
                  end_frac=end_frac, subcycle=subcycle, skip=skip, ftout=ftout, x2d=x2d, y2d=y2d, $
                  z2d=z2d, t1=t1, get_times=get_times,den_norm=den_norm,full_array=full_array, $
                  dims=dims
;
; NAME:
;       READ_H5
;
; PURPOSE:
;       To read in HDF5 data of EPPIC runs from a series of datasets and
;       combine it into a single array or size [nx*nsubdomains,ny,nz,nt/nout+1)
;       (equivalent of read_domains binary routine)
;
; CATEGORY:
;       Eppic I/O
;
; CALLING SEQUENCE:
;       result = READ_H5('dataname_str'[,dataType=dataType, startt=startt, endt=endt, $
;                         start_frac=start_frac, end_frac=end_frac,subcycle=subcycle, $
;                         skip=skip, ftout=ftout, x2d=x2d, y2d=y2d, z2d=z2d, t1=t1,   $
;                         get_times=get_times])
;
; INPUTS:
;       dataname_str:   The name of the datasets containing the data.
;       HDF datasets are given the same names as their binary equivalents.
;      
; KEYWORDS: 
;      dataType:   1=h5 non-FT data
;                  2=ph5 non-FT data
;                  3=h5 FT data
;                  4=ph5 FT data
;                  This keyword only needs to be set if multiple
;                  forms of the data exist in the same directory; otherwise the routine
;                  should identify the correct type
;
;    startt,endt:  timestep to start/end read in (default to 0, nt)
;
;  start_frac/end_frac:  fraction of data to start/end read in
;                  (default to 0, 1)
;
;       subcycle:  if known, output subcycle of data (defaut assumes
;                  data is output at every write)
;
;           skip:  reads in every skip-th output (defaults to 1) -
;                  in units of output timesteps, not simulation
;                  timesteps (i.e. skip=8 reads in every 8th output,
;                  not every 8th simulation timestep) 
;
;          ftout:  if data is in FT form, do not take inverse fourier
;                  transform and return spectral output (off by default)
;
;        x/y/z2d:  if data is 3D, reads in slice of data along the
;                  center of x/y/z (/x2d reads in [nx/2,*,*,*]) (off by default)
;
;             t1:  reads in only one timestep (off by default)
;
;      get_times:  set if an older data set could be in the same
;                  directory or if data
;                  has unusal/unknown subcycling to manually check which
;                  output timesteps have good data (off by default)
;
;       den_norm:  automatically renormalize density output (off by default)
;
;     full_array:  set if using FT output with infrequent non-FT output
;
;           dims:  optionally specify array dimensions
;
; OUTPUT:
;     data_array: A 4-dimensional array containing the data

  @eppic.i
  @params_in.pro
  ndim = ndim_space
  ndomains = nsubdomains

  ; check keywords & set defaults for unassigned values
  if (n_elements(skip) eq 0) then skip = 1
  if (n_elements(get_times) eq 0) then get_times = 0
  if (n_elements(full_array) eq 0) then full_array = 0
  if (ndim lt 3) then nz = 1
  if (n_elements(start_frac) ne 0) then startt = start_frac*(nt/nout+1)*nout $
  else if (n_elements(startt) eq 0) then startt = 0.
  if (n_elements(end_frac) ne 0) then endt = end_frac*(nt/nout+1)*nout $
  else if (n_elements(endt) eq 0) then endt = nt
  if (n_elements(t1) ne 0) then begin
     startt = long(t1)
     endt = long(t1)
  endif
  if (startt mod nout ne 0) or (endt mod nout ne 0) then begin
     startt=long(startt)
     endt=long(endt)
     while startt mod nout ne 0 do startt++
     while endt mod nout ne 0 do endt++
  endif
  if (endt gt nt) then endt=nt
  if (startt ne 0) then print,'Starting read-in at timestep ',startt
  if (endt ne nt) then print,'Ending read-in at timestep ',endt  
  if (n_elements(den_norm) eq 0) then begin
     den_norm = 0
     if (dataname_str eq 'den0') or (dataname_str eq 'den1') then den_norm = 1
  endif

  ; identify data type (parallel vs. non-parallel & ft vs. non ft)
  if (n_elements(hdf_output_arrays) eq 0) then hdf_output_arrays = 0
  if (n_elements(ft_output_arrays) eq 0) then ft_output_arrays = 0
  if (n_elements(dataType) eq 0) then dataType = hdf_output_arrays + 2*ft_output_arrays
  if (dataType gt 4) or (dataType eq 0) then stop, 'Could not identify valid dataType.  HDF5 data may not exist in this directory; try specifying dataType in function call.'

  ; check subcycle
  if (n_elements(flux_out_subcycle0) eq 0) then flux_out_subcycle0 = 1
  if (n_elements(vdist_out_subcycle0) eq 0) then vdist_out_subcycle0 = 1
  if (n_elements(nvsqr_out_subcycle0) eq 0) then nvsqr_out_subcycle0 = 1
  dataname_base = strmid(dataname_str,0,4)
  if n_elements(subcycle) eq 0 then begin
     case dataname_base of
        'flux': begin
           subcycle=flux_out_subcycle0*nout
        end
        'vdis': begin
           subcycle=vdist_out_subcycle0*nout
        end           
        'nvsq': begin
           subcycle=nvsqr_out_subcycle0*nout
        end
        else: subcycle = nout
     endcase
  endif

  ; reset dimensions if reading in a 2d slice
  if (dataType ge 3) then nout_avg = 1 ;spectral data is not averaged spatially
  if (n_elements(x2d) eq 0) then x2d = 0
  if (x2d eq 0) then nx2 = long(nx*ndomains/nout_avg) > 1 else begin
     nx2 = 1
     print,'Reading in y-z slices'
  endelse
  if (n_elements(y2d) eq 0) then y2d = 0
  if (y2d eq 0) then ny2 = long(ny/nout_avg) > 1 else begin
     ny2 = 1
     print,'Reading in x-z slices'
  endelse
  if (n_elements(z2d) eq 0) then z2d = 0
  if (z2d eq 0) then nz2 = long(nz/nout_avg) > 1 else begin
     nz2 = 1
     print,'Reading in x-y slices'
  endelse

  ; adjust nout for infrequent full array output when using spectral output
  if (full_array gt 0) then nout = full_array_nout

  ; change data names for spectral data
  if (dataType ge 3) then begin
     len = strlen(dataname_str)
     dist_n = strmid(dataname_str,len-1)
     if (dist_n eq '0') or (dist_n eq '1') then dat_n = strmid(dataname_str,0,len-1) $
     else begin
        dat_n = dataname_str
        dist_n = ''
     endelse
     dataname_str = dat_n + 'ft' + dist_n
  endif

  ; determine array dimensions
  if (n_elements(dims) eq 0) then dims = [nx2, ny2, nz2]
  if (dataname_base eq 'vdis') then begin
     if (n_elements(pnvy0) eq 0 or ndim_space eq 1) then pnvy0=1
     if (n_elements(pnvz0) eq 0  or ndim_space eq 1) then pnvz0=1
     dims = [pnvx0/nout_avg, pnvy0/nout_avg, pnvz0/nout_avg]
  endif

  ;; ;determine which timesteps to read in
  ;; if (n_elements(t1) eq 0) then  ntimes = get_timesteps_h5(dataname_str, dataType, $
  ;;                                get_times=get_times) else ntimes=[1,1,t1/nout]
  ;;    ;ntimes array: first element: # of data output timesteps
  ;;    ;second element: subcycle
  ;;    ;then array of indices where data is found within H5 file
  ;; asize = size(ntimes)
  ;; indices = ntimes[1:asize[1]-1]
  ;; if (n_elements(subcycle) eq 0) then subcycle = (ntimes[2]-ntimes[1])*nout

  tsize = (endt-startt)/(subcycle*skip) + 1 ;number of timesteps that will be read in
  
  ; create permanent array to store data
  data_array = fltarr(dims[0],dims[1],dims[2],tsize)
  data_size = size(data_array)

  i = 0
  for tstep = startt, endt, subcycle*skip do begin

     ; define name/location within HDF file
     time_str = strtrim(string(round(tstep)),1)

     ; do some basic error handling. attempts to return any useable data if run has been
     ; unexpectedly cut off
     catch, error_status
     if error_status eq -1008 then begin ; data doesn't exist
        print, 'Error on read-in.  Data does not exist at t=' + time_str + $
               '.  Returned array contains available data.'
        if i gt 0 then return, data_array[*,*,*,0:i-1] else return,0
        catch, /cancel
     endif else if error_status eq -170 then begin ; data is the wrong size
        print, 'Error on read in.  Data has dimensions different than specified in eppic.i.'
        return,0
        catch, /cancel
     endif else if (error_status ne 0) then begin ; something else is wrong
        print, 'Unknown error on read in.  Attempting to return any useable data.'
        if i gt 0 then return, data_array[*,*,*,0:i-1] else return,0
        catch, /cancel
     endif

    ; read in data
    data_array[*,*,*,i] = read_tstep_h5(dataname_str, time_str, dataType, ndim=ndim, $
                          ndomains=ndomains, data_size=data_size, ftout=ftout,x2d=x2d, $
                          y2d=y2d, z2d=z2d) 
    i++
  endfor

  ; correct density normalization:
  if (den_norm ne 0) then begin
     data_array+=1
     data_array*=n0d1*param6_1
  endif

  print, dataname_str + ' read-in complete' 
  return, data_array
  
end
