pro plot_h5, filename_str, data_array=data_array, savename_str=savename_str, step=step, $
             ptitle=ptitle, xtitle=xtitle,ytitle=ytitle, ztitle=ztitle, xy=xy,xz=xz,yz=yz,$
             ftout=ftout,startt=startt,endt=endt,dataType=dataType
;
; NAME:
;       PLOT_H5
;
; PURPOSE:
;       Create plots of HDF5 data
;
; CATEGORY:
;       Eppic I/O
;
; CALLING SEQUENCE:
;       FT_PLOT_H5,'filename_str'[,'savename_str',step=step,xtitle='xtitle',ytitle='ytitle']
;
; INPUTS:
;       filename_str:   The name of the dataset containing the data
;       
;       savename_str:   The name of the PostScript file the plots will be saved to
;                       Defaults to 'filename_str_image'
;       
;       step:           The fraction of timesteps to be plotted
;                       Defaults to 10
;       
;     x/y/z title:      Titles of the plot axes
;                       Defaults to 'x' 'y' 'z'

@eppic.i
  ndim = ndim_space
  if (n_elements(xy) eq 0) && (n_elements(yz) eq 0) && (n_elements(xz) eq 0) then xy=1
  if not keyword_set(step)   then step = 1;0
  if not keyword_set(savename_str) then savename_str = filename_str + '_image'
  if not keyword_set(ptitle) then ptitle = filename_str + ' t='
  if ft_output_arrays eq 1 then nout_avg = 1

  xvec = findgen(nx*nsubdomains/nout_avg)*dx*nout_avg
  xtitle = 'x (m)'
  yvec = findgen(ny/nout_avg)*dy *nout_avg
  ytitle = 'y (m)'
  if (ndim eq 3) then begin
     zvec = findgen(nz/nout_avg)*dy*nout_avg
     ztitle = 'z (m)'
  endif

  if (n_elements(data_array) eq 0) then data_array=read_h5(filename_str, ftout=ftout,skip=step,$
                                             startt=startt,endt=endt,dataType=dataType,/den_norm) 
  dsize = size(data_array)
  ntimes = n_elements(data_array[0,0,0,*])

  ps,savename_str + '.ps', /LANDSCAPE, /COLOR
  !P.MULTI = [0,2,2,0,0]

  case ndim of
     2: for tstep = 0, ntimes[0]-1 do image_plot,data_array[*,*,0,tstep],xvec,yvec,skip=step, $
           xtitle=xtitle,ytitle=ytitle,title=ptitle+string(tstep*nout*dt),/aspect,/legend
     3: begin
        for tstep = 0, ntimes[0]-1 do image_plot,data_array[*,*,(dsize[3])/2,tstep],xvec,yvec, $
           skip=step,xtitle=xtitle, ytitle=ytitle,title=ptitle+string(tstep*nout*dt),/aspect,/legend
        for tstep = 0, ntimes[0]-1 do image_plot,data_array[*,(dsize[2])/2,*,tstep],xvec,zvec, $
           skip=step,xtitle=xtitle,ytitle=ztitle,title=ptitle+string(tstep*nout*dt),/aspect,/legend
        for tstep = 0, ntimes[0]-1 do image_plot,data_array[(dsize[1])/2,*,*,tstep],yvec,zvec, $
           skip=step,xtitle=ytitle,ytitle=ztitle,title=ptitle+string(tstep*nout*dt),/aspect,/legend
        end
  endcase

  stop_ps
  print, filename_str + ' plotted'

end
