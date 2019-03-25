;Make density movie using image_plot_mpeg

if (n_elements(den0) le 2) then begin
@den3d
endif

nden1=size(den1)

image_plot_mpeg, reform(den1(*,*,0,*)), xv,yv,/legend,$
  title='Ion Density' , $
  xtitle='ExB Drift Direction',ytitle='Eo Direction' 

end
