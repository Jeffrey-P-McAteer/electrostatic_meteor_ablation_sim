; Generate x-t plot of A2 for the maximum iy.

if (n_elements(A2) le 2) then message, 'Please run envelope first'

a2size=size(A2)

xv=findgen(A2size(1))*dx*nout_avg
time=findgen(A2size(A2size(0)))*dt*nout

A2max=max(A2,i)
A2maxindex=calc_index(A2,i)
if a2size(0) eq 3 then a2mag=reform(A2(*,A2maxindex(1),*)) else $
a2mag=reform(A2(*,A2maxindex(1),A2maxindex(2),*))

label='A!E2!N: max y= '+ strcompress(string(A2maxindex(1)*dy*nout_avg))
if a2size(0) eq 4 then label=label +' z=' + $
  strcompress(string(A2maxindex(2)*dz*nout_avg))

image_plot,a2mag, xv, time, $
  /zlog, /legend, nlabels=4, $
  title=label, $
  xtitle='x (m) debye=1', $
  ytitle='time (s) wp=1'

date_plot,title
  
end
