pro date_plot, label, charsize=charsize, _EXTRA = e

;Put the label, date and time on the lower left corner of the plot

if n_elements(label) eq 0 then label=""
if n_elements(charsize) eq 0 then charsize=0.5

xyouts, 0, 0, label + " " + systime(), /device, charsize=charsize, _EXTRA = e

end
