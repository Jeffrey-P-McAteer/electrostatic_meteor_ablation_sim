@eppic.i

; Set defaults (read in all density data & renormalize)
if (n_elements(startt) eq 0) then startt=0
if (n_elements(endt) eq 0) then endt=0
if (n_elements(start_frac) eq 0) then start_frac=0
if (n_elements(end_frac) eq 0) then end_frac=1
if (n_elements(skip) eq 0) then skip=1
if (n_elements(den_norm) eq 0) then den_norm=1
if (n_elements(subcycle) eq 0) then subcycle=nout

den0 = read_h5('den0',startt=startt,endt=endt,start_frac=start_frac,end_frac=end_frac,skip=skip, $
              den_norm=den_norm)

den1 = read_h5('den1',startt=startt,endt=endt,start_frac=start_frac,end_frac=end_frac,skip=skip, $
              den_norm=den_norm)

