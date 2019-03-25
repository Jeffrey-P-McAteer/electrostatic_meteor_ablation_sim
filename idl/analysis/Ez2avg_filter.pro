;  Routine to produce plots of filtered E2avg

@phiep
.r E2avg3d

Ex2avg0=Ex2avg
E2avg0=E2avg

; Now filter phi and repeat

fwidth=6
fsteep=3
.r E2avg3d
Ex2avg6=Ex2avg
E2avg6=E2avg
fwidth=9
.r E2avg3d


Ex2avg9=Ex2avg
E2avg9=E2avg
save,Ex2avg0,E2avg0,Ex2avg6,E2avg6,Ex2avg9,E2avg9,filename='EavgE0105.pro'

ps,'Epar2avg_filter.ps'
plot,tv,Ex2avg0/E2avg0,  title='Ex!E2!N/E!E2!N Avg',  xtitle='time (s)', ytitle='Ex!E2!N/E!E2!N',/ylog,yrange=[.001,1]*max(Ex2avg0/E2avg0)

oplot,tv,Ex2avg6/E2avg6
oplot,tv,Ex2avg9/E2avg9
;Linear scale
plot,tv,sqrt(Ex2avg0/E2avg0),  title='<Ex>/<E> Avg',  xtitle='time (s)', ytitle='<Ex>/<E>'

oplot,tv,sqrt(Ex2avg6/E2avg6)
oplot,tv,sqrt(Ex2avg9/E2avg9)

date_plot,title
stop_ps
