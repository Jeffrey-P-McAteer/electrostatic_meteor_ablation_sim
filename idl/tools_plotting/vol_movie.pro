pro vol_movie,phi,angx=angx,angz=angz,threshold=threshold,LOW=low, $
              relative=relative,skip=iskip,winsize=winsize

; phi is a 4D array (x,y,z,time)
; angx,angz are rotation angles around x and z (default: 30 and 50)
; threshold is relative threshold (see code) (0.-1.)
; skip skips every (iskip-1) image (default iskip=1 - no skipping)
; winsize sets the window size (default 600x600)
; there is an option called LOW that might be useful but I have not
; included it in this argument list.
; This calls on pro cube in order to draw a coordinate cube

svol = SIZE(phi);	Get the dimensions of the volume.

fx=1 & fy=1 & fz=1
vx=svol(1)/fx 
vy=svol(2)/fy
vz=svol(3)/fz

if n_elements(angx) eq 0 then angx=30 
if n_elements(angz) eq 0 then angz=50
if n_elements(threshold) eq 0 then threshold=.001
if n_elements(iskip) eq 0 then iskip=1
if n_elements(winsize) ne 2 then winsize=[600,600]

mxv=max(phi)
mnv=min(phi)
thresh=threshold*(mxv-mnv)

lx=0.8*vx & ly=0.8*vy & lz=0.8*vz


create_view,ax=angx,az=angz,xmax=lx,ymax=ly,zmax=lz ,xmin=-lx,ymin=-ly,zmin=-lz

window,2,xsize=winsize(0),ysize=winsize(1)
xinteranimate,set=[long(winsize),svol(4)]

for time=0,svol(4)-1, iskip do begin
    print,'time=',time
    vol=phi(*,*,*,time)
    s=size(vol)
    
    mxv=max(vol) & mnv=min(vol)
;    print,' min and max of volume are ',mnv,' and ',mxv

    if keyword_set(relative) then thresh=threshold*(mxv-mnv)

    IF s(0) NE 3 THEN print,'this is not a volume array'
;   Use SCALE3 to establish the 3D transformation and coordinate
;   ranges.
    SCALE3, XRANGE=[0, S(1)], YRANGE=[0, S(2)], ZRANGE=[0, S(3)],ax=angx,az=angz

;   Default = view high side of contour surface.
    IF N_ELEMENTS(low) EQ 0 THEN low = 0

    SHADE_VOLUME, vol, thresh, v, p, LOW = low

;   Produce vertices and polygons.
    if (n_elements(p) gt 0) then begin
        TVscl, POLYSHADE(v,p,/data,/T3D) ;Produce image of surface and display.
    endif else $
      print, 'VOL_MOVIE skipping element ',time,' - threshold not reached'

    cube,s(1),s(2),s(3)

    xinteranimate,frame=time,window=2

endfor    

xinteranimate

END


