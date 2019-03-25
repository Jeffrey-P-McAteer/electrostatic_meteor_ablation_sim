FUNCTION int_theta, z, x, y, GRID = grid, POINTS = points, $
	QUINTIC = quintic, RECALL = recall, BOUNDS = bounds
;+
; NAME:
;	INT_THETA
;
; PURPOSE:
;	This function interpolates a surface from rectangular coordinates
;       (X, Y, Z) to polar coordinates (R, Theta, Z) and integrates along
;	Theta.
;
; CATEGORY:
;	Gridding.
;
; CALLING SEQUENCE:
;	Result = INT_THETA(Z, X, Y)
;
; INPUTS:
;	Z:	 An array containing the surface value at each point.
;		 If the data are regularly gridded (GRID=1) in X and 
;		 Y, Z is a two dimensional array, where Z[i,j] has a
;		 x-coordinate of X[i] and an y-coordinate of Y[j].
;                If GRID is not set, X[i] and Y[i] contain the coordinates
;		 of each Z[i].
;	X:	 The x-coordinate. If GRID is set, Z[i,j] has a
;                x-coordinate of X[i]. If GRID is not set, R must have
;                the same number of elements as Z, and contains the
;		 x-coordinate of each point.
;	Y:	 The y-coordinate. If GRID is set, Z[i,j] has a
;                y-coordinate of Y[i]. If GRID is not set, R must have
;                the same number of elements as Z, and contains the
;		 y-coordinate of each point.
;
; KEYWORD PARAMETERS:
;	GRID:    Set GRID to indicate that Z is regularly gridded in
;		 X and Y.
;	POINTS:  The number of elements in the output array.
;	QUINTIC: If set, the function uses quintic interpolation, which is
;		 slower but smoother than the default linear interpolation.
;	RECALL:	 If set, reuse the coordinate transform from the last call.
;	BOUNDS:	 Set this keyword to four-element vector which returns
;		 the bounds of the returned array
;		 [Radius0, Theta0, Radius1, Theta1]
;
; OUTPUTS:
;	This function returns a one-dimensional array of the same type as Z.
;
; COMMON BLOCKS:
;	int_theta_common: 	Keeps the triangulation, to be reused
;				when RECALL is set.
; PROCEDURE:
;	First, each data point is transformed to (X, Y, Z). Then
;	the TRIANGULATE and TRIGRID procedures are used to interpolate
;	the surface over the polar grid.
;
; EXAMPLE:
;	n = 100
;	x = FINDGEN(n) - n/2
;	y = FINDGEN(n) - n/2
;	z = DIST(n)
;	PLOT, int_theta(z, x, y, /grid, points=n)
;
; MODIFICATION HISTORY:
;	AvE, 	June 1999, 	Written
;-

common int_theta_common, $
	nx, ny, np, theta, radius, binsize

; if problems occur, return to prompt
on_error, 1

; check the input for consistency
xsize = size(x)
ysize = size(y)
zsize = size(z)
IF (KEYWORD_SET(grid)) THEN BEGIN
    IF (zsize[0] NE 2) THEN MESSAGE, $
			'Z must be a two dimensional array, if GRID is set'
    IF (xsize[0] NE 1) THEN MESSAGE, $
			'X must be a one dimensional array, if GRID is set'
    IF (ysize[0] NE 1) THEN MESSAGE, $
			'Y must be a one dimensional array, if GRID is set'
    IF (xsize[1] NE zsize[1]) THEN MESSAGE, $
			'First dimension of Z must match dimension of X'
    IF (ysize[1] NE zsize[2]) THEN MESSAGE, $
			'Second dimension of Z must match dimension of Y'
ENDIF ELSE BEGIN
    IF (zsize[0] NE 1) THEN MESSAGE, 'Z must be a one dimensional array'
    IF (xsize[0] NE 1) THEN MESSAGE, 'X must be a one dimensional array'
    IF (ysize[0] NE 1) THEN MESSAGE, 'Y must be a one dimensional array'
    IF (xsize[1] NE zsize[1] OR ysize[1] NE zsize[1]) THEN MESSAGE, $
			'X, Y and Z must have the same dimensions'
ENDELSE


IF (KEYWORD_SET(grid)) THEN BEGIN
    IF (KEYWORD_SET(recall)) THEN BEGIN
	; check array sizes
	IF (xsize[1] NE nx OR ysize[1] NE ny) THEN MESSAGE, $
			'Input array size must not be changed if RECALL is set'
	IF KEYWORD_SET(points) THEN MESSAGE, $
		        'POINTS must not be set if RECALL is set'
    ENDIF ELSE BEGIN
	IF NOT KEYWORD_SET(points) THEN BEGIN
	    np = 50
	ENDIF ELSE BEGIN
	    IF N_ELEMENTS(points) NE 1 THEN MESSAGE, $
		   'POINTS must be a scalar'
	    np = points
	ENDELSE
	nx = xsize[1]
	ny = ysize[1]
	; triangulate the mesh
	xa = reform( rebin(          x , nx, ny, /sample), nx * ny)
	ya = reform( rebin(transpose(y), nx, ny, /sample), nx * ny)
	pol_coord=cv_coord(from_rect=transpose([ [[xa]],[[ya]] ], [2,0,1]), $
		           /to_polar)
	theta  = reform(pol_coord[0,*], nx, ny)
	radius = reform(pol_coord[1,*], nx, ny)
	binsize = (max(radius) - min(radius)) / np
    ENDELSE
ENDIF ELSE BEGIN
;--- needs to be completed ----
    MESSAGE, 'irregular gridding not supported yet'
ENDELSE

r0 = min(radius)
r1 = max(radius)
t0 = min(theta)
t1 = max(theta)

IF KEYWORD_SET(bounds) then bounds = [r0, t0, r1, t1]

result = Fltarr(np)
FOR i = 0, np-1 DO BEGIN
    index = Where(radius GE i*binsize+r0 AND radius LE (i+1)*binsize+r0, count)
    IF count GT 0 THEN result[i] = total(z[index]) / count * 2.0*!Pi $
		  ELSE result[i] = 0.0
ENDFOR



RETURN, result

END
