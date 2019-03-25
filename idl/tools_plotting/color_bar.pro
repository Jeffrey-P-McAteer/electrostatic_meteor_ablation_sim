PRO COLOR_BAR,THICK=THICK,ROTATION=ROTATION,XPOS=XPOS,YPOS=YPOS,LENGTH=LENGTH,UPPER=UPPER,LOWER=LOWER
;+
;
; NAME:
;	COLOR_BAR
; PURPOSE:
;	To prvide a colour bar on the graphics device.
; CATEGORY:
;	Displays
; CALLING SEQUENCE:
;	COLOR_BAR [,THICK=THICK,ROTATION=ROTATION,XPOS=XPOS,YPOS=YPOS,LENGTH=LENGTH,UPPER=UPPER,LOWER=LOWER]
; INPUTS:
;	NONE
; OPTIONAL INPUT PARAMETERS:
;	NONE
; OUTPUTS:
;	To window or to PS file.
; OPTIONAL OUTPUT PARAMETERS:
;	NONE
; KEYWORDS:
;	THICK = a factor for multiplying the thickness of the bar
;               (the default size is 5%  of the plot size)
;	LENGTH =  a multiplication factor controlling the length (the
;                 default length is the plot size)
;	ROTATION = controls whether the colourbar lies down is vertical etc (default is lying down left to right).  
;	XPOS,YPOS = controls the position in the window (default = 0,0)
;	UPPER and LOWER specify the maximum and minimum byte values for the colour bar (in case
;		you are reserving other colours for overlays etc.) (default =255,0)
; COMMON BLOCKS:
;	NONE
; SIDE EFFECTS:
;	NONE
; RESTRICTIONS:
;	NONE
; PROCEDURE:
;	TRIVIAL
; EXAMPLE:	
;	COLOR_BAR,THICK=1.,ROTATION=1,XPOS=200,YPOS=0,LENGTH=1
;
; MODIFICATION HISTORY:
;	31-1-1993 NT
;       13-3-1997 Meers Oppenheim
;-


STATUS=STRARR(3)
STATUS(0)='S'
STATUS(1)='SUCCESS'
STATUS(2)='COLOR_BAR'
VERSION='V1.0'


IF (N_PARAMS(0) LT 0) OR (N_PARAMS(0) GT 0) THEN BEGIN
PRINT,' '
PRINT,STRTRIM(STRUPCASE(STATUS(2)),2)+' :  BAD INPUT PARAMETERS.  EXIT.'
STATUS(0)='F'
STATUS(1)='FAILURE'
GOTO,CLOSING
ENDIF

IF N_ELEMENTS(LENGTH) EQ 0 THEN LENGTH=1
IF N_ELEMENTS(THICK) EQ 0 THEN THICK=1
IF N_ELEMENTS(ROTATION) EQ 0 THEN ROTATION=0
IF N_ELEMENTS(XPOS) EQ 0 THEN XPOS=0
IF N_ELEMENTS(YPOS) EQ 0 THEN YPOS=0
IF N_ELEMENTS(UPPER) EQ 0 THEN UPPER=!D.TABLE_SIZE
IF N_ELEMENTS(LOWER) EQ 0 THEN LOWER=0

IF UPPER LT LOWER THEN BEGIN
PRINT,' '
PRINT,STRTRIM(STRUPCASE(STATUS(2)),2)+' :  BAD INPUT UPPER AND LOWER VALUES.  EXIT.'
STATUS(0)='F'
STATUS(1)='FAILURE'
GOTO,CLOSING
ENDIF

ROT=FIX(ROTATION)
IF (ROT LT 0) OR (ROT GT 3) THEN BEGIN
PRINT,' '
PRINT,STRTRIM(STRUPCASE(STATUS(2)),2)+' :  ROTATION MUST BE AN INTEGER >0<3 .  EXIT.'
STATUS(0)='F'
STATUS(1)='FAILURE'
GOTO,CLOSING
ENDIF

xwin=!x.window*!d.x_size
ywin=!y.window*!d.y_size


;Define array of colors
A=BYTE(FINDGEN(UPPER-LOWER,2) MOD (UPPER-LOWER) ) 
IF ROT GE 2 THEN A=REVERSE(A)

; If not scalable pixels:
if (!d.flags and 1) eq 0 then begin
    if rot eq 0 or 2 then $
      A=CONGRID(A,LENGTH*(xwin(1)-xwin(0)),THICK*(ywin(1)-ywin(0))/20.) $
    else $
      A=CONGRID(A,LENGTH*(ywin(1)-ywin(0)),THICK*(xwin(1)-xwin(0))/20.)
endif

IF ROT EQ 0 OR ROT EQ 2 THEN $
  TV,A,XPOS,YPOS, /device, $
  XSIZE=LENGTH*(xwin(1)-xwin(0)),YSIZE=THICK*(ywin(1)-ywin(0))/20. $
ELSE BEGIN
    A=TRANSPOSE(A)
    TV,A,XPOS,YPOS, /device, $
      YSIZE=LENGTH*(ywin(1)-ywin(0)),XSIZE=THICK*(xwin(1)-xwin(0))/20.
ENDELSE

CLOSING:

RETURN
END
