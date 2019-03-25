
PRO PWD
  
;+
;
; NAME:
;	PWD
; PURPOSE:
;	To emulate the UNIX print path (pwd) command.
; CATEGORY:
;	Screwing around 
; CALLING SEQUENCE:
;	PWD
; INPUTS:
;	NONE
; OPTIONAL INPUT PARAMETERS:
;	NONE
; OUTPUTS:
;	NONE
; OPTIONAL OUTPUT PARAMETERS:
;	NONE
; KEYWORDS:
;	NONE
; COMMON BLOCKS:
;	NONE
; SIDE EFFECTS:
;	NONE
; RESTRICTIONS:
;	NONE
; EXAMPLE:
;	
; PROCEDURE:
;	Use the CD command.
;
; MODIFICATION HISTORY:
;	20-11-1994 NT
;
;----------------------------------------------------------------------------
;		        D I S R S O F T    (C)  1 9 9 4      
;-----------------------------------------------------------------------------
;-

cd,curr=curr&print,curr
 
end
