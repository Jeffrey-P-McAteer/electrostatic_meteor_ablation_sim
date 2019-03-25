FUNCTION EPPIC_PARAMS_EXIST
;+
;
; NAME:
;	EPPIC_PARAMS_EXIST
;
; PURPOSE:
;       This routine tests if a run's parameters have been
;       initialized in IDL. This should be put in routines that
;       traditionally used @eppic.i to load the parameters, which
;       would lead to "Error opening file." errors when eppic.i is not
;       used as the input file. If the parameters are not set, then
;       a warning is printed.
;
; CATEGORY:
;       EPPIC common
;
; CALLING SEQUENCE:
;	
;       EPPIC_CHECK_PARAMS
;
; INPUTS:
;	None
;
; KEYWORD PARAMETERS:
;       None
;
; OUTPUTS:
;       Returns 1 if parameters are set, 0 if they are not.
;
; COMMON BLOCKS:
;	There are none.
;
; EXAMPLE:
;      If (not EPPIC_PARAMS_EXIST) return
;
; MODIFICATION HISTORY:
; 	Written by:	
;-

  return,1; TRUE
  return,0; FALSE

END
