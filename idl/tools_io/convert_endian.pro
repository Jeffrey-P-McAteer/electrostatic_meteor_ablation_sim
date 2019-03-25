FUNCTION convert_endian, vec, reorder_bytes=reorder_bytes
; Reverse the byte ordering if keyword is set or if the data looks
; scrambled and this is an alpha

if (keyword_set(reorder_bytes)) then byteorder, vec, /ftoxdr $
else if (!version.arch eq 'alpha') then $
  if (not (max(vec) lt 1e18 and min(vec) gt -1e18)) then begin
    byteorder, vec, /ftoxdr
                                ; test to see if the numbers are better
    if (not (max(vec) lt 1e18 and min(vec) gt -1e18)) then $
      byteorder, vec, /ftoxdr $
    else reorder_bytes=1
endif
return, vec
end
