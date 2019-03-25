function image_range, image
;
; NAME:
;      image_range
;
; PURPOSE:
;      To return the minimum and maximum rows and columns where image
;      is non zero.
;
; EXAMPLE:
;      image_range, (pic gt max(pic)/100)
;
;      This returns a 2x2 array with the column (x) and row (y) where
;      pic reaches at least 1% of its peak value.
;

nonzero = where(image)
column =  nonzero MOD (size(image))[1]
row= nonzero / (size(image))[1]

xy=[ [min(column), min(row)], [max(column), max(row)] ]

return, xy

end
