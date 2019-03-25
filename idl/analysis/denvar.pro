
; get density
@denep
@phiep
densize = size(den1,/dimensions)
nt = densize[3]
; initialize
denvar = findgen(nt)
denmax = denvar
phivar = denvar
phimax = denvar
for i=0,nt-1 do denvar[i] = sqrt(total((den1[*,*,0,i])^2)/((nx*nsubdomains+0.0)*ny*nz))
for i=0,nt-1 do phivar[i] = sqrt(total((phi[*,*,0,i])^2)/((nx*nsubdomains+0.0)*ny*nz))

for i=0,nt-1 do denmax[i] = max(den1[*,*,0,i])
for i=0,nt-1 do phimax[i] = max(phi[*,*,0,i])

;end
