;; A small program to workup the results of the fluid_flux.i test


threshold=1.0E-7
status=0
@fluid_flux.i
@denep
@fluxep
;.run ../../idl/domains/denep.pro
;.run ../../idl/domains/fluxep.pro

den0=n0d0*(den0+1)
difference=total(fluxx0/den0)/nx/ny/nsubdomains-vx0d0
print,"Difference x: ",difference
if (abs(difference) GT threshold) then status=1
difference=total(fluxy0/den0)/nx/ny/nsubdomains-vy0d0
print,"Difference y: ",difference
if (abs(difference) GT threshold) then status=1

exit,status=status
