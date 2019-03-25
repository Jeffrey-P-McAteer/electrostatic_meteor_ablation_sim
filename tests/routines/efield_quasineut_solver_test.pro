; compile require routines, incase they've changed
.compile from_flexpde
.compile readarray
.compile writearray
.compile read_domains
.compile write_domains

;---------------FLEX PDE DATA---------------------;
; read in FlexPDE's initial conditions
from_flexpde,"density.out",den_flexpde
den0=den_flexpde[0:63,0:63] ; while periodic, den output adds extra row, col
from_flexpde,"jx.out",jx_flexpde
jx0=jx_flexpde[0:63,0:63] ; while periodic, den output adds extra row, col
from_flexpde,"jy.out",jy_flexpde
jy0=jy_flexpde[0:63,0:63] ; while periodic, den output adds extra row, col
from_flexpde,"phi.out",phi_flexpde
phi0=phi_flexpde[0:63,0:63] ; while periodic, den output adds extra row, col
                                phi_avg = total(phi0)/(64*64)
for i=0,63 do for j=0,63 do phi0[i,j] = phi0[i,j]-phi_avg

;---------------IDL DATA--------------------------;
;; den0=findgen(64,64)
;; PI=3.14
;; for i=0,63 do for j=0,63 do den0[i,j]=1e11*(1.0-0.03*sin(2.0*PI*i/64.0)-0.01*sin(2.0*PI*j/64.0))

;; jx0 = findgen(64,64)*0+1.0
;; jy0 = findgen(64,64)*0+1.0

; FLEXPDE overwrite

;for i=0,63 do for j=0,63 do den0[i,j]=den0[i,j]
;for i=0,63 do for j=0,63 do jx0[i,j]=jx_flexpde[i,j]
;for i=0,63 do for j=0,63 do jy0[i,j]=jy_flexpde[i,j]

; write them to binary format to feed the test program
writearray,"bin_density.out",den0,order=[1,0],/binary
writearray,"bin_jx.out",jx0,order=[1,0],/binary
writearray,"bin_jy.out",jy0,order=[1,0],/binary
write_domains,"bin_density_domain*.out",den0,2,order=[1,0],/binary,basepath="./"
write_domains,"bin_jx_domain*.out",jx0,2,order=[1,0],/binary,basepath="./"
write_domains,"bin_jy_domain*.out",jy0,2,order=[1,0],/binary,basepath="./"

; run the test program
$make efield_quasineut_solver_test
;$./efield_quasineut_solver_test
$mpirun -np 2 ./efield_quasineut_solver_test

; get the data from the test program
;den_eppic = readarray("bin_density.in",[64,64,1],order=[2,1,0],/binary)
den_ranged_eppic = read_domains("bin_density_domain*.in",[32,64,1],ndomains=2,order=[2,1,0],/binary,basepath="./")
phi_ranged_eppic = read_domains("bin_phi_domain*.in",[32,64,1],ndomains=2,order=[2,1,0],/binary,basepath="./")
;phi_eppic = readarray("bin_phi.out",[64,64,1],order=[2,1,0],/binary)
;phi = phi_eppic[2:65,0:63]
; for i=0,63 do for j=0,63 do if (i=0) then $
;    phi[i,j] = phi_eppic[65,j] $
;for i=0,63 do for j=0,63 do phitmp[i,j] = phi[63-i,63-j]
;phi = phitmp

; process and plot
;window,0
;image_plot,[[den0,den0],[den0,den0]],/legend
;image_plot,phi0,/legend
;image_plot,[[phi0,phi0],[phi0,phi0]],/legend
;window,1
;image_plot,[[den_eppic,den_eppic],[den_eppic,den_eppic]],/legend
;image_plot,phi,/legend
;image_plot,[[phi,phi],[phi,phi]],/legend
;window,2
;image_plot,abs((phi-phi0)/phi),/legend,/zlog
;
; else  phi[i,j] = phi_eppic[3+i,j] 
      
;phitmp = phi
;for i=0,63 do for j=0,63 do phitmp[i,j] = phi[63-i,63-j]
;phi = phitmp

; FROM EPPIC 
window,0
image_plot,phi_ranged_eppic,/legend
;image_plot,den_ranged_eppic,/legend

; FROM FLEX PDE
window,1
image_plot,phi0,/legend
;image_plot,den0,/legend

; DIFFERENCE
window,2
image_plot,abs(phi_ranged_eppic-phi0)/max(phi0),/legend,/zlog




