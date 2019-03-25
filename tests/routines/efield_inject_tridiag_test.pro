; compile require routines, incase they've changed
.compile readarray
.compile writearray
.compile read_domains
.compile write_domains
.compile laplacian

eps = 8.8542e-12


;loadct,1
;------------- READ IN PHI/RHO DATA --------------;
nsubdomains=2
nx=48/nsubdomains
ny=24
nz=24
cd,"~/src/build_dirs/yann_experimental_3d/tests/routines/"
$pwd
print,"nx=",nx," ny=",ny," nx=",nx," nsubdomains=",nsubdomains

!p.multi=[0,2,2]
rho_inject = read_domains("efield_inject_tridiag_rho*.out",[nx,ny,nz],ndomains=nsubdomains,order=[2,1,0],/binary,basepath="./")
if (n_elements(rho_inject) gt 1) then $
  image_plot,rho_inject,/legend,TITLE='RHO INJECT',/interp
phi_inject = read_domains("efield_inject_tridiag_phi*.out",[nx,ny,nz],ndomains=nsubdomains,order=[2,1,0],/binary,basepath="./") ;,non_periodic_x=[1,1])
if (n_elements(phi_inject) gt 1) then $
  image_plot,phi_inject,/legend,TITLE='PHI INJECT',/interp
;; rho_spect = read_domains("efield_test_rho*.out",[nx,ny,1],ndomains=nsubdomains,order=[2,1,0],/binary,basepath="./")
;; image_plot,rho_spect,/legend,TITLE='RHO SPECTRAL',/interp
;; phi_spect = read_domains("efield_test_phi*.out",[nx,ny,1],ndomains=nsubdomains,order=[2,1,0],/binary,basepath="./");,non_periodic_x=[1,1])
;; image_plot,phi_spect,/legend,TITLE='PHI SPECTRAL',/interp
cd,"~/src/yann_experimental/tests/routines"
$pwd
