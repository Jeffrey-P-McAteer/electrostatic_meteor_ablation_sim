; EPPIC parameter file (also readable by IDL)
title = 'Small run to test fluid flux output'
; Mesh size and Grid size
ndim_space =        2
nsubdomains =        2
nx =       16
dx =       1.0
ny =      32
dy =       1.0
; Time step, number of time steps to run, and time-steps between output:
dt =    0.5
nt =     1
nout =        1
nout_avg =        1
; Fraction of particles to output
npout =     4096
; Define the dielectric of free space (default epsilon=8.8542e-12)
eps = 1.0
; Bz in natural units
Bz =       2.0
; Damping width: (try fwidth=3)
fwidth =        3
; Damping steepness: (try fsteep=3)
fsteep =        3
; Time steps beteen dumps
iwrite =    32000
; Start from t=0 (iread=0) or dump (iread=1)
iread =        0
; Number of distributions to initialize
ndist =        2
; Neutral parameters:
vth_neutral =       0.025
m_neutral =       400.0	
; External  Electric field in y-direction
Ey0_external =   0.001
;
; First distribution: electrons
dist = 0
method = -3
init_dist0 =        1
;method = 0
;npd0 =       204850
;init_dist0 =        3
;param3 =       1.00000
;param4 =       5.0
; number of particles, density, mass, charge, thermal velocity and drift speed
n0d0 = 1.0
; Make electrons too massive
md0 =       1.0
qd0 =      -1.00000
; Increase the coll_rate 
coll_rate0 = 0.001
; Make the mass of the neutrals 100 x lighter(for e-) 
massd_neutral0 =       4.0
vthd_neutral0 =   0.5
; Thermal velocities 
vxthd0 =       1.0
vythd0 =       1.0
vx0d0 =       1.00000
vy0d0 =       2.00000
vdist_out_subcycle0 =        8
part_out_subcycle0 =        8
flux_out_subcycle0 =        8
nvsqr_out_subcycle0 =        8
;
;Second distribution - ions - 1
dist = 1
method = 0
npd1 =       10000
md1 =       400.0
n0d1 =       1.00000
qd1 =       1.00000
coll_rate1 =    0.001
vx0d1 =       0.00000
vy0d1 =       0.00000
vxthd1 =       0.05
vythd1 =       0.05
init_dist1 =        3
param3 =       1.00000
param4 =       10.0
subcycle1 =        1
vdist_out_subcycle1 =        8
part_out_subcycle1 =        8
flux_out_subcycle1 =        8
nvsqr_out_subcycle1 =        8
;
;
;