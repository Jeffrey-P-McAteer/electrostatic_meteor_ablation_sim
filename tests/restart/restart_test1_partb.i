;; This is an input example that should lead to a Farley-Buenamen instability
;; formation. It is testing the efield_quasineut algorithm. 

;;----------------------------------------------------------------------
;; Begin System Parameters
;;----------------------------------------------------------------------
;;--- grid ---
nx = 64
ny = 64
nsubdomains = 1
dx = 0.125
dy = 0.125

;;--- time ---
dt = 1.0e-5
nt = 4

;;--- efield ---
Ey0_external = -2.0e-2
efield_algorithm = 2
efield_petsc_algorithm = 0
;       fwidth = 3
;       fsteep = 3

;;--- misc ---
Bz = 2.5e-5
eps = 8.8542e-12

;;----------------------------------------------------------------------
;; Begin I/O Parameters
;;----------------------------------------------------------------------
nout = 1
nout_avg = 1
npout = 512
iwrite = 1
iread = 1
divj_out_subcycle = 8
outdir = "./restart_test1"

;;----------------------------------------------------------------------
;; Begin Collision Parameters
;;----------------------------------------------------------------------
vth_neutral = 249.85950
m_neutral = 4.6e-26

;;----------------------------------------------------------------------
;; Begin Distribution Definitions
;;----------------------------------------------------------------------
ndist = 2

;;----------------------------------------------------------------------
;; Begin Distribution (0) Parameters
;;----------------------------------------------------------------------
dist=0
method=-3                                                         ;fluid
md0 = 9.1094e-31
qd0 = -1.602e-19
vx0 = 0
vy0 = 0

;;--- Temp/Collisions ---
vxthd0 =5.61474e4
vythd0 =5.61474e4
coll_rate0 = 4.0e4

;;--- Initialization ---
init_dist0 = 1

;;--- I/O ---
vdist_out_subcycle0 = 2
part_out_subcycle0 = 4

;;----------------------------------------------------------------------
;; Begin Distribution (1) Parameters
;;----------------------------------------------------------------------
dist=1
method=0                                                            ;PIC
n0d1 = 1.0e11
npd1 = 921600
md1 = 4.6e-26
qd1 = 1.602e-19

;;--- Temp/Collisions ---
vxthd1 = 249.85950
vythd1 = 249.85950
coll_rate1 = 2.8e3

;;--- Initialization ---
init_dist1 = 0
pnvx1 = 64
pnvy1 = 64
pvxmin1 = -2000
pvxmax1 = 2000
pvymin1 = -2000
pvymax1 = 2000
pvzmin1 = -2000
pvzmax1 = 2000

;;--- I/O ---
vdist_out_subcycle1 = 2
part_out_subcycle1 = 4



;; THE END

