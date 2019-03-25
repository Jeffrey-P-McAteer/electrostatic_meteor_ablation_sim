;; This is an input example that should lead to a Farley-Buenamen instability
;; formation. It is testing the efield_quasineut algorithm. 
;; The parameters match those used in Meers's Disesrtation.

;;----------------------------------------------------------------------
;; Begin System Parameters
;;----------------------------------------------------------------------
;;--- grid ---
nx = 64
ny = 256
nsubdomains = 4
dx = 0.08
dy = 0.08

;;--- time ---
dt = 1.0e-6
nt = 50000

;;--- efield ---
boundary_type =1
Ex0_external = -2.0e-3
efield_algorithm = 0
efield_petsc_algorithm = 0
fwidth = 0
fsteep = 0

;;--- misc ---
Bz = 2.5e-5
eps = 8.8542e-12

;;----------------------------------------------------------------------
;; Begin I/O Parameters
;;----------------------------------------------------------------------
nout = 100
nout_avg = 1
npout = 512
iwrite = 50000
iread = 0
divj_out_subcycle = 1
outdir = "./non_periodic_data"

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
method=0
n0d0 = 1.0e8
n0lhsd0 = 0.999e8
n0rhsd0 = 1.0e8
npd0 = 250000
md0 = 9.1094e-31
qd0 = -1.602e-19
vx0d0 = 0
vy0d0 = 0

;;--- Temp/Collisions ---
vxthd0 = 5.61474e4
vythd0 = 5.61474e4
coll_rate0 = 4.0e4

;;--- Initialization ---
init_dist0 = 9

;;--- I/O ---
vdist_out_subcycle0 = 1
part_out_subcycle0 = 1

;;----------------------------------------------------------------------
;; Begin Distribution (1) Parameters
;;----------------------------------------------------------------------
dist=1
method=0                                                            ;PIC
n0d1 = 1.0e8
n0lhsd1 = 0.999e8
n0rhsd1 = 1.0e8
npd1 = 10000000
md1 = 4.6e-26
qd1 = 1.602e-19

;;--- Temp/Collisions ---
vxthd1 = 249.85950
vythd1 = 249.85950
coll_rate1 = 2.8e3
subcycle = 100

;;--- Initialization ---
init_dist1 = 9
pnvx1 = 64
pnvy1 = 64
pvxmin1 = -2000
pvxmax1 = 2000
pvymin1 = -2000
pvymax1 = 2000
pvzmin1 = -2000
pvzmax1 = 2000

;;--- I/O ---
vdist_out_subcycle1 = 1
part_out_subcycle1 = 1

;; THE END

