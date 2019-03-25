{ Fill in the following sections (removing comment marks ! if necessary),
  and delete those that are unused.}
TITLE 'Efield Quasi-neutral Solver'     { the problem identification }
COORDINATES cartesian2  { coordinate system, 1D,2D,3D, etc }
VARIABLES        { system variables }
  phi            { choose your own names }

SELECT         { method controls }
aspect 2.0
curvegrid
 ngrid 64
gridarc 30.0
gridlimit 8
initgridlimit 5
 nodelimit 100000
regrid
DEFINITIONS    { parameter definitions }
{grid params}
    xmin=-37.5
    xmax=37.5
    lenx = xmax-xmin
    ymin=-37.5
    ymax=37.5
    leny = ymax-ymin
{fields}
	Bz=2.8e-5
	E0y=2.8e-3
	E0=vector(0.0,E0y)
	zdir=vector(0.0,0.0,1.0)
{temperature}
Te=660
kb=1.3806504e-23
qe=-1.6e-19
qi=abs(qe)
me=9.1e-31
mi=4.6e-26
nue=4.0e4
nui=2.5e3
omegae=4.92476e+06 !Bz*qe/me
omegai=90
kappa = omegae/nue
gamma = (1+kappa*kappa)*me*nue/abs(qe) !(me/qe^2)*(omegae^2+nue^2)/nue
L=0
{density}
n0=6e7
nPeriodsX=1
nPeriodsY=1
n = n0*(1-0.03*sin(nPeriodsX*2*pi*x/lenx)-0.01*sin(nPeriodsY*2*pi*y/lenx) + L*(y-ymin))
{Current}
Jx= n0*qi*(sin(nPeriodsX*2*pi*x/lenx)+sin(nPeriodsY*2*pi*y/lenx))*0.1
Jy= n0*qi*(sin(nPeriodsX*2*pi*x/lenx)+sin(nPeriodsY*2*pi*y/lenx))*0.1
vel = vector(Jx/(n*qi),Jy/(n*qi))
 phinot=nue*nui/(omegae*omegai)

! INITIAL VALUES
EQUATIONS        { PDE's, one for each variable }

! equation taken from Oppenheim et al 1996

  n*div(grad(phi))+dot(grad(phi),grad(n)) +kappa*(dx(phi)*dy(n)-dx(n)*dy(phi)) =Te*kb/abs(qe)*div(grad(n))+dot(E0,grad(n))-kappa*(dx(n)*E0y)+gamma/qi*div(vector(Jx,Jy))



! CONSTRAINTS    { Integral constraints }
BOUNDARIES       { The domain definition }

  REGION 1       { For each material region }

    { Periodic bottom boundary }
	!start(-1,-1)
	start(xmin,ymin)
	!periodic(x,y+2)  line to(0.95,-1)
	periodic(x,y+leny)  line to(0.95*xmax,ymin)
      { New "line" spec breaks periocity }
	!   optional:   natural(u) = normal(K*grad(u))
	!line to (1,-1)
 	line to (xmax,ymin)
      { Periodic right boundary }
      	!periodic(x-2,y) line to (1,1)
	periodic(x-lenx,y) line to (xmax,ymax)

      { Images of non-periodic stub and periodic bottom boundry }
	!line to (0.95,1) to (-1,1)
	line to (0.95*xmax,ymax) to (xmin,ymax)

      { Image of periodic right boundary }
      	line  to close



 ! TIME 0 TO 1    { if time dependent }
MONITORS         { show progress }
PLOTS            { save result displays }
  vector(((-grad(phi)+E0)/Bz)/((1+phinot)*sqrt(1e5))) fixed range (0,1)
!vector(((E0)/Bz)/((1+phinot)*sqrt(1e5))) fixed range (0,1)
  contour(n) painted export (65) file "density.out" export format "#x#b#y#b#1"  noheader
  contour(jx) painted export (65) file "jx.out" export format "#x#b#y#b#1"  noheader
  contour(jy) painted export (65) file "jy.out" export format "#x#b#y#b#1"  noheader
  contour(phi) painted export (65) file "phi.out" export format "#x#b#y#b#1"  noheader
  Summary as "parameters"
	report E0y/Bz as "Drift Speed/ m*s-1"
	report gamma
	report omegae
	report n0
	report kappa
END

