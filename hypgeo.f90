FUNCTION hypgeo(a,b,c,z)
  
  implicit none
  
COMPLEX hypgeo,a,b,c,z
REAL EPS
PARAMETER (EPS=1.e-6) !Accuracy parameter.
! !..C USES bsstep,hypdrv,hypser,odeint
! Complex hypergeometric function 2F1 for complex a, b, c, and z, by direct integration of
! the hypergeometric equation in the complex plane. The branch cut is taken to lie along
! the real axis, Re z > 1.
INTEGER kmax,nbad,nok
EXTERNAL bsstep,hypdrv
COMPLEX z0,dz,aa,bb,cc,y(2)
COMMON /hypg/ aa,bb,cc,z0,dz
COMMON /path/ kmax !..Used by odeint.
kmax=0
if (real(z)**2+aimag(z)**2.le.0.25) then !..Use series...
call hypser(a,b,c,z,hypgeo,y(2))
return
else if (real(z).lt.0.) then !...or pick a starting point for the path integration
z0=cmplx(-0.5,0.)
else if (real(z).le.1.0) then
z0=cmplx(0.5,0.)
else
z0=cmplx(0.,sign(0.5,aimag(z)))
endif
aa=a !..Load the common block, used to pass parameters “over the head” of odeint to hypdrv.
bb=b
cc=c
dz=z-z0
call hypser(aa,bb,cc,z0,y(1),y(2)) !..Get starting function and derivative.
call odeint(y,4,0.,1.,EPS,.1,.0001,nok,nbad,hypdrv,bsstep)
! The arguments to odeint are the vector of independent variables, its length, the starting and
! ending values of the dependent variable, the accuracy parameter, an initial guess for stepsize,
! a minimum stepsize, the (returned) number of good and bad steps taken, and the names of
! the derivative routine and the (here Bulirsch-Stoer) stepping routine.
hypgeo=y(1)
return
END FUNCTION hypgeo

