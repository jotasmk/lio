!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%%%%%%%%%%%%%%    Effective Core Potential Module    %%%%%%%%%%%%%%%!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
! Contains commmon variables and functions for ECP subroutines       !
!                                                                    !
! V0.9 sept 2015, first functional version for energy calculations   !
!                                                                    !
! Nicolas Foglia                                                     !
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!


      MODULE ECP_mod
      IMPLICIT NONE

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%%%%%%%%%%%%%%%%%%%%%    Namelist Variables    %%%%%%%%%%%%%%%%%%%%%!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!

      logical :: ecpmode
!activa los pseudopotenciales
      integer :: ecptypes
!cantidad de atomos con ECP
      character (len=30) :: tipeECP
!tipo de ECP usado, tiene que estar en $LIOHOME/dat/ECP/
      integer, dimension(128) :: ZlistECP
!Z de atomos con ECP
      logical :: cutECP
!activa cuts en las integrales de ECP
      double precision :: cutecp2, cutecp3, distcutECP, distcutECP_bohr2, cut2_0
!valores de corte para las integrales de 2 y 3 centros (AAB y BAC)
      double precision , dimension (128) :: Raction
!arrays para cut 2.0

      logical :: FOCK_ECP_read, FOCK_ECP_write
!activan lectura y grabado de Fock para evitar recalcularlo en
!dinamicas a nucleo fijo
      logical :: Fulltimer_ECP
!activa los timers para int radiales
      double precision :: tlocal,tsemilocal,tQ1,tQ2,Tiempo, Taux
!auxiiares para timers

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%%%%%%%%%%%%%%%%%%    Dbug & Verbose Variables    %%%%%%%%%%%%%%%%%%!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
      logical :: ecp_debug,ecp_full_range_int
! ecp_debug activa el modo de debugueo, ecp_full_range_int activa el
! calculo de integrales radiales en todo el rango disponible por los
! arrays
      integer :: local_nonlocal
! =1 solo calcula terminos locales <xi|V|xj>
! =2 solo calcula terminos no locales <xi|Ylm>V<Ylm|xj>
! default =0
      integer :: verbose_ECP
! controla la impresion de variables
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!



!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%%%%%%%%%%%%%     Effective Core potential Data    %%%%%%%%%%%%%%%%%!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!

        integer, dimension(118) :: Zcore, Lmax
! Zcore(Z) carga del core para el ECP elegido del atomo con carga
! nuclear Z, Lmax(Z) L maximo del ECP elegido para el atomo con
! carga nuclear Z
        integer, dimension(118,0:5) :: expnumbersECP
!expnumbersECP(Z,l) cantidad de terminos del ECP para el atomo con
! carga nuclear Z y l del ECP
        integer, dimension(118,0:5,10) :: nECP
	double precision, dimension(118,0:5,10) :: bECP, aECP

! Los pseudopotenciales vienen dados por:

!                     _LMAX-1         _l
! \  / = \  /     +  \    \  /       \    |l,m\/l,m|
!  \/     \/LMAX     /_    \/l-LMAX  /_   |   /\   |
!                    l=0             m=-l

!donde Vl= aECP * r^nECP * exp(-bECP r^2)

! Los  xECP estan escritos como: xECP(Z,l,i) Z carga del nucleo, l 
!momento angular de la expansion del ecp, i numero de funcion del
! ecp con Z,l

	integer, dimension (:), ALLOCATABLE :: IzECP
!cargas nucleares sin corregir por la carga del core (Zcore)

	integer, dimension (:,:), ALLOCATABLE :: Lxyz
! Lxyz(i,j) contiene los exponentes de la parte angular de la funcion
! de base i
!|xi> = ci x^lx y^ly z^lz *e^(-a * r^2)
!j=1 lx, j=2 ly, j=3 lz para la funcion i de la base

        double precision, dimension(:,:), Allocatable :: VAAAcuadrada, VAABcuadrada, VBACcuadrada
        double precision, dimension(:), Allocatable :: VAAA, VAAB, VBAC,term1e
! VXXX contiene los coeficientes de Fock del pseudopotencial.
! VAAA integrales de un centro (base y ecp en el mismo atomo)
! VAAB integrales de 2 centros (1 base y ecp en el mismo atomo)
! VBAC integrales de 3 centros (ninguna base en el atomo con ecp)
! term1e contiene una copia de los terminos de 1e- sin la modificaion
! por agregar los terminos de los pseudopotenciales

	double precision, dimension(:,:), Allocatable :: distx, disty, distz
!guarda la distancia en x, y, z entre los atomos i y j  dist(i,j)=xi-xj
!Cuidado, esta en unidades atomicas
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!



!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%%%%%%%%%%%%%%%%    Normalized Basis Coeficients    %%%%%%%%%%%%%%%%!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
	double precision, dimension(:,:), Allocatable :: Cnorm
! Lio guarda en el array c los coeficientes de la base 
! |xi> = ci x^lx y^ly z^lz *e^(-a * r^2)
! En el caso de las funciones d el coeficiente ci esta normalizado 
! para todos los terminos como si fueran xy, xz o yz.
! para los terminos xx, yy y zz hay que dividir por 3^0.5
! en los calculos de lio este factor ya se considera, mientras que en 
! los calculos con pseudo potenciales se modica la base copiandola 
! a Cnorm para hacer a la rutina mas facil de adaptar a otras 
! implementaciones
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%%%%%%%%%%%%%    Parameters For Radial Integration    %%%%%%%%%%%%%%!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!

	integer, parameter, dimension (4,4) :: alpha = (/1,0,1,0,0,-3,0,-10,0,0,15,0,0,0,0,-105/)
	integer, parameter, dimension (5,5) :: betha = (/1,0,1,0,1,0,-1,0,-6,0,0,0,3,0,45,0,0,0,-15,0,0,0,0,0,105/)

!alpha y betha contains coeficients for expantion of modified
! spherical bessel function of the first kind Mk(x) in terms of 
!sinh(x)/x^i (betha(k+1,i) and cosh(x)/x^i (alpha(k,i)) for k between
! 0 and 4, it is enought for energy calculations of functions s to g

        DOUBLE PRECISION, dimension (-12:14) :: Bn1, Bn2, Cn1, Cn2, rho, tau, sigma, sigmaR
        DOUBLE PRECISION, dimension (-12:14) :: Bn,Cn

!Bni(j) contain the value of  int t^(n-2-j)*(exp)-c(t-ai)^2 + exp(-c(t+ai)^2) from 0 to inf
!Cni(j) contain the value of  int t^(n-2-j)*(exp)-c(t-ai)^2 - exp(-c(t+ai)^2) from 0 to inf

! rho(n) = int exp(-cr^2) * sinh(Ka*r)* sinh(Kb*r) r^n  dr  from 0 to inf
! sigma(n) = int exp(-cr^2) * sinh(Ka*r)* cosh(Kb*r) r^n  dr  from 0 to inf
! sigmaR(n) = int exp(-cr^2) * cosh(Ka*r)* sinh(Kb*r) r^n  dr  from 0 to inf
! tau(n) = int exp(-cr^2) * cosh(Ka*r)* cosh(Kb*r) r^n  dr  from 0 to inf

        DOUBLE PRECISION, dimension (0:10,0:4,0:4) :: Qnl1l2
!Qnl1l2(n,l1,l2) = int Ml1(kA*r)* Ml2(kB*r)*r^n * exp(-cr^2) dr from 0 to inf

        DOUBLE PRECISION, dimension (0:10,0:4) :: Qnl
!Qnl(n,l) = int Ml(k*r)*r^n * exp(-cr^2) dr from 0 to inf

!Parameters
	DOUBLE PRECISION, parameter :: pi=3.14159265358979312D0, pi12=1.77245385090552D0 !pi12 = pi^0.5
!factorial
	integer, dimension(0:15) :: fac = (/1,1,2,6,24,120,720,5040,40320,362880,3628800,39916800,479001600,1932053504,1278945280,2004310016/)
!double factorial
        integer*8, parameter, dimension (-1:33) :: doublefactorial=(/1,1,1,2,3,8,15,48,105,384,945,3840,10395,46080,135135,645120,2027025,10321920,34459425,185794560,654729075,3715891200,13749310575,81749606400,316234143225,1961990553600,7905853580625,51011754393600,213458046676875,1428329123020800,6190283353629375,42849873690624000,191898783962510625,1371195958099968000,6332659870762850625/)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%%%%%%%%%%%%%    Parameters for Angular Integrals    %%%%%%%%%%%%%%%!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!

!li contains the coeficiente ul,m,lx,ly,lz fof expansion of 
! Normalized real spherical harmonics Slm in terms of unitary sphere
! polinomials
!S(l,m)= sum ul(m,lx,ly,lz) * (x/r)^lx * (y/r)^ly * (z/r)^lz
! lx+ly+lz=l
!ul=F(f(lx,ly,lz),m)
        DOUBLE PRECISION, dimension (1) :: l0 = (/0.5d0/pi12/)
	Double precision, parameter :: aux1=sqrt(3.d0)/(2d0*pi12)
	DOUBLE PRECISION, dimension (3,-1:1) :: l1=(/0.d0,aux1,0.d0,aux1,0.d0,0.d0,0.d0,0.d0,aux1/)
        Double precision, parameter :: aux2= 0.5d0 * sqrt(15.d0/pi)
        Double precision, parameter :: aux3= 0.25d0 * sqrt(5.d0/pi)
	DOUBLE PRECISION, dimension (6,-2:2) :: l2=(/0.d0,0.d0,0.d0,0.d0,aux2,0.d0,0.d0,aux2,0.d0,0.d0,0.d0,0.d0,2*aux3,0.d0,-aux3,0.d0,0.d0,-aux3,0.d0,0.d0,0.d0,aux2,0.d0,0.d0,0.d0,0.d0,-0.5d0*aux2,0.d0,0.d0,0.5d0*aux2/)
	Double precision, parameter :: aux4=0.25d0 * sqrt(17.5d0/pi)
	Double precision, parameter :: aux5=0.5d0 * sqrt(105.d0/pi)
	Double precision, parameter :: aux6=sqrt(10.5d0/pi)
	Double precision, parameter :: aux7=sqrt(7.0d0/pi)
	DOUBLE PRECISION, dimension (10,-3:3) :: l3=(/0.d0,0.d0,0.d0,-aux4,0.d0,0.d0,0.d0,0.d0,3.d0*aux4,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,aux5,0.d0,0.d0,0.d0,0.d0,0.d0,aux6,0.d0,-0.25d0*aux6,0.d0,0.d0,0.d0,0.d0,-0.25d0*aux6,0.d0,0.5d0*aux7,0.d0,-0.75d0*aux7,0.d0,0.d0,0.d0,0.d0,-0.75d0*aux7,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,aux6,0.d0,-0.25*aux6,0.d0,0.d0,-0.25d0*aux6,0.d0,0.d0,-0.5d0*aux5,0.d0,0.d0,0.d0,0.d0,0.5d0*aux5,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,-3.d0*aux4,0.d0,0.d0,aux4/)
	Double precision, parameter :: aux8=sqrt(35.d0/pi)
        Double precision, parameter :: aux9=sqrt(35.d0/(2*pi))
        Double precision, parameter :: aux10=sqrt(5.d0/pi)
        Double precision, parameter :: aux11=sqrt(5.d0/(2*pi))
        Double precision, parameter :: aux12=1.d0/pi12
        DOUBLE PRECISION, dimension (15,-4:4) :: l4=(/0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,-0.75d0*aux8,0.d0,0.d0,0.d0,0.d0,0.75d0*aux8,0.d0,0.d0,0.d0,0.d0,-0.75d0*aux9,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,2.25d0*aux9,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,4.5d0*aux10,0.d0,-0.75d0*aux10,0.d0,0.d0,0.d0,0.d0,-0.75d0*aux10,0.d0,0.d0,3.d0*aux11,0.d0,-2.25d0*aux11,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,-2.25d0*aux11,0.d0,0.d0,0.d0,0.d0,1.5d0*aux12,0.d0,-4.5d0*aux12,0.d0,0.5625d0*aux12,0.d0,0.d0,0.d0,0.d0,-4.5d0*aux12,0.d0,1.125d0*aux12,0.d0,0.d0,0.5625d0*aux12,0.d0,0.d0,0.d0,0.d0,0.d0,3.d0*aux11,0.d0,-2.25d0*aux11,0.d0,0.d0,0.d0,0.d0,-2.25d0*aux11,0.d0,0.d0,0.d0,0.d0,-2.25d0*aux10,0.d0,0.375d0*aux10,0.d0,0.d0,0.d0,0.d0,2.25*aux10,0.d0,0.d0,0.d0,0.d0,-0.375*aux10,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,-2.25d0*aux9,0.d0,0.d0,0.d0,0.d0,0.75d0*aux9,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.1875d0*aux8,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,-1.125d0*aux8,0.d0,0.d0,0.1875d0*aux8/)


        double precision, dimension (0:9,0:9,0:9) :: angularint
! contiene los valores tabulados de  la rutina del mismo nombre (comentada mas abajo)

	Contains
	
	subroutine defineparams()
!define parametros de variables del modulo
	implicit none
!parametros auxiliares para integrales angulares
!angular int
	angularint=0.D0
        angularint(0,0,0)=12.5663706143592D0
        angularint(0,0,2)=4.18879020478639D0
        angularint(0,0,4)=2.51327412287183D0
        angularint(0,0,6)=1.79519580205131D0
        angularint(0,0,8)=1.39626340159546D0
        angularint(0,2,0)=4.18879020478639D0
        angularint(0,2,2)=0.837758040957278D0
        angularint(0,2,4)=0.359039160410262D0
        angularint(0,2,6)=0.199466200227923D0
        angularint(0,2,8)=0.126933036508679D0
        angularint(0,4,0)=2.51327412287183D0
        angularint(0,4,2)=0.359039160410262D0
        angularint(0,4,4)=0.119679720136754D0
        angularint(0,4,6)=5.439987278943365D-2
        angularint(0,4,8)=2.929223919431043D-2
        angularint(0,6,0)=1.79519580205131D0
        angularint(0,6,2)=0.199466200227923D0
        angularint(0,6,4)=5.439987278943364D-2
        angularint(0,6,6)=2.092302799593602D-2
        angularint(0,6,8)=9.764079731436807D-3
        angularint(0,8,0)=1.39626340159546D0
        angularint(0,8,2)=0.126933036508679D0
        angularint(0,8,4)=2.929223919431043D-2
        angularint(0,8,6)=9.764079731436809D-3
        angularint(0,8,8)=4.020503418826921D-3
        angularint(2,0,0)=4.18879020478639D0
        angularint(2,0,2)=0.837758040957278D0
        angularint(2,0,4)=0.359039160410262D0
        angularint(2,0,6)=0.199466200227923D0
        angularint(2,0,8)=0.126933036508679D0
        angularint(2,2,0)=0.837758040957278D0
        angularint(2,2,2)=0.119679720136754D0
        angularint(2,2,4)=3.989324004558467D-2
        angularint(2,2,6)=1.813329092981121D-2
        angularint(2,2,8)=9.764079731436809D-3
        angularint(2,4,0)=0.359039160410262D0
        angularint(2,4,2)=3.989324004558467D-2
        angularint(2,4,4)=1.087997455788673D-2
        angularint(2,4,6)=4.184605599187203D-3
        angularint(2,4,8)=1.952815946287362D-3
        angularint(2,6,0)=0.199466200227923D0
        angularint(2,6,2)=1.813329092981121D-2
        angularint(2,6,4)=4.184605599187203D-3
        angularint(2,6,6)=1.394868533062401D-3
        angularint(2,6,8)=5.743576312609886D-4
        angularint(2,8,0)=0.126933036508679D0
        angularint(2,8,2)=9.764079731436809D-3
        angularint(2,8,4)=1.952815946287362D-3
        angularint(2,8,6)=5.743576312609887D-4
        angularint(2,8,8)=2.116054430961537D-4
        angularint(4,0,0)=2.51327412287183D0
        angularint(4,0,2)=0.359039160410262D0
        angularint(4,0,4)=0.119679720136754D0
        angularint(4,0,6)=5.439987278943365D-2
        angularint(4,0,8)=2.929223919431043D-2
        angularint(4,2,0)=0.359039160410262D0
        angularint(4,2,2)=3.989324004558467D-2
        angularint(4,2,4)=1.087997455788673D-2
        angularint(4,2,6)=4.184605599187203D-3
        angularint(4,2,8)=1.952815946287362D-3
        angularint(4,4,0)=0.119679720136754D0
        angularint(4,4,2)=1.087997455788673D-2
        angularint(4,4,4)=2.510763359512322D-3
        angularint(4,4,6)=8.369211198374407D-4
        angularint(4,4,8)=3.446145787565932D-4
        angularint(4,6,0)=5.439987278943365D-2
        angularint(4,6,2)=4.184605599187203D-3
        angularint(4,6,4)=8.369211198374408D-4
        angularint(4,6,6)=2.461532705404238D-4
        angularint(4,6,8)=9.068804704120875D-5
        angularint(4,8,0)=2.929223919431043D-2
        angularint(4,8,2)=1.952815946287362D-3
        angularint(4,8,4)=3.446145787565932D-4
        angularint(4,8,6)=9.068804704120875D-5
        angularint(4,8,8)=3.022934901373625D-5
        angularint(6,0,0)=1.79519580205131D0
        angularint(6,0,2)=0.199466200227923D0
        angularint(6,0,4)=5.439987278943364D-2
        angularint(6,0,6)=2.092302799593602D-2
        angularint(6,0,8)=9.764079731436807D-3
        angularint(6,2,0)=0.199466200227923D0
        angularint(6,2,2)=1.813329092981121D-2
        angularint(6,2,4)=4.184605599187203D-3
        angularint(6,2,6)=1.394868533062401D-3
        angularint(6,2,8)=5.743576312609886D-4
        angularint(6,4,0)=5.439987278943364D-2
        angularint(6,4,2)=4.184605599187203D-3
        angularint(6,4,4)=8.369211198374406D-4
        angularint(6,4,6)=2.461532705404237D-4
        angularint(6,4,8)=9.068804704120873D-5
        angularint(6,6,0)=2.092302799593602D-2
        angularint(6,6,2)=1.394868533062401D-3
        angularint(6,6,4)=2.461532705404238D-4
        angularint(6,6,6)=6.477717645800625D-5
        angularint(6,6,8)=2.159239215266875D-5
        angularint(6,8,0)=9.764079731436807D-3
        angularint(6,8,2)=5.743576312609886D-4
        angularint(6,8,4)=9.068804704120873D-5
        angularint(6,8,6)=2.159239215266875D-5
        angularint(6,8,8)=6.571597611681793D-6
        angularint(8,0,0)=1.39626340159546D0
        angularint(8,0,2)=0.126933036508679D0
        angularint(8,0,4)=2.929223919431043D-2
        angularint(8,0,6)=9.764079731436809D-3
        angularint(8,0,8)=4.020503418826921D-3
        angularint(8,2,0)=0.126933036508679D0
        angularint(8,2,2)=9.764079731436809D-3
        angularint(8,2,4)=1.952815946287362D-3
        angularint(8,2,6)=5.743576312609887D-4
        angularint(8,2,8)=2.116054430961537D-4
        angularint(8,4,0)=2.929223919431043D-2
        angularint(8,4,2)=1.952815946287362D-3
        angularint(8,4,4)=3.446145787565932D-4
        angularint(8,4,6)=9.068804704120875D-5
        angularint(8,4,8)=3.022934901373625D-5
        angularint(8,6,0)=9.764079731436809D-3
        angularint(8,6,2)=5.743576312609887D-4
        angularint(8,6,4)=9.068804704120875D-5
        angularint(8,6,6)=2.159239215266875D-5
        angularint(8,6,8)=6.571597611681793D-6
        angularint(8,8,0)=4.020503418826921D-3
        angularint(8,8,2)=2.116054430961537D-4
        angularint(8,8,4)=3.022934901373625D-5
        angularint(8,8,6)=6.571597611681793D-6
        angularint(8,8,8)=1.840047331270902D-6
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
	end subroutine defineparams


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%%%%%%%%%%%%%%%%%%%    Other Functions/Subs    %%%%%%%%%%%%%%%%%%%%%!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!

        subroutine asignacion(Z,simb)
!Toma el numero atomico (Z)  y devuelve el simbolo del elemento en mayusculas (simb)
        implicit none
        integer,intent(in) :: Z
        character, intent(out) :: simb*3
        character (len=3), dimension(118) :: vec
        vec=(/'H', 'HE', 'LI', 'BE', 'B', 'C', 'N', 'O', 'F', 'NE', 'NA', 'MG', 'AL', 'SI', 'P', 'S', 'CL', 'AR', 'K', 'CA', 'SC', 'TI', 'V', 'CR', 'MN', 'FE', 'CO', 'NI', 'CU', 'ZN', 'GA', 'GE', 'AS', 'SE', 'BR', 'KR', 'RB', 'SR', 'Y', 'ZR', 'NB', 'MO', 'TC', 'RU', 'RH', 'PD', 'AG', 'CD', 'IN', 'SN', 'SB', 'TE', 'I', 'XE', 'CS', 'BA', 'LA', 'CE', 'PR', 'ND', 'PM', 'SM', 'EU', 'GD', 'TB', 'DY', 'HO', 'ER', 'TM', 'YB', 'LU', 'HF', 'TA', 'W', 'RE', 'OS', 'IR', 'PT', 'AU', 'HG', 'TL', 'PB', 'BI', 'PO', 'AT', 'RN', 'FR', 'RA', 'AC', 'TH', 'PA', 'U', 'NP', 'PU', 'AM', 'CM', 'BK', 'CF', 'ES', 'FM', 'MD', 'NO', 'LR', 'RF', 'DB', 'SG', 'BH', 'HS', 'MT', 'DS', 'UUU', 'UUB', 'UUT', 'UUQ', 'UUP', 'UUH', 'UUS', 'UUO'/)
        simb=vec(Z)
        return
        end subroutine asignacion


        DOUBLE PRECISION FUNCTION NEXTCOEF(sgn,n,cados,expo,c0coef,coefn1, coefn2)
!calcula el coef Bn y Cn, con n<2
!si sgn=1 calcula Bn; si sgn=-1 calcula Cn
        integer, intent(in) :: sgn,n
        DOUBLE PRECISION, intent(in) :: cados,expo,c0coef,coefn1, coefn2
	if ( -n-1 .lt. 0) stop " se pide fac(n), n<0 en NEXTCOEF"
        NEXTCOEF=(1+sgn*(-1)**(-n-1))*cados**(-n-1)*expo/fac(-n-1)-2*c0coef*coefn2 +cados*coefn1
	NEXTCOEF=NEXTCOEF/(-n-1)
        RETURN
        END Function NEXTCOEF



!        recursive double precision function facrecursive(N)
! rutina reemplazada por el array fac
!        implicit none
!        integer, intent(in) :: N
!        if ( n .gt. 170) then
!        stop "n greather than 170 in facrecursive function"
!        end if
!        if (n.eq.-1 .or. n.eq.0 .or. n.eq.1) then
!        facrecursive=1.d0
!        else
!        facrecursive=N*facrecursive(N-1)
!        end if
!        return
!        end function facrecursive

!        recursive integer*8 function doublefacrecursive(N)
! rutina reemplazada por el array doublefactorial
!la rutina funciona hasta N=19
!        implicit none
!        integer, intent(in) :: N
!        if (n.eq.-1 .or. n.eq.0 .or. n.eq.1) then
!        doublefacrecursive=1.d0
!	else if (n .gt. 33) then
!		write(*,*) "Error in doble factorial routine, n greater than 33"
!		write(*,*) "insuficient presicion"
!        else
!        doublefacrecursive=N*doublefacrecursive(N-2)
!        end if
!        return
!        end function doublefacrecursive


!        integer*8 function doublefac(N)
!esta rutina fue reemplazada por un array
!        implicit none
!        integer, intent(in) :: N
!        integer :: i
!        doublefac=1
!        if (n .gt. 33) then
!                stop "error in doublefac function N greater than 33"
!        end if
!        doublefac=1.d0
!        if (n.gt.1) then
!                do i=N,2,-2
!                        doublefac=doublefac*i
!                end do
!        end if
!        return
!        end function doublefac

!	DOUBLE PRECISION function angularint(i,j,k)
!	implicit none
!	integer, intent(in) :: i,j,k
!	DOUBLE PRECISION, parameter :: pi=3.14159265358979312D0
!esta rutina calcula la integral (x/r)^i (y/r)^j (z/r)^k, luego sera reemplazada por un array
!	if (mod(i,2) .eq. 0 .and. mod(j,2) .eq. 0 .and.mod(k,2) .eq. 0 ) then
!	   angularint=4.d0*pi*doublefactorial(i-1)*doublefactorial(j-1)*doublefactorial(k-1)/doublefactorial(i+j+k+1)
!        else
!           angularint=0.d0
!        end if
!	return
!	end function angularint

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%%%%%%%%%%%%%%%%%%%    Other Functions/Subs    %%%%%%%%%%%%%%%%%%%%%!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!

!      double precision  function Ran()
!      implicit none
!      integer :: count1, count_rate, count_max
!      CALL SYSTEM_CLOCK(count1, count_rate, count_max)
!      Ran=2*RANDOM(count1)-1
!      return
!      end function Ran

!      Function RANDOM(Seed)
!      Implicit None
!      Real :: Random
!      Integer :: Seed
!      Integer :: OldSeed = 0
!      Integer, Parameter :: C1 = 19423
!      Integer, Parameter :: C2 = 811
!      Save OldSeed

!      If (OldSeed .EQ. 0) OldSeed = Seed
!      OldSeed = Mod(C1 * OldSeed, C2)
!      RANDOM = 1.0 * OldSeed / C2
!      End Function RANDOM



!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%%%%%%%%%%%%%%%%%%%    Not Made By Nicolas    %%%%%%%%%%%%%%%%%%%%%%!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!


	FUNCTION DAW(XX) RESULT(fn_val)
! Code converted using TO_F90 by Alan Miller
! Date: 2003-01-14  Time: 15:25:00
 
!----------------------------------------------------------------------

! This function program evaluates Dawson's integral,

!                       2  / x   2
!                     -x   |    t
!             F(x) = e     |   e    dt
!                          |
!                          / 0

!   for a real argument x.

!   The calling sequence for this function is

!                   Y=DAW(X)

!   The main computation uses rational Chebyshev approximations published
!   in Math. Comp. 24, 171-178 (1970) by Cody, Paciorek and Thacher.
!   This transportable program is patterned after the machine-dependent
!   FUNPACK program DDAW(X), but cannot match that version for efficiency or
!   accuracy.  This version uses rational approximations that are
!   theoretically accurate to about 19 significant decimal digits.
!   The accuracy achieved depends on the arithmetic system, the compiler,
!   the intrinsic functions, and proper selection of the machine-dependent
!   constants.

!*******************************************************************

! Explanation of machine-dependent constants.  Let

!   XINF   = largest positive machine number
!   XMIN   = the smallest positive machine number.
!   EPS    = smallest positive number such that 1+eps > 1.
!            Approximately  beta**(-p), where beta is the machine radix
!            and p is the number of significant base-beta digits in a
!            floating-point number.

! Then the following machine-dependent constants must be declared
!   in DATA statements.  IEEE values are provided as a default.

!   XMAX   = absolute argument beyond which DAW(X) underflows.
!            XMAX = min(0.5/xmin, xinf).
!   XSMALL = absolute argument below DAW(X)  may be represented
!            by X.  We recommend XSMALL = sqrt(eps).
!   XLARGE = argument beyond which DAW(X) may be represented by
!            1/(2x).  We recommend XLARGE = 1/sqrt(eps).

!     Approximate values for some important machines are

!                        beta  p     eps     xmin       xinf

!  CDC 7600      (S.P.)    2  48  7.11E-15  3.14E-294  1.26E+322
!  CRAY-1        (S.P.)    2  48  7.11E-15  4.58E-2467 5.45E+2465
!  IEEE (IBM/XT,
!    SUN, etc.)  (S.P.)    2  24  1.19E-07  1.18E-38   3.40E+38
!  IEEE (IBM/XT,
!    SUN, etc.)  (D.P.)    2  53  1.11D-16  2.23E-308  1.79D+308
!  IBM 3033      (D.P.)   16  14  1.11D-16  5.40D-79   7.23D+75
!  VAX 11/780    (S.P.)    2  24  5.96E-08  2.94E-39   1.70E+38
!                (D.P.)    2  56  1.39D-17  2.94D-39   1.70D+38
!   (G Format)   (D.P.)    2  53  1.11D-16  5.57D-309  8.98D+307

!                         XSMALL     XLARGE     XMAX

!  CDC 7600      (S.P.)  5.96E-08   1.68E+07  1.59E+293
!  CRAY-1        (S.P.)  5.96E-08   1.68E+07  5.65E+2465
!  IEEE (IBM/XT,
!    SUN, etc.)  (S.P.)  2.44E-04   4.10E+03  4.25E+37
!  IEEE (IBM/XT,
!    SUN, etc.)  (D.P.)  1.05E-08   9.49E+07  2.24E+307
!  IBM 3033      (D.P.)  3.73D-09   2.68E+08  7.23E+75
!  VAX 11/780    (S.P.)  2.44E-04   4.10E+03  1.70E+38
!                (D.P.)  3.73E-09   2.68E+08  1.70E+38
!   (G Format)   (D.P.)  1.05E-08   9.49E+07  8.98E+307

!*******************************************************************

! Error Returns

!  The program returns 0.0 for |X| > XMAX.

! Intrinsic functions required are:

!     ABS


!  Author: W. J. Cody
!          Mathematics and Computer Science Division
!          Argonne National Laboratory
!          Argonne, IL 60439

!  Latest modification: March 9, 1992

!----------------------------------------------------------------------
	IMPLICIT NONE
	INTEGER, PARAMETER  :: dp = SELECTED_REAL_KIND(12, 60)

	REAL (dp), INTENT(IN)  :: xx
	REAL (dp)              :: fn_val

! Local variables

	INTEGER    :: i
	REAL (dp)  :: frac, sump, sumq, w2, x, y
!----------------------------------------------------------------------
!  Mathematical constants.
!----------------------------------------------------------------------
	REAL (dp), PARAMETER  :: zero = 0.0_dp, half = 0.5_dp, one = 1.0_dp,  &
                         six25 = 6.25_dp, one225 = 12.25_dp, two5 = 25.0_dp
!----------------------------------------------------------------------
!  Machine-dependent constants
!----------------------------------------------------------------------
	REAL (dp), PARAMETER  :: XSMALL = 1.05D-08, XLARGE = 9.49D+07,   &
                         XMAX = 2.24D+307
!----------------------------------------------------------------------
!  Coefficients for R(9,9) approximation for  |x| < 2.5
!----------------------------------------------------------------------
	REAL (dp), PARAMETER  :: P1(10) = (/  &
        -2.69020398788704782410D-12, 4.18572065374337710778D-10,  &
        -1.34848304455939419963D-08, 9.28264872583444852976D-07,  &
        -1.23877783329049120592D-05, 4.07205792429155826266D-04,  &
        -2.84388121441008500446D-03, 4.70139022887204722217D-02,  &
        -1.38868086253931995101D-01, 1.00000000000000000004D+00 /)
	REAL (dp), PARAMETER  :: Q1(10) = (/  &
         1.71257170854690554214D-10, 1.19266846372297253797D-08,  &
         4.32287827678631772231D-07, 1.03867633767414421898D-05,  &
         1.78910965284246249340D-04, 2.26061077235076703171D-03,  &
         2.07422774641447644725D-02, 1.32212955897210128811D-01,  &
         5.27798580412734677256D-01, 1.00000000000000000000D+00 /)
!----------------------------------------------------------------------
!  Coefficients for R(9,9) approximation in J-fraction form
!     for  x in [2.5, 3.5)
!----------------------------------------------------------------------
	REAL (dp), PARAMETER  :: P2(10) = (/  &
        -1.70953804700855494930D+00,-3.79258977271042880786D+01,  &
         2.61935631268825992835D+01, 1.25808703738951251885D+01,  &
        -2.27571829525075891337D+01, 4.56604250725163310122D+00,  &
        -7.33080089896402870750D+00, 4.65842087940015295573D+01,  &
        -1.73717177843672791149D+01, 5.00260183622027967838D-01 /)
	REAL (dp), PARAMETER  :: Q2(9) = (/  &
         1.82180093313514478378D+00, 1.10067081034515532891D+03,  &
        -7.08465686676573000364D+00, 4.53642111102577727153D+02,  &
         4.06209742218935689922D+01, 3.02890110610122663923D+02,  &
         1.70641269745236227356D+02, 9.51190923960381458747D+02,  &
         2.06522691539642105009D-01 /)
!----------------------------------------------------------------------
!  Coefficients for R(9,9) approximation in J-fraction form
!     for  x in [3.5, 5.0]
!----------------------------------------------------------------------
	REAL (dp), PARAMETER  :: P3(10) = (/  &
        -4.55169503255094815112D+00,-1.86647123338493852582D+01,  &
        -7.36315669126830526754D+00,-6.68407240337696756838D+01,  &
         4.84507265081491452130D+01, 2.69790586735467649969D+01,  &
        -3.35044149820592449072D+01, 7.50964459838919612289D+00,  &
        -1.48432341823343965307D+00, 4.99999810924858824981D-01 /)
	REAL (dp), PARAMETER  :: Q3(9) = (/  &
         4.47820908025971749852D+01, 9.98607198039452081913D+01,  &
         1.40238373126149385228D+01, 3.48817758822286353588D+03,  &
        -9.18871385293215873406D+00, 1.24018500009917163023D+03,  &
        -6.88024952504512254535D+01,-2.31251575385145143070D+00,  &
         2.50041492369922381761D-01 /)
!----------------------------------------------------------------------
!  Coefficients for R(9,9) approximation in J-fraction form
!     for  |x| > 5.0
!----------------------------------------------------------------------
	REAL (dp), PARAMETER  :: P4(10) = (/  &
        -8.11753647558432685797D+00,-3.84043882477454453430D+01,  &
        -2.23787669028751886675D+01,-2.88301992467056105854D+01,  &
        -5.99085540418222002197D+00,-1.13867365736066102577D+01,  &
        -6.52828727526980741590D+00,-4.50002293000355585708D+00,  &
        -2.50000000088955834952D+00, 5.00000000000000488400D-01 /)
	REAL (dp), PARAMETER  :: Q4(9) = (/  &
         2.69382300417238816428D+02, 5.04198958742465752861D+01,  &
         6.11539671480115846173D+01, 2.08210246935564547889D+02,  &
         1.97325365692316183531D+01,-1.22097010558934838708D+01,  &
        -6.99732735041547247161D+00,-2.49999970104184464568D+00,  &
         7.49999999999027092188D-01 /)
!----------------------------------------------------------------------

	x = xx
	IF (ABS(x) > xlarge) THEN
	  IF (ABS(x) <= xmax) THEN
	    fn_val = half / x
	  ELSE
	    fn_val = zero
	  END IF
	ELSE IF (ABS(x) < xsmall) THEN
	  fn_val = x
	ELSE
	  y = x * x
	  IF (y < six25) THEN
!----------------------------------------------------------------------
!  ABS(X) .LT. 2.5
!----------------------------------------------------------------------
	    sump = p1(1)
	    sumq = q1(1)
	    DO  i = 2, 10
	      sump = sump * y + p1(i)
	      sumq = sumq * y + q1(i)
	    END DO
	    fn_val = x * sump / sumq
	  ELSE IF (y < one225) THEN
!----------------------------------------------------------------------
!  2.5 .LE. ABS(X) .LT. 3.5
!----------------------------------------------------------------------
	    frac = zero
	    DO  i = 1, 9
	      frac = q2(i) / (p2(i)+y+frac)
	    END DO
	    fn_val = (p2(10)+frac) / x
	  ELSE IF (y < two5) THEN
!----------------------------------------------------------------------
!  3.5 .LE. ABS(X) .LT. 5.0
!---------------------------------------------------------------------
	    frac = zero
	    DO  i = 1, 9
	      frac = q3(i) / (p3(i)+y+frac)
	    END DO
	    fn_val = (p3(10)+frac) / x
	  ELSE
!----------------------------------------------------------------------
!  5.0 .LE. ABS(X) .LE. XLARGE
!------------------------------------------------------------------
	    w2 = one / x / x
	    frac = zero
	    DO  i = 1, 9
	      frac = q4(i) / (p4(i)+y+frac)
	    END DO
	    frac = p4(10) + frac
	    fn_val = (half + half*w2*frac) / x
	  END IF
	END IF
	Return
!---------- Last line of DAW ----------
	END FUNCTION DAW

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!


	end module ECP_mod
