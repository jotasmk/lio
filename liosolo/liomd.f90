!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
  program liomd
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
  use garcha_mod
  use ehrensubs
  use basis_data
  use liosubs
#ifdef CUBLAS
  use cublasmath
#endif
  implicit none
  real*8  :: escf,dipxyz(3)
  integer :: ii,kk
  integer :: narg,nqnuc,ios
  integer :: ng2,ng3,ngdDyn,ngdnu,ngDyn,ngnu

  character(len=20)::argument,inpfile,inpbasis,inpcoords
  integer::charge
  logical::filexist,writeforces
  REAL*8, dimension (:,:), ALLOCATABLE   :: dxyzqm,dxyzcl

  namelist /lio/ natom,nsol,charge,OPEN,NMAX,Nunp,VCINP,frestartin, &
  GOLD,told,rmax,rmaxs,predcoef, &
  idip,writexyz,intsoldouble,DIIS,ndiis,dgtrig, &
  Iexch,integ,dens,igrid,igrid2,timedep, tdstep, ntdstep, &
  propagator,NBCH, &
  field,a0,epsilon,exter,Fx,Fy,Fz, tdrestart, writedens, &
  writeforces,basis_set,fitting_set,int_basis, &
  cubegen_only,cube_res, &
  cube_dens,cube_dens_file, &
  cube_orb,cube_sel,cube_orb_file,cube_elec,cube_elec_file

  integer             :: nn
  integer             :: md_steps
  real*8              :: time_step,Kenergy,Uenergy
  real*8,allocatable  :: oldpos(:,:),nowpos(:,:),newpos(:,:)
  real*8,allocatable  :: nowvel(:,:),nowfza(:,:)
!  real*8,allocatable  :: atom_mass(:)

 !defaults
  basis='input'  ! name of the base file
  basis_set='DZVP'
  fitting_set='DZVP Coulomb Fitting'
  int_basis=.false.
  cubegen_only=.false.
  cube_res=40
  cube_dens=.false.
  cube_dens_file='dens.cube'
  cube_orb=.false.
  cube_sel=0
  cube_orb_file="orb.cube"
  cube_elec=.false.
  cube_elec_file="field.cube"
  restart_freq=1
  energy_freq=1
  output='output'
  fcoord='qm.xyz'
  fmulliken='mulliken'
  frestart='restart.out'
  frestartin='restart.in'
  verbose=.false.
  OPEN=.false.
  NMAX= 100
  NUNP= 0
  VCINP= .false.
  GOLD= 10.
  told=1.0D-6
  rmax=16
  rmaxs=5
  predcoef=.false.
  idip=1
  writexyz=.true.
  intsoldouble=.true.
  DIIS=.true.
  ndiis=30
  dgtrig=100.
  Iexch=9
  integ=.true.
  DENS = .true.
  IGRID = 2
  IGRID2 = 2
  timedep = 0
  tdstep = 2.D-3
  field=.false.
  a0=1000.0
  epsilon=1.D0
  Fx=0.05
  Fy=0.05
  Fz=0.05
  NBCH=10
  propagator=1
  tdrestart=.false.
  writedens=.true.
  writeforces=.false.
  narg=command_argument_count()

  do ii=1, narg
    call get_command_argument(ii,argument)
    select case(adjustl(argument))
      case("-i")
    call get_command_argument(ii+1,inpfile)
      case("-b")
    call get_command_argument(ii+1,basis)
      case("-bs")
    call get_command_argument(ii+1,basis_set)
      case("-fs")
    call get_command_argument(ii+1,fitting_set)
      case("-ib")
    int_basis=.true.
    !call get_command_argument(i+1,int_basis)
      case("-c")
    call get_command_argument(ii+1,inpcoords)
      case("-v")
    verbose=.true.
      case default
    end select
  enddo

  call g2g_timer_sum_start("Total")

  inquire(file=inpfile,exist=filexist)

  if(filexist) then
    open(unit=100,file=inpfile,iostat=ios)
  else
    write(*,*) 'input file ',adjustl(inpfile),' not found'
    stop
  endif
  read(100,nml=lio,iostat=ios)

  if(ios.gt.0) stop 'input error in lio namelist'

  inquire(file=inpcoords,exist=filexist)
  if(filexist) then
    open(unit=101,file=inpcoords,iostat=ios)
  else
    write(*,*) 'input file ',adjustl(inpcoords),' not found'
    stop
  endif
  write(*,nml=lio)

  ntatom=natom+nsol
  ngnu=natom*ng0
  ngdnu=natom*ngd0
  ngDyn=ngnu
  ngdDyn=ngdnu
  ng3=4*ngDyn
  ng2=5*ngDyn*(ngDyn+1)/2+3*ngdDyn*(ngdDyn+1)/2+&
        ngDyn+ngDyn*norbit+Ngrid

  allocate(X(ngDyn,ng3),XX(ngdDyn,ngdDyn))
  allocate(RMM(ng2))

  allocate (c(ngnu,nl),a(ngnu,nl),Nuc(ngnu),ncont(ngnu) &
   ,cx(ngdnu,nl),ax(ngdnu,nl),Nucx(ngdnu),ncontx(ngdnu) &
  ,cd(ngdnu,nl),ad(ngdnu,nl),Nucd(ngdnu),ncontd(ngdnu) &
  ,indexii(ngnu),indexiid(ngdnu))

  allocate (r(ntatom,3),v(ntatom,3),rqm(natom,3),Em(ntatom) &
  ,Rm(ntatom),pc(ntatom),Iz(natom),af(natom*ngd0), &
   B(natom*ngd0,3))
  allocate (nnat(100))
  allocate(d(natom,natom))

  do ii=1,natom
    read(101,*) iz(ii),r(ii,1:3)
    rqm(ii,1:3)=r(ii,1:3)
  enddo
  do ii=natom+1,ntatom
    read(101,*) pc(ii),r(ii,1:3)   ! o es pc(i-natom)???
   enddo
   r=r/0.529177D0
   rqm=rqm/0.529177D0

   call g2g_init()   !initialize g2g

    nqnuc=0
   do ii=1,natom
     nqnuc=nqnuc+Iz(ii)
   enddo

  nco=((nqnuc - charge)-Nunp)/2

  call drive(ng2,ngDyn,ngdDyn)
!  call lio_init()   !initialize lio
  call liomain()

! ALLOCATION AND INITIALIZATION
!--------------------------------------------------------------------!
  if (.not.allocated(Smat))    allocate(Smat(M,M))
  if (.not.allocated(RealRho)) allocate(RealRho(M,M))
  if (allocated(dxyzqm)) deallocate (dxyzqm)
  allocate(dxyzqm(3,natom))
  allocate(oldpos(natom,3),nowpos(natom,3),newpos(natom,3))
  allocate(nowvel(natom,3),nowfza(natom,3))

  if (allocated(atom_mass)) deallocate(atom_mass)
  allocate(atom_mass(natom))
  call set_masses(natom,Iz,atom_mass)

! SETUP
!--------------------------------------------------------------------!
  do kk=1,3
  do nn=1,natom
    oldpos(nn,kk)=r(nn,kk)
    nowpos(nn,kk)=r(nn,kk)
    newpos(nn,kk)=r(nn,kk)
    nowvel(nn,kk)=0.0d0
    nowfza(nn,kk)=0.0d0
  enddo
  enddo

! RUN DYNAMIC
!--------------------------------------------------------------------!
  md_steps=5000
  time_step=0.0001               ! time in ps
  time_step=time_step*41341.0d0  ! time from ps to au
  open(unit=501,file='liomd-trays.xyz')
  open(unit=502,file='liomd-energy.dat')
  call write_energy(-1.0d0,escf,Kenergy,escf+Kenergy,502)

  do nn=1,md_steps+1
    call SCF(escf,dipxyz)
    Uenergy=escf
    call dft_get_qm_forces(dxyzqm)
    do kk=1,3
    do ii=1,natom
      nowfza(ii,kk)=-dxyzqm(kk,ii)
    enddo
    enddo
    call nuclear_verlet(natom,time_step,atom_mass,nowfza, &
                        oldpos,nowpos,newpos,nowvel,Kenergy)
    call write_geom(natom,Iz,nowpos,501)
    call write_energy(nn*tdstep,escf,Kenergy,escf+Kenergy,502)

    oldpos=nowpos
    nowpos=newpos
    do kk=1,3
    do ii=1,natom
      r(ii,kk)   = nowpos(ii,kk)
      rqm(ii,kk) = nowpos(ii,kk)
    enddo
    enddo
  enddo
!
!
!
! Write Energy
!--------------------------------------------------------------------!
  write(*,*) 'SCF ENRGY=',escf
!
!
!
! Write Forces
!--------------------------------------------------------------------!
  if(writeforces) then
    open(unit=123,file='fuerzas')
    if (allocated(dxyzqm)) deallocate(dxyzqm)
    allocate (dxyzqm(3,natom))
    dxyzqm=0.0

    if(nsol.gt.0) then
      allocate (dxyzcl(3,natom+nsol))
      dxyzcl=0.
    endif

    call dft_get_qm_forces(dxyzqm)
    if (nsol.gt.0) then
      call dft_get_mm_forces(dxyzcl,dxyzqm)
    endif

    do kk=1,natom
      write(123,100) kk,dxyzqm(kk,1),dxyzqm(kk,2),dxyzqm(kk,3)
    enddo

    if(nsol.gt.0) then
      do kk=natom,natom+nsol
        write(123,100) kk,dxyzcl(kk,1),dxyzcl(kk,2),dxyzcl(kk,3)
      enddo
    endif
    deallocate (dxyzqm)
    if(nsol.gt.0) deallocate(dxyzcl)
  endif

  call lio_finalize()
100    format (I5,2x,f10.6,2x,f10.6,2x,f10.6)
  stop;end program
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
