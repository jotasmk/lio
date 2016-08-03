      program liosolo
c MAIN SUBROUTINE ----------------------------------------------------
C DFT calculation with gaussian basis sets
c---------------------------------------------------------------------
      use garcha_mod
      use ECP_mod, only : ecpmode, ecptypes, tipeECP, ZlistECP
     > ,cutECP,local_nonlocal, ecp_debug,ecp_full_range_int
     > ,verbose_ECP,Cnorm,FOCK_ECP_read, FOCK_ECP_write,Fulltimer_ECP
     > ,cut2_0,cut3_0
#ifdef CUBLAS
      use cublasmath
#endif
      implicit real*8 (a-h,o-z)

      character(len=20)::argument,inpfile,inpbasis,inpcoords
      integer::charge
      logical::filexist,writeforces
      REAL*8, dimension (:,:), ALLOCATABLE   :: dxyzqm,dxyzcl
      namelist /lio/ natom,nsol,charge,OPEN,NMAX,Nunp,VCINP,frestartin,
     > GOLD,told,Etold,rmax,rmaxs,predcoef,
     > idip,writexyz,intsoldouble,DIIS,ndiis,dgtrig,
     > Iexch,integ,dens,igrid,igrid2,timedep, tdstep, ntdstep,
     > propagator,NBCH,
     > field,a0,epsilon,exter,Fx,Fy,Fz, tdrestart, writedens,
     > writeforces,basis_set,fitting_set,int_basis,

!%% Effective Core Potential Variables %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     > ecpmode,ecptypes,tipeECP,ZlistECP,
     > cutECP,local_nonlocal, ecp_debug,ecp_full_range_int,verbose_ECP,
     > verbose,FOCK_ECP_read, FOCK_ECP_write,
     > Fulltimer_ECP,cut2_0,cut3_0,
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!%% Hybrid damping-diis Variables %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     > hybrid_converg,good_cut,
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!%% Output format Variables %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     > style, allnml,
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

     > cubegen_only,cube_res,
     > cube_dens,cube_dens_file,
     > cube_orb,cube_sel,cube_orb_file,cube_elec,cube_elec_file



      integer :: ifind, ierr


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
      Etold=1.0D-4
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

!%% Hybrid damping-diis Variables %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      hybrid_converg=.false.
      good_cut=1D-5
      omit_bas=.false.

!%% Output format Variables %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      style=.true.
      allnml=.true.

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%%%%%%%%%%%%%    Effective Core Potential Variables    %%%%%%%%%%%%%!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
      ecpmode=.false.
      tipeECP='NOT-DEFINED'
      ZlistECP=0
      ecptypes=0
      cutECP=.true.
      local_nonlocal=0
      ecp_debug=.false.
      ecp_full_range_int=.false.
      verbose_ECP=0
      FOCK_ECP_read=.false.
      FOCK_ECP_write=.false.
      Fulltimer_ECP=.false.
      cut2_0=15.d0
      cut3_0=12.d0
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!

	call LIO_LOGO2()

      do i=1, narg
        call get_command_argument(i,argument)
        select case(adjustl(argument))
          case("-i")
            call get_command_argument(i+1,inpfile)
          case("-b")
            call get_command_argument(i+1,basis)
          case("-bs")
            omit_bas=.true.
            call get_command_argument(i+1,basis_set)
          case("-fs")
            call get_command_argument(i+1,fitting_set)
          case("-ib")
            int_basis=.true.
            !call get_command_argument(i+1,int_basis)
          case("-c")
            call get_command_argument(i+1,inpcoords)
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
      read(100,nml=lio,iostat=ierr)
      if(ierr.gt.0) stop 'input error in lio namelist'



      inquire(file=inpcoords,exist=filexist)
      if(filexist) then
        open(unit=101,file=inpcoords,iostat=ios)
      else

        write(*,*) 'input file ',adjustl(inpcoords),' not found'
        stop
      endif

c%%%%%%%%%%%%%%%%%%   Namelist write   %%%%%%%%%%%%%%%%%%
c        write(*,*) natom,nsol
      if (allnml) write(*,nml=lio)
      if (style) call NEW_WRITE_NML(charge)
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      ntatom=natom+nsol
      ngnu=natom*ng0
      ngdnu=natom*ngd0
      ngDyn=ngnu
      ngdDyn=ngdnu
      ng3=4*ngDyn

      ng2=5*ngDyn*(ngDyn+1)/2+3*ngdDyn*(ngdDyn+1)/2+
     >           ngDyn+ngDyn*norbit+Ngrid


c aca hay un bug para sistemas muy grandes ng2 da <0, Nick
c      write(*,*)ng2,ngDyn,ngdDyn,norbit,Ngrid

      allocate(X(ngDyn,ng3),XX(ngdDyn,ngdDyn))
      allocate(RMM(ng2))

      allocate (c(ngnu,nl),a(ngnu,nl),Nuc(ngnu),ncont(ngnu)
     > ,cx(ngdnu,nl),ax(ngdnu,nl),Nucx(ngdnu),ncontx(ngdnu)
     > ,cd(ngdnu,nl),ad(ngdnu,nl),Nucd(ngdnu),ncontd(ngdnu)
     > ,indexii(ngnu),indexiid(ngdnu))

      if (ecpmode) then
         allocate (Cnorm(ngnu,nl)) ! Cnorm contains normalized coefficient to 1 of basis, diff C for  x^2,y^2,z^2 and  xy,xz,yx (3^0.5 factor)
      end if


      allocate (r(ntatom,3),v(ntatom,3),rqm(natom,3),Em(ntatom)
     >,Rm(ntatom),pc(ntatom),Iz(natom),af(natom*ngd0),
     >  B(natom*ngd0,3))
      allocate (nnat(100))
      allocate(d(natom,natom))

      do i=1,natom
        read(101,*) iz(i),r(i,1:3)
        rqm(i,1:3)=r(i,1:3)
c       write(*,*) iz(i),r(i,1:3)
      enddo
      do i=natom+1,ntatom
        read(101,*) pc(i),r(i,1:3)   ! o es pc(i-natom)???
c       write(*,*) pc(i),r(i,1:3)
       enddo
       r=r/0.529177D0
       rqm=rqm/0.529177D0


       call g2g_init()   !initialize g2g

        nqnuc=0
       do i=1,natom
         nqnuc=nqnuc+Iz(i)
       enddo

       nco=((nqnuc - charge)-Nunp)/2


c       write(*,*) 'NCO=',NCO
c       write(*,*) natom,ntatom,ngDyn,ngdDyn,ng0,ngd0
c       write(*,*) ng2,ngDyn,ngdDyn
c--------------------------------------------------------
       call drive(ng2,ngDyn,ngdDyn)   !en esta rutina le pasa variables a g2g
       call lio_init()   !initialize lio
       call liomain()       !no hace nada!!!!!!
       if (.not.allocated(Smat))    allocate(Smat(M,M))
       if (.not.allocated(RealRho)) allocate(RealRho(M,M))
c--------------------------------------------------------

       if(OPEN) then
         if (ecpmode) stop "Lio havent got ECP for open shell"         
         call SCFOP(escf,dipxyz)
       else
         call SCF(escf,dipxyz)
       endif
c--------------------------------------------------------

c       write(*,*) 'SCF ENRGY=',escf

      if(writeforces) then
	if (ecpmode) stop "Lio havent got ECP forces"
       open(unit=123,file='fuerzas')
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
c       call g2g_solve_groups(3, Exc, dxyzqm)
c       write(*,*) dxyzqm

       do k=1,natom
       write(123,100)
     >     k,dxyzqm(k,1),dxyzqm(k,2),dxyzqm(k,3)
       enddo
         if(nsol.gt.0) then
          do k=natom,natom+nsol
!         write(123,'("fuerza",I,D,D,D)')
           write(123,100)
     >     k,dxyzcl(k,1),dxyzcl(k,2),dxyzcl(k,3)
          enddo

         endif
       deallocate (dxyzqm)
       if(nsol.gt.0) deallocate(dxyzcl)
       endif
       call lio_finalize()
100    format (I5,2x,f10.6,2x,f10.6,2x,f10.6)
       end program

