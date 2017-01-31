!This subroutine is a copy of int2.f but without use RMM array


!-------------------------------------------------------------------
! Integrals subroutine -Second part
! 2 e integrals, 2 index : density fitting functions
! All of them are calculated
! using the Obara-Saika recursive method.
!
!
! loop over all basis functions
! now the basis is supposed to be ordered according to the type,
! all s, then all p, then all d, .....
! inside each type, are ordered in shells
! px,py,pz , dx2,dxy,dyy,dzx,dzy,dzz, .....
!
! ns ... marker for end of s
! np ... marker for end of p
! nd ... marker for end of d
!
! r(Nuc(i),j) j component of position of nucleus i , j=1,3
! Input :  density basis
! Output: G matrix
! G matrix should be inverted,
! later on, for evaluating  Coulomb terms
!-----------------------------------------------------------------
      subroutine int2_new()
       use garcha_mod
!
      implicit real*8 (a-h,o-z)
       real*8, dimension(:), allocatable :: dgelss_temp
       real*8, dimension(Md) :: inv_work
       integer XXX(8*Md)
!
! aux . things
      dimension Q(3),aux(ngd),Det(2)
      real, dimension (:), ALLOCATABLE :: trabajo
!
      if (NORM) then
      sq3=sqrt(3.D0)
      else
      sq3=1.D0
      endif
!
!      do 50 i=1,natom
!      do 50 j=1,natom
!       d(i,j)=(r(i,1)-r(j,1))**2+(r(i,2)-r(j,2))**2+
!     >        (r(i,3)-r(j,3))**2
! 50   continue
!
      nsd=nshelld(0)
      npd=nshelld(1)
      ndd=nshelld(2)
      Md2=2*Md
      M2=2*M
      MM=M*(M+1)/2
      MMd=Md*(Md+1)/2

!
! pointers
      M1=1 ! first P
      M3=M1+MM ! now Pnew
      M5=M3+MM ! now S
      M7=M5+MM ! now G
      M9=M7+MMd ! now Gm
      M11=M9+MMd ! now H
      M13=M11+MM ! now F
      M15=M13+M ! auxiliar things

!
!
! end ------------------------------------------------
      do 1 k=1,MMd
 1      Density_fitting_G(k)=0.D0
!
!--- 2 index electron repulsion for density basis set
!
! first loop (s|s) case -------------------------------------------
!
      do 200 i=1,nsd
      do 200 j=1,i
!
      dd=d(Nucd(i),Nucd(j))
!
      do 200 ni=1,ncontd(i)
      do 200 nj=1,ncontd(j)
!
! (0|0) calculation
      zij=ad(i,ni)+ad(j,nj)
      t0=ad(i,ni)*ad(j,nj)
      alf=t0/zij
      t1=sqrt(zij)*t0

      u=alf*dd
      ccoef=cd(i,ni)*cd(j,nj)
!
      s0s=pi5/t1*FUNCT(0,u)
!
      k=i+((Md2-j)*(j-1))/2
      Density_fitting_G(k)=Density_fitting_G(k)+ccoef*s0s
!
 200  continue
!
!------------------------------------------------------------------
!
! second loop  (p|s) case
!
!
      do 300 i=nsd+1,nsd+npd,3
      do 300 j=1,nsd
!
      dd=d(Nucd(i),Nucd(j))
!
      do 300 ni=1,ncontd(i)
      do 300 nj=1,ncontd(j)
!
      zij=ad(i,ni)+ad(j,nj)
      t0=ad(i,ni)*ad(j,nj)
      alf=t0/zij
      t1=sqrt(zij)*t0
      t2=pi5/t1
!
      Q(1)=(ad(i,ni)*r(Nucd(i),1)+ad(j,nj)*r(Nucd(j),1))/zij
      Q(2)=(ad(i,ni)*r(Nucd(i),2)+ad(j,nj)*r(Nucd(j),2))/zij
      Q(3)=(ad(i,ni)*r(Nucd(i),3)+ad(j,nj)*r(Nucd(j),3))/zij
!
      u=alf*dd
      s1s=t2*FUNCT(1,u)
!
!
      ccoef=cd(i,ni)*cd(j,nj)
!
! l2: different p in the p shell ( x,y,z respectively)
!
      do 305 l1=1,3
        t1=Q(l1)-r(Nucd(i),l1)
        tn=t1*s1s
!
        iii=i+l1-1
! ii index , taking into account different components of the shell
!
        k=iii+((Md2-j)*(j-1))/2
        Density_fitting_G(k)=Density_fitting_G(k)+tn*ccoef
 305   continue
 300   continue
!
!------------------------------------------------------------
!
! (p|p) case
!
      do 400 i=nsd+1,nsd+npd,3
      do 400 j=nsd+1,i,3
!
      dd=d(Nucd(i),Nucd(j))
!
      do 400 ni=1,ncontd(i)
      do 400 nj=1,ncontd(j)
!
      zij=ad(i,ni)+ad(j,nj)
      z2=2.D0*zij
      t0=ad(i,ni)*ad(j,nj)
      alf=t0/zij
      t1=sqrt(zij)*t0
      t2=pi5/t1
!
      Q(1)=(ad(i,ni)*r(Nucd(i),1)+ad(j,nj)*r(Nucd(j),1))/zij
      Q(2)=(ad(i,ni)*r(Nucd(i),2)+ad(j,nj)*r(Nucd(j),2))/zij
      Q(3)=(ad(i,ni)*r(Nucd(i),3)+ad(j,nj)*r(Nucd(j),3))/zij
!
      u=alf*dd
      s1s=t2*FUNCT(1,u)
      s2s=t2*FUNCT(2,u)

      ccoef=cd(i,ni)*cd(j,nj)
!
      do 405 l1=1,3
       t1=Q(l1)-r(Nucd(i),l1)
       ps=t1*s2s
!
       do 405 l2=1,3
!
       t1=Q(l2)-r(Nucd(j),l2)
       tn=t1*ps
!
       if (l1.eq.l2) then
        tn=tn+s1s/z2
       endif
!
       iii=i+l1-1
       jj=j+l2-1
!
!      this to convert to 1 dimensional array, in diagonal case
!      we calculate more things than necessary . They should be
!      eliminated
       if(iii.ge.jj) then
       k=iii+((Md2-jj)*(jj-1))/2
       Density_fitting_G(k)=Density_fitting_G(k)+tn*ccoef
       endif
 405  continue
!
 400  continue
!-------------------------------------------------------------------
! (d|s) case
      do 500 i=nsd+npd+1,Md,6
      do 500 j=1,nsd
!
      dd=d(Nucd(i),Nucd(j))
!
      do 500 ni=1,ncontd(i)
      do 500 nj=1,ncontd(j)
!
      zij=ad(i,ni)+ad(j,nj)
      t0=ad(i,ni)*ad(j,nj)
      alf=t0/zij
      roz=ad(j,nj)/zij
      t1=sqrt(zij)*t0
      t2=pi5/t1
!
      Q(1)=(ad(i,ni)*r(Nucd(i),1)+ad(j,nj)*r(Nucd(j),1))/zij
      Q(2)=(ad(i,ni)*r(Nucd(i),2)+ad(j,nj)*r(Nucd(j),2))/zij
      Q(3)=(ad(i,ni)*r(Nucd(i),3)+ad(j,nj)*r(Nucd(j),3))/zij
!
      u=alf*dd
      s0s=t2*FUNCT(0,u)
      s1s=t2*FUNCT(1,u)
      s2s=t2*FUNCT(2,u)

      ccoef=cd(i,ni)*cd(j,nj)
!
      do 505 l1=1,3
       t1=Q(l1)-r(Nucd(i),l1)
       ps=t1*s2s
!
       do 505 l2=1,l1
!
       t1=Q(l2)-r(Nucd(i),l2)
       tn=t1*ps
!
       f1=1.
       if (l1.eq.l2) then
        tn=tn+(s0s-roz*s1s)/(2.D0*ad(i,ni))
        f1=sq3
       endif
!
       l12=l1*(l1-1)/2+l2
       iii=i+l12-1
!
       cc=ccoef/f1
       k=iii+((Md2-j)*(j-1))/2
       Density_fitting_G(k)=Density_fitting_G(k)+tn*cc
 505  continue
!
 500  continue
!-------------------------------------------------------------------
! (d|p) case
      do 600 i=nsd+npd+1,Md,6
      do 600 j=nsd+1,nsd+npd,3
!
      dd=d(Nucd(i),Nucd(j))
!
      do 600 ni=1,ncontd(i)
      do 600 nj=1,ncontd(j)
!
      zij=ad(i,ni)+ad(j,nj)
      z2=2.D0*zij
      t0=ad(i,ni)*ad(j,nj)
      alf=t0/zij
      t1=sqrt(zij)*t0
      t2=pi5/t1
!
      Q(1)=(ad(i,ni)*r(Nucd(i),1)+ad(j,nj)*r(Nucd(j),1))/zij
      Q(2)=(ad(i,ni)*r(Nucd(i),2)+ad(j,nj)*r(Nucd(j),2))/zij
      Q(3)=(ad(i,ni)*r(Nucd(i),3)+ad(j,nj)*r(Nucd(j),3))/zij
!
      u=alf*dd
      s1s=t2*FUNCT(1,u)
      s2s=t2*FUNCT(2,u)
      s3s=t2*FUNCT(3,u)
!
      ccoef=cd(i,ni)*cd(j,nj)
!
      do 605 l1=1,3
       t1=Q(l1)-r(Nucd(i),l1)
       pis=t1*s2s
       pi2s=t1*s3s
!
       do 605 l2=1,l1
!
       t1=Q(l2)-r(Nucd(i),l2)
       pjs=t1*s2s
!
       ds=t1*pi2s
!
       f1=1.D0
       if (l1.eq.l2) then
        f1=sq3
        ds=ds+(s1s-alf*s2s/ad(i,ni))/(2.D0*ad(i,ni))
       endif

! index of p
!
       do 605 l3=1,3
!
       t0=Q(l3)-r(Nucd(j),l3)
       tn=t0*ds
!
       if (l1.eq.l3) then
        tn=tn+pjs/z2
       endif
!
       if (l2.eq.l3) then
        tn=tn+pis/z2
       endif
!
!
       l12=l1*(l1-1)/2+l2
       iii=i+l12-1
       jj=j+l3-1
!
       cc=ccoef/f1
!
       k=iii+((Md2-jj)*(jj-1))/2
       Density_fitting_G(k)=Density_fitting_G(k)+tn*cc
 605  continue
!
 600  continue
!
!-------------------------------------------------------------------
! (d|d) case
      do 700 i=nsd+npd+1,Md,6
      do 700 j=nsd+npd+1,i,6
!
      dd=d(Nucd(i),Nucd(j))
!
      do 700 ni=1,ncontd(i)
      do 700 nj=1,ncontd(j)
!
      zij=ad(i,ni)+ad(j,nj)
      z2=2.D0*zij
      za=2.D0*ad(i,ni)
      zc=2.D0*ad(j,nj)
      t0=ad(i,ni)*ad(j,nj)
      alf=t0/zij
      t1=sqrt(zij)*t0
      t2=pi5/t1
!
      ti=ad(i,ni)/zij
      tj=ad(j,nj)/zij
      Q(1)=ti*r(Nucd(i),1)+tj*r(Nucd(j),1)
      Q(2)=ti*r(Nucd(i),2)+tj*r(Nucd(j),2)
      Q(3)=ti*r(Nucd(i),3)+tj*r(Nucd(j),3)
!
      u=alf*dd
      s0s=t2*FUNCT(0,u)
      s1s=t2*FUNCT(1,u)
      s2s=t2*FUNCT(2,u)
      s3s=t2*FUNCT(3,u)
      s4s=t2*FUNCT(4,u)
!
      t3=(s0s-tj*s1s)/za
      t4=(s1s-tj*s2s)/za
      t5=(s2s-tj*s3s)/za
      ccoef=cd(i,ni)*cd(j,nj)
!
      do 705 l1=1,3
       t1=Q(l1)-r(Nucd(i),l1)
       pis=t1*s2s
       pi2s=t1*s3s
       pi3s=t1*s4s
!
       do 705 l2=1,l1
!
       t1=Q(l2)-r(Nucd(i),l2)
       pjs=t1*s2s
       pj2s=t1*s3s
!
       ds=t1*pis
       d1s=t1*pi2s
       d2s=t1*pi3s
!
       f1=1.D0
       if (l1.eq.l2) then
        ds=ds+t3
        d1s=d1s+t4
        d2s=d2s+t5
        f1=sq3
       endif

!
       t6=(ds-ti*d1s)/zc
!
      lij=3
      if (i.eq.j) then
       lij=l1
      endif
       do 705 l3=1,lij
!
       t0=Q(l3)-r(Nucd(j),l3)
       dp=t0*d2s
       pip=t0*pi2s
       pjp=t0*pj2s
!
       if (l1.eq.l3) then
        dp=dp+pj2s/z2
        pip=pip+s2s/z2
       endif
!
       if (l2.eq.l3) then
        dp=dp+pi2s/z2
        pjp=pjp+s2s/z2
       endif
!
!
      lk=l3
      if (i.eq.j) then
       lll=l1*(l1-1)/2-l3*(l3-1)/2+l2
       lk=min(l3,lll)
      endif
       do 705 l4=1,lk
!
       t0=Q(l4)-r(Nucd(j),l4)
       tn=t0*dp
!
       if (l1.eq.l4) then
        tn=tn+pjp/z2
       endif
!
       if (l2.eq.l4) then
        tn=tn+pip/z2
       endif
!
       f2=1.D0
       if (l3.eq.l4) then
        tn=tn+t6
        f2=sq3
       endif
!
       l12=l1*(l1-1)/2+l2
       l34=l3*(l3-1)/2+l4
       iii=i+l12-1
       jj=j+l34-1
!
       cc=ccoef/(f1*f2)
!
!      this to convert to 1 dimensional array, in diagonal case
       k=iii+(Md2-jj)*(jj-1)/2
!
       Density_fitting_G(k)=Density_fitting_G(k)+tn*cc
 705  continue
!
 700  continue
!
!
!
       do 216 i=1,Md
       do 216 j=1,Md
!
        if(i.ge.j) then
         k=i+(Md*2-j)*(j-1)/2
         else
         k=j+(Md*2-i)*(i-1)/2
        endif
        XX(i,j)=Density_fitting_G(k)
!	write(25,*) Density_fitting_G(k) hasta aca esta bien, Nick
 216   continue
!
!
       kk=0
      do 112 j=1,Md
      do 112 i=j,Md
       kk=kk+1
       tmp=XX(i,j)
       Density_fitting_G(kk)=tmp
 112  continue
!
      MMp=Md*(Md+1)/2
      do 199 k=1,MMp
 199   Density_fitting_Gm(k)=0.0D0


!Test, Nick hasta aca ok
!	do k=1,MMp
!	   write(65,*) Density_fitting_G(k), Density_fitting_Gm(k)
!	end do


!
!     M10=M9+Md
      M10=M9+MMd+MM+1  !esto seria la posicion M13 definida en SCF contiene W, Nick
      M12=M10+Md !no cae en ninguna posicion definida en SCF
      Md3=3*Md
! ESSL OPTION ------------------------------
#ifdef essl
      CALL DGESVF(10,XX,Md,Density_fitting_Gm,Md,1,RMM(M10), &
                   Md,Md,RMM(M12),Md3)
	write(*,*) "int2 test1"
!por aca no esta pasando, testear luego, Nick

       ss=RMM(M10)/RMM(M10+Md-1)
!
#endif
!
! LAPACK OPTION ------------------------------
!
#ifdef pack

!
       do i=1,Md
        aux(i)=0.0D0
       enddo
       Md5=5*Md
      rcond=1.0D-07

!
! CH - why call dgelss here? We only want the singular values - couldn't just
! something like dgesvd be called without calculating singular vectors?
!
!      call dgelss(Md,Md,1,XX,Md,aux,Md,RMM(M9),rcond,irank,RMM(M10),
!     >            -1,info)
!      Md5=RMM(M10)
!      allocate(dgelss_temp(Md5))
!      call dgelss(Md,Md,1,XX,Md,aux,Md,RMM(M9),rcond,irank,dgelss_temp,
!     >            Md5,info)
!      deallocate(dgelss_temp)


!Test, Nick aca pincho en Density_fitting_Gm
       do k=1,MMp
          write(63,*) Density_fitting_G(k), Density_fitting_Gm(k)
       end do

      call g2g_timer_sum_start('G condition') 
#ifdef  magma
	write(*,*) "int2 mal"
      call magmaf_dgesdd('N',Md,Md,XX,Md,RMM(M9),0,1,0,1, &
                  RMM(M10),-1,XXX,info)
#else
	write(*,*) "int2 test2"
      call dgesdd('N',Md,Md,XX,Md,Density_fitting_Gm(1),0,1,0,1, &
                  RMM(M10),-1,XXX,info) 

#endif

!Test, Nick aca pincho en Density_fitting_Gm
       do k=1,MMp
          write(64,*) Density_fitting_G(k), Density_fitting_Gm(k)
       end do



      Md5=RMM(M10)
      allocate(dgelss_temp(Md5)) 
#ifdef  magma
	write(*,*) "int2 mal"
      call magmaf_dgesdd('N',Md,Md,XX,Md,RMM(M9),0,1,0,1, &
                  dgelss_temp,Md5,XXX,info)
#else
	write(*,*) "int2, test3"

      call dgesdd('N',Md,Md,XX,Md,Density_fitting_Gm(1),0,1,0,1, &
                  dgelss_temp,Md5,XXX,info) 
#endif
      deallocate(dgelss_temp)

!Test, Nick aca pincho en Density_fitting_Gm
       do k=1,MMp
          write(65,*) Density_fitting_G(k), Density_fitting_Gm(k)
       end do


!      ss=RMM(M9)/RMM(M9+Md-1)
       ss=Density_fitting_Gm(1)/Density_fitting_Gm(Md-1)
	write(*,*) "ss vale", ss
!
#endif
!      write (*,*) ss, "criterio ajuste base auxiliar, Nick"
       if (ss.gt.1.D14) then
        SVD=.true.
      stop "trata de usar SVD"
       endif

      call g2g_timer_sum_stop('G condition')
!
!------------------------------
! inversion of G matrix , kept in Gm
!
      if (SVD) then
       write(*,900) ss
       call aint_query_gpu_level(igpu)
       if (igpu.eq.5) then
         write(*,*) "G IS ILL-CONDITIONED"
         write(*,*) "THE SVD AUXILIARY DENSITY FIT IS NOT SUPPORTED"
         write(*,*) "IN THE GPU VERSION OF LIO"
         stop
       endif

      else
!
!
!
!
! LINPACK OPTION
#ifdef pack
!
      call g2g_timer_sum_start('G invert')

       do i=1,Md
       do j=1,Md

        if(i.ge.j) then
         k=i+(Md*2-j)*(j-1)/2
         else
         k=j+(Md*2-i)*(i-1)/2
        endif
        XX(i,j)=Density_fitting_G(k)
      enddo
      enddo

!      kk=0
!      do 313 j=1,Md
!       do 313 i=1,j
!       kk=kk+1
!       kx=M7+j+(2*Md-i)*(i-1)/2-1
!       RMM(M9+kk-1)=RMM(kx)
!
! 313  continue
!

!      call dppco(RMM(M9),Md,rcond,aux,info)
      call dsytrf('U',Md,XX,Md,XXX,RMM(M10),-1,info)
      Md5=RMM(M10)
      allocate(dgelss_temp(Md5))
      call dsytrf('U',Md,XX,Md,XXX,dgelss_temp,Md5,info)
      deallocate(dgelss_temp)

      call dsytri('U',Md,XX,Md,XXX,inv_work,info)

      do i=1,Md
      do j=1,i
        k=i+(Md*2-j)*(j-1)/2
        Density_fitting_G(k) = XX(j,i)
      enddo
      enddo

      call g2g_timer_sum_stop('G invert') 
#endif
!
      endif
 900  format('SWITCHING TO SVD rcond=',D10.3)
!
!-------------------------------------------------------------------

       do k=1,MMp
          write(70,*) Density_fitting_G(k), Density_fitting_Gm(k)
       end do


      return
      end subroutine int2_new

