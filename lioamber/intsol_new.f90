!Esta rutina es una copia de intsol.f, pero reemplaza el array RMM por arrayas mas pequeños

!-------------------------------------------------------------------
! INTEGRAL CALCULATIONS FOR THE SOLVENT POINT CHARGES ELECTROSTATIC
! INTERACTIONS WITH THE SOLUTE ELECTRONIC DENSITY
!
! Dario Estrin
! Buenos Aires, August. 1994.
!
! 1 e integrals
! using the Obara-Saika recursive method.
! are added into the Fock matrix
!
! It's the same as int.f, but using the solvent atoms partial charges
!
!
!-------------------------------------------------------------------
      subroutine intsol_new(E1s,Ens,elec)

       use garcha_mod
!
      implicit real*8 (a-h,o-z)
      logical elec
      dimension xi(3)
!
!      real*8, dimension (:,:), ALLOCATABLE :: d
      real*8, dimension (:), ALLOCATABLE :: s0s,s1s,s2s,s3s,s4s
      dimension Q(3)
!,d(natom,natom),s0s(ntatom),s1s(ntatom),
!     > s2s(ntatom),s3s(ntatom),s4s(ntatom)
! distance between pairs of centers
!
!
!      allocate(d(natom,natom))
      allocate(s0s(ntatom),s1s(ntatom),s2s(ntatom),s3s(ntatom) &
       ,s4s(ntatom))
      if (NORM) then
      sq3=sqrt(3.D0)
      else
      sq3=1.D0
      endif
!
      do 1 l=1,3
 1     Ll(l)=l*(l-1)/2
!
      ns=nshell(0)
      np=nshell(1)
      nd=nshell(2)
      MM=M*(M+1)/2
      MMd=Md*(Md+1)/2
      M2=2*M
!
! Pointers
! first P
      M1=1
! now Pnew
      M3=M1+MM
! now S, F also uses the same position after S was used
      M11=M3+MM
! now G
      M7=M11+MM
! now Gm
      M9=M7+MMd
! now H
      M11=M9+MMd
!
      E1s=0.0D0
      Ens=0.0D0
      Ese=0.0D0
!      do 50 i=1,natom
!      do 50 j=1,natom
!       d(i,j)=(r(i,1)-r(j,1))**2+(r(i,2)-r(j,2))**2+
!     >        (r(i,3)-r(j,3))**2
! 50   continue

! SOLVENT-SOLUTE ELECTROSTATIC
!
!        ncerca2=0
!        nlejos=0
!        if(watermod.eq.1) nvecin=2 !tre sites
!        if(watermod.eq.2) nvecin=3 !four sites
!        if(watermod.eq.1.or.watermod.eq.2) then
!      do i=natom+1,natom+nsol
!          ncerca=0
!        do j=natom+1,natom+nsol

!
!          distx=(r(i,1)-r(j,1))**2
!          disty=(r(i,2)-r(j,2))**2
!          distz=(r(i,3)-r(j,3))**2
!          distint=distx+disty+distz
!          if (distint.lt.8.45D0) then
!            ncerca=ncerca+1
!          endif

!         enddo
!          if(ncerca.le.nvecin) then
!            pc(i)=0.D0
!             nlejos=nlejos+1
!           else
!           ncerca2=ncerca2+1
!           endif
!         enddo
!        write(*,*) 'ncerca2=',ncerca2,nlejos,watermod
!         endif


       do 125 j1=1,natom
       do 125 j2=natom+1,nsol+natom
!
!
       tx=r(j1,1)-r(j2,1)
       ty=r(j1,2)-r(j2,2)
       tz=r(j1,3)-r(j2,3)
       dd2=tx**2+ty**2+tz**2
       dd2=sqrt(dd2)
!
        Ens=Ens+Iz(j1)*pc(j2)/dd2

!
 125   continue


        if (elec)  then
!
!---------------------------------------------------------------
!
! first loop (s|s) case -------------------------------------------
!
      do 200 i=1,ns
      do 200 j=1,i
!
      dd=d(Nuc(i),Nuc(j))
!
      do 200 ni=1,ncont(i)
      do 200 nj=1,ncont(j)
!
! (0|0) calculation
      zij=a(i,ni)+a(j,nj)
      alf=a(i,ni)*a(j,nj)/zij
      Q(1)=(a(i,ni)*r(Nuc(i),1)+a(j,nj)*r(Nuc(j),1))/zij
      Q(2)=(a(i,ni)*r(Nuc(i),2)+a(j,nj)*r(Nuc(j),2))/zij
      Q(3)=(a(i,ni)*r(Nuc(i),3)+a(j,nj)*r(Nuc(j),3))/zij
      ccoef=c(i,ni)*c(j,nj)
!        write(88,333) i,j,c(i,ni),c(j,nj),Nuc(i),Nuc(j)
       rexp=alf*dd
      if(rexp.lt.rmax) then


!
      ss=pi32*exp(-alf*dd)/(zij*sqrt(zij))
      temp=2.D0*sqrt(zij/pi)*ss
!
      k=i+((M2-j)*(j-1))/2
!
       tna=0.D0
      do 202 j1=natom+1,nsol+natom
!
!
!
       tx=Q(1)-r(j1,1)
       ty=Q(2)-r(j1,2)
       tz=Q(3)-r(j1,3)
!
       u=tx**2 + ty**2 + tz**2
       u=u*zij

       s0s(j1)=pc(j1)*temp*FUNCT(0,u)
       tna=tna-s0s(j1)
 202   continue

!
      term=ccoef*tna
      Fock_Hcore(k)=Fock_Hcore(k)+ term

      E1s=E1s+P_density(k)*term
       if(E1s.ne.E1s) then
       write(*,*) 'E1s NaN 1'
        stop
       endif
        endif
 200  continue
!
!------------------------------------------------------------------
!
! second loop  (p|s) case
!
!
      do 300 i=ns+1,ns+np,3
      do 300 j=1,ns
!
      dd=d(Nuc(i),Nuc(j))
!
      do 300 ni=1,ncont(i)
      do 300 nj=1,ncont(j)
!
      zij=a(i,ni)+a(j,nj)
      Q(1)=(a(i,ni)*r(Nuc(i),1)+a(j,nj)*r(Nuc(j),1))/zij
      Q(2)=(a(i,ni)*r(Nuc(i),2)+a(j,nj)*r(Nuc(j),2))/zij
      Q(3)=(a(i,ni)*r(Nuc(i),3)+a(j,nj)*r(Nuc(j),3))/zij
      alf=a(i,ni)*a(j,nj)/zij
       rexp=alf*dd
      if(rexp.lt.rmax) then
      ss=pi32*exp(-alf*dd)/(zij*sqrt(zij))
       temp=2.D0*sqrt(zij/pi)*ss
!
! loop over nuclei, part common for all shell
      do 302 j1=natom+1,natom+nsol
!
!
!
!
       tx=Q(1)-r(j1,1)
       ty=Q(2)-r(j1,2)
       tz=Q(3)-r(j1,3)
!
!
!
       u= tx**2 +ty**2 +tz**2
       u=u*zij

       s0s(j1)=temp*FUNCT(0,u)
       s1s(j1)=temp*FUNCT(1,u)
!
 302  continue
!
      ccoef=c(i,ni)*c(j,nj)
!
! l2: different p in the p shell ( x,y,z respectively)
!
      do 305 l2=1,3
        ii=i+l2-1
        t1=Q(l2)-r(Nuc(i),l2)
! ii index , taking into account different components of the shell
        k=ii+((M2-j)*(j-1))/2
!
! loop over nuclei, specific part
       tna=0.D0
      do 303 j1=natom+1,natom+nsol
!
!
       t2=Q(l2)-r(j1,l2)
       term=t1*s0s(j1)-t2*s1s(j1)
       tna=tna-pc(j1)*term
!
 303  continue

        term=ccoef*tna
        Fock_Hcore(k)=Fock_Hcore(k)+term
        E1s=E1s+P_density(k)*term

 305    continue
!
      endif
 300  continue
!-------------------------------------------------------------------
!
! (p|p) case
!
      do 400 i=ns+1,ns+np,3
      do 400 j=ns+1,i,3
!
      dd=d(Nuc(i),Nuc(j))
!
      do 400 ni=1,ncont(i)
      do 400 nj=1,ncont(j)
!
      zij=a(i,ni)+a(j,nj)
      z2=2.D0*zij
      Q(1)=(a(i,ni)*r(Nuc(i),1)+a(j,nj)*r(Nuc(j),1))/zij
      Q(2)=(a(i,ni)*r(Nuc(i),2)+a(j,nj)*r(Nuc(j),2))/zij
      Q(3)=(a(i,ni)*r(Nuc(i),3)+a(j,nj)*r(Nuc(j),3))/zij
      alf=a(i,ni)*a(j,nj)/zij

       rexp=alf*dd
      if(rexp.lt.rmax) then
      ss=pi32*exp(-alf*dd)/(zij*sqrt(zij))
       temp=2.D0*sqrt(zij/pi)*ss
!
! loop over nuclei, part common for all shell
      do 402 j1=natom+1,natom+nsol
!
!

       xi(1)=Q(1)-r(j1,1)
       xi(2)=Q(2)-r(j1,2)
       xi(3)=Q(3)-r(j1,3)
!
       u=xi(1)**2+xi(2)**2+xi(3)**2
       u=u*zij

       s0s(j1)=temp*FUNCT(0,u)
       s1s(j1)=temp*FUNCT(1,u)
       s2s(j1)=temp*FUNCT(2,u)
!
 402  continue
!
!
      ccoef=c(i,ni)*c(j,nj)
!
! loop over partial charges ( specific part)
      do 403 j1=natom+1,natom+nsol
!
!
!
      do 406 l1=1,3
       t1=Q(l1)-r(Nuc(i),l1)
!
       t2=Q(l1)-r(j1,l1)
       p0s=t1*s0s(j1)-t2*s1s(j1)
       p1s=t1*s1s(j1)-t2*s2s(j1)

!
      lij=3
      if (i.eq.j) then
       lij=l1
      endif
!
       do 406 l2=1,lij
       t1=Q(l2)-r(Nuc(j),l2)
       t2=Q(l2)-r(j1,l2)
       tna=t1*p0s-t2*p1s
!
       if (l1.eq.l2) then
        tna=tna+(s0s(j1)-s1s(j1))/z2
       endif
!
        tna1=tna*pc(j1)
!
       ii=i+l1-1
       jj=j+l2-1
       k=ii+((M2-jj)*(jj-1))/2
       term=-tna1*ccoef
       Fock_Hcore(k)=Fock_Hcore(k)+term
       E1s=E1s+P_density(k)*term
 406  continue
!

 403   continue
! ---------------
      endif
 400  continue
!
!-------------------------------------------------------------------
! (d|s) case
!
      do 500 i=ns+np+1,M,6
      do 500 j=1,ns
!
      dd=d(Nuc(i),Nuc(j))
!
      do 500 ni=1,ncont(i)
      do 500 nj=1,ncont(j)
!
      zij=a(i,ni)+a(j,nj)
      z2=2.D0*zij
      Q(1)=(a(i,ni)*r(Nuc(i),1)+a(j,nj)*r(Nuc(j),1))/zij
      Q(2)=(a(i,ni)*r(Nuc(i),2)+a(j,nj)*r(Nuc(j),2))/zij
      Q(3)=(a(i,ni)*r(Nuc(i),3)+a(j,nj)*r(Nuc(j),3))/zij
      alf=a(i,ni)*a(j,nj)/zij
       rexp=alf*dd
      if(rexp.lt.rmax) then
      ss=pi32*exp(-alf*dd)/(zij*sqrt(zij))
       temp=2.D0*sqrt(zij/pi)*ss
!
! loop over partial charges, part common for all shell
      do 502 j1=natom+1,natom+nsol
!
!
       xi(1)=Q(1)-r(j1,1)
       xi(2)=Q(2)-r(j1,2)
       xi(3)=Q(3)-r(j1,3)
!
       u=xi(1)**2+xi(2)**2+xi(3)**2
       u=u*zij
!
       s0s(j1)=temp*FUNCT(0,u)
       s1s(j1)=temp*FUNCT(1,u)
       s2s(j1)=temp*FUNCT(2,u)
 502  continue
!
!
      ccoef=c(i,ni)*c(j,nj)
!
      do 503 j1=natom+1,natom+nsol
!
!
      do 506 l1=1,3
       t1=Q(l1)-r(Nuc(i),l1)
       t2=Q(l1)-r(j1,l1)
       p0s=t1*s0s(j1)-t2*s1s(j1)
       p1s=t1*s1s(j1)-t2*s2s(j1)
!
      do 506 l2=1,l1
!
       t1=Q(l2)-r(Nuc(i),l2)
       t2=Q(l2)-r(j1,l2)
       tna=t1*p0s-t2*p1s
!
       f1=1.D0
       if (l1.eq.l2) then
        tna=tna+(s0s(j1)-s1s(j1))/z2
        f1=sq3
       endif
!
       l12=l1*(l1-1)/2+l2
! ordering of d shell should be:
! xx,yx,yy,zx,zy,zz ( 11, 21, 22, 31, 32, 33 )
!
       ii=i+l12-1
!
       k=ii+((M2-j)*(j-1))/2
       cc=ccoef/f1
       term=-cc*tna*pc(j1)
       Fock_Hcore(k)=Fock_Hcore(k)+term
       E1s=E1s+P_density(k)*term
 506  continue
!
 503  continue
      endif
 500  continue
!-----------------------------------------------------------------
!
! (d|p) case
!
      do 600 i=ns+np+1,M,6
      do 600 j=ns+1,ns+np,3
!
      dd=d(Nuc(i),Nuc(j))
!
      do 600 ni=1,ncont(i)
      do 600 nj=1,ncont(j)
!
      zij=a(i,ni)+a(j,nj)
      z2=2.D0*zij
      Q(1)=(a(i,ni)*r(Nuc(i),1)+a(j,nj)*r(Nuc(j),1))/zij
      Q(2)=(a(i,ni)*r(Nuc(i),2)+a(j,nj)*r(Nuc(j),2))/zij
      Q(3)=(a(i,ni)*r(Nuc(i),3)+a(j,nj)*r(Nuc(j),3))/zij
      alf=a(i,ni)*a(j,nj)/zij
       rexp=alf*dd
       if(rexp.lt.rmax) then
       ss=pi32*exp(-alf*dd)/(zij*sqrt(zij))
       temp=2.D0*sqrt(zij/pi)*ss
!
! loop over nuclei, part common for all shell
      do 602 j1=natom+1,natom+nsol
!
!
       xi(1)=Q(1)-r(j1,1)
       xi(2)=Q(2)-r(j1,2)
       xi(3)=Q(3)-r(j1,3)
!
       u=xi(1)**2+xi(2)**2+xi(3)**2
       u=u*zij
!
       s0s(j1)=temp*FUNCT(0,u)
       s1s(j1)=temp*FUNCT(1,u)
       s2s(j1)=temp*FUNCT(2,u)
       s3s(j1)=temp*FUNCT(3,u)
 602  continue
!
!
      ccoef=c(i,ni)*c(j,nj)
!
      do 603 j1=natom+1,natom+nsol
!
!
      do 606 l1=1,3
       t1=Q(l1)-r(Nuc(i),l1)
       t2=Q(l1)-r(j1,l1)
       p0s=t1*s0s(j1)-t2*s1s(j1)
       p1s=t1*s1s(j1)-t2*s2s(j1)
       p2s=t1*s2s(j1)-t2*s3s(j1)
!
      do 606 l2=1,l1
       t1=Q(l2)-r(Nuc(i),l2)
       t2=Q(l2)-r(j1,l2)
       pj0s=t1*s0s(j1)-t2*s1s(j1)
       pj1s=t1*s1s(j1)-t2*s2s(j1)
!
       f1=1.D0
       d0s=t1*p0s-t2*p1s
       d1s=t1*p1s-t2*p2s
!
       if (l1.eq.l2) then
        f1=sq3
        d0s=d0s+(s0s(j1)-s1s(j1))/z2
        d1s=d1s+(s1s(j1)-s2s(j1))/z2
       endif
!
!
      do 606 l3=1,3
!
       t1=Q(l3)-r(Nuc(j),l3)
       t2=Q(l3)-r(j1,l3)
       tna=t1*d0s-t2*d1s
!
       if (l1.eq.l3) then
        tna=tna+(pj0s-pj1s)/z2
       endif
!
       if (l2.eq.l3) then
        tna=tna+(p0s-p1s)/z2
       endif
!
!
!
       l12=l1*(l1-1)/2+l2
       ii=i+l12-1
       jj=j+l3-1
!
       k=ii+((M2-jj)*(jj-1))/2
       cc=ccoef/f1
       term=-cc*tna*pc(j1)
       Fock_Hcore(k)=Fock_Hcore(k)+term
       E1s=E1s+P_density(k)*term
 606  continue
!
 603  continue
!
      endif
 600  continue
!
!-------------------------------------------------------------------
!
! (d|d) case
!
      do 700 i=ns+np+1,M,6
      do 700 j=ns+np+1,i,6
!
      dd=d(Nuc(i),Nuc(j))
!
      do 700 ni=1,ncont(i)
      do 700 nj=1,ncont(j)
!
      zij=a(i,ni)+a(j,nj)
      z2=2.D0*zij
      Q(1)=(a(i,ni)*r(Nuc(i),1)+a(j,nj)*r(Nuc(j),1))/zij
      Q(2)=(a(i,ni)*r(Nuc(i),2)+a(j,nj)*r(Nuc(j),2))/zij
      Q(3)=(a(i,ni)*r(Nuc(i),3)+a(j,nj)*r(Nuc(j),3))/zij
      alf=a(i,ni)*a(j,nj)/zij
       rexp=alf*dd
      if(rexp.lt.rmax) then
      ss=pi32*exp(-alf*dd)/(zij*sqrt(zij))
      temp=2.D0*sqrt(zij/pi)*ss
!
! loop over nuclei, part common for all shell
      do 702 j1=natom+1,natom+nsol
!
       xi(1)=Q(1)-r(j1,1)
       xi(2)=Q(2)-r(j1,2)
       xi(3)=Q(3)-r(j1,3)
!
       u=xi(1)**2+xi(2)**2+xi(3)**2
       u=u*zij
!
       s0s(j1)=temp*FUNCT(0,u)
       s1s(j1)=temp*FUNCT(1,u)
       s2s(j1)=temp*FUNCT(2,u)
       s3s(j1)=temp*FUNCT(3,u)
       s4s(j1)=temp*FUNCT(4,u)
 702  continue
!
!
      ccoef=c(i,ni)*c(j,nj)
!
! Loop over partial charges
      do 703 j1=natom+1,natom+nsol
!
!
!
      do 706 l1=1,3
       t1=Q(l1)-r(Nuc(i),l1)
       t2=Q(l1)-r(j1,l1)
       p0s=t1*s0s(j1)-t2*s1s(j1)
       p1s=t1*s1s(j1)-t2*s2s(j1)
       p2s=t1*s2s(j1)-t2*s3s(j1)
       p3s=t1*s3s(j1)-t2*s4s(j1)
!
      do 706 l2=1,l1
       f1=1.D0
       t1=Q(l2)-r(Nuc(i),l2)
       t2=Q(l2)-r(j1,l2)
       pj0s=t1*s0s(j1)-t2*s1s(j1)
       pj1s=t1*s1s(j1)-t2*s2s(j1)
       pj2s=t1*s2s(j1)-t2*s3s(j1)
!
       d0s=t1*p0s-t2*p1s
       d1s=t1*p1s-t2*p2s
       d2s=t1*p2s-t2*p3s
!
!
       if (l1.eq.l2) then
        f1=sq3
        d0s=d0s+(s0s(j1)-s1s(j1))/z2
        d1s=d1s+(s1s(j1)-s2s(j1))/z2
        d2s=d2s+(s2s(j1)-s3s(j1))/z2
       endif
!
!
      lij=3
      if (i.eq.j) then
       lij=l1
      endif
!
       do 706 l3=1,lij
!
       t1=Q(l3)-r(Nuc(j),l3)
       t2=Q(l3)-r(j1,l3)
!
       d0p=t1*d0s-t2*d1s
       d1p=t1*d1s-t2*d2s
!
       pi0p=t1*p0s-t2*p1s
       pi1p=t1*p1s-t2*p2s
       pj0p=t1*pj0s-t2*pj1s
       pj1p=t1*pj1s-t2*pj2s
!
       if (l1.eq.l3) then
        d0p=d0p+(pj0s-pj1s)/z2
        d1p=d1p+(pj1s-pj2s)/z2
        pi0p=pi0p+(s0s(j1)-s1s(j1))/z2
        pi1p=pi1p+(s1s(j1)-s2s(j1))/z2
       endif
!
       if (l2.eq.l3) then
        d0p=d0p+(p0s-p1s)/z2
        d1p=d1p+(p1s-p2s)/z2
        pj0p=pj0p+(s0s(j1)-s1s(j1))/z2
        pj1p=pj1p+(s1s(j1)-s2s(j1))/z2
       endif
!
!
      lk=l3
      if (i.eq.j) then
       lk=min(l3,Ll(l1)-Ll(l3)+l2)
      endif
!
      do 706 l4=1,lk
!
       f2=1.D0
       t1=Q(l4)-r(Nuc(j),l4)
       t2=Q(l4)-r(j1,l4)
       tna=t1*d0p-t2*d1p
!
       if (l4.eq.l1) then
        tna=tna+(pj0p-pj1p)/z2
       endif
!
       if (l4.eq.l2) then
        tna=tna+(pi0p-pi1p)/z2
       endif
!
       if (l4.eq.l3) then
        f2=sq3
        tna=tna+(d0s-d1s)/z2
       endif
!
       cc=ccoef/(f1*f2)
       term=-cc*pc(j1)*tna
!
       l12=l1*(l1-1)/2+l2
       l34=l3*(l3-1)/2+l4
       ii=i+l12-1
       jj=j+l34-1
!
       k=ii+((M2-jj)*(jj-1))/2
       Fock_Hcore(k)=Fock_Hcore(k)+term
       E1s=E1s+P_density(k)*term
!
 706  continue
 703  continue
      endif
 700  continue


       endif
      deallocate(s0s,s2s,s3s,s4s)

!
 333  format(2(I4,2x),2(F10.4,2x),2(I4,2x))
      return
      end subroutine intsol_new
!-------------------------------------------------------------------
