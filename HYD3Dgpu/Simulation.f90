      module commons
      implicit none
      integer::nhy
      integer,parameter:: nhymax =600000
      real(8)::time,dt
      data time / 0.0d0 /
      real(8),parameter:: timemax=5.0d0
      real(8),parameter:: dtout=5.0d0/600
      
      integer,parameter::ngrid=128
      integer,parameter::mgn=2
      integer,parameter::in=ngrid+2*mgn+1 &
     &                  ,jn=ngrid+2*mgn+1 &
     &                  ,kn=ngrid+2*mgn+1
      integer,parameter::is=mgn+1 &
     &                  ,js=mgn+1 &
     &                  ,ks=mgn+1
      integer,parameter::ie=ngrid+mgn &
     &                  ,je=ngrid+mgn &
     &                  ,ke=ngrid+mgn

      real(8),parameter:: x1min=-0.5d0,x1max=0.5d0
      real(8),parameter:: x2min=-0.5d0,x2max=0.5d0
      real(8),parameter:: x3min=-0.5d0,x3max=0.5d0
      real(8),dimension(in)::x1a,x1b
      real(8),dimension(jn)::x2a,x2b
      real(8),dimension(kn)::x3a,x3b

      real(8),dimension(in,jn,kn)::d,et,mv1,mv2,mv3
      real(8),dimension(in,jn,kn)::p,ei,v1,v2,v3,cs

      real(8),parameter::gam=5.0d0/3.0d0

!$acc declare create(ngrid,mgn)
!$acc declare create(in,jn,kn)
!$acc declare create(is,js,ks)
!$acc declare create(ie,je,ke)
      
!$acc declare create(dt)
!$acc declare create(x1a,x1b)
!$acc declare create(x2a,x2b)
!$acc declare create(x3a,x3b)
      
!$acc declare create(d,et,mv1,mv2,mv3)
!$acc declare create(p,ei,v1,v2,v3,cs)

!$acc declare create(gam)
      end module commons
      
      module fluxmod
      use commons, only : in,jn,kn
      implicit none
      integer,parameter::nden=1,nve1=2,nve2=3,nve3=4,nene=5,npre=6
      integer,parameter::nhyd=6
      real(8),dimension(nhyd,in,jn,kn):: svc

      integer,parameter::mudn=1,muvu=2,muvv=3,muvw=4,muet=5  &
     &                  ,mfdn=6,mfvu=7,mfvv=8,mfvw=9,mfet=10 &
     &                  ,mcsp=11,mvel=12,mpre=13
      integer,parameter:: mflx=5,madd=3

      integer,parameter:: mden=1,mrv1=2,mrv2=3,mrv3=4,meto=5 &
     &                          ,mrvu=muvu,mrvv=muvv,mrvw=muvw
      real(8),dimension(mflx,in,jn,kn):: nflux1,nflux2,nflux3

!$acc declare create(svc,nflux1,nflux2,nflux3)
      end module fluxmod

      program main
      use commons
      implicit none
      write(6,*) "setup grids and fields"
      call GenerateGrid
      call GenerateProblem
      call ConsvVariable
      write(6,*) "entering main loop"
! main loop
      mloop: do nhy=1,80000
         call TimestepControl
         if(mod(nhy,100) .eq. 0 )write(6,*)nhy,time,dt
         call BoundaryCondition
         call StateVevtor
         call NumericalFlux1
         call NumericalFlux2
         call NumericalFlux3
         call UpdateConsv
         call PrimVariable
         time=time+dt
         call Output
         if(time > timemax) exit mloop         
      enddo mloop

      write(6,*) "program has been finished"
      end program main

      subroutine GenerateGrid
      use commons
      implicit none
      real(8)::dx,dy,dz
      integer::i,j,k

      dx=(x1max-x1min)/ngrid
      do i=1,in
         x1a(i) = dx*(i-(mgn+1))+x1min
      enddo
      do i=1,in-1
         x1b(i) = 0.5d0*(x1a(i+1)+x1a(i))
      enddo

      dy=(x2max-x2min)/ngrid
      do j=1,jn
         x2a(j) = dy*(j-(mgn+1))+x2min
      enddo
      do j=1,jn-1
         x2b(j) = 0.5d0*(x2a(j+1)+x2a(j))
      enddo

      dz=(x3max-x3min)/ngrid
      do k=1,kn
         x3a(k) = dz*(k-(mgn+1))+x3min
      enddo
      do k=1,kn-1
         x3b(k) = 0.5d0*(x3a(k+1)+x3a(k))
      enddo

!$acc update device (x1a,x1b)
!$acc update device (x2a,x2b)
!$acc update device (x3a,x3b)

      return
      end subroutine GenerateGrid

      subroutine GenerateProblem
      use commons
      implicit none
      integer::i,j,k

      real(8)::pi

      real(8)::Ahl,Bhl,Chl
      real(8),parameter:: k_ini=4.0d0
      real(8),parameter:: ekin = 2.0d0
      real(8),parameter:: emag = 0.0d0
      real(8),parameter:: eint = 1.0d0
      real(8),parameter:: v0 = sqrt(ekin*2.0)
      real(8),parameter:: b0 = sqrt(emag*2.0)
      real(8),parameter:: p0 = eint/(gam-1.0d0)
      real(8),parameter:: eps = 1.0d-1

      integer::seedsize
      integer,allocatable:: seed(:)
      real(8)::x

      call random_seed(size=seedsize)
      write(6,*)"seed size",seedsize
      allocate(seed(seedsize))  
      call random_seed(get=seed)

      pi=acos(-1.0d0)

      Ahl = 0.5d0
      Bhl = 0.5d0
      Chl = 0.5d0

      d(:,:,:) = 1.0d0

      do k=ks,ke
      do j=js,je
      do i=is,ie
         v1(i,j,k) = v0*(  Ahl*sin(k_ini*x3b(k)*2.0d0*pi/(x3max-x3min)) &
   &                     + Chl*cos(k_ini*x2b(j)*2.0d0*pi/(x2max-x2min)))
         v2(i,j,k) = v0*(  Bhl*sin(k_ini*x1b(i)*2.0d0*pi/(x1max-x1min)) &
   &                     + Ahl*cos(k_ini*x3b(k)*2.0d0*pi/(x3max-x3min)))
         v3(i,j,k) = v0*(  Chl*sin(k_ini*x2b(j)*2.0d0*pi/(x2max-x2min)) &
   &                     + Bhl*cos(k_ini*x1b(i)*2.0d0*pi/(x1max-x1min)))
          p(i,j,k) = p0

         call random_number(x)
         v1(i,j,k) = v1(i,j,k)*(1.0d0+eps*(x-0.5d0))
         call random_number(x)
         v2(i,j,k) = v2(i,j,k)*(1.0d0+eps*(x-0.5d0))
         call random_number(x)
         v3(i,j,k) = v3(i,j,k)*(1.0d0+eps*(x-0.5d0))
      enddo
      enddo
      enddo
   
      do k=ks,ke
      do j=js,je
      do i=is,ie
         call random_number(x)
          ei(i,j,k) = p(i,j,k)/(gam-1.0d0)
          cs(i,j,k) = sqrt(gam*p(i,j,k)/d(i,j,k))
      enddo
      enddo
      enddo
      

      write(6,*)"initial profile is set"
      call BoundaryCondition

!$acc update device (d,v1,v2,v3)
!$acc update device (p,ei,cs)

      return
      end subroutine GenerateProblem

      subroutine BoundaryCondition
      use commons
      implicit none
      integer::i,j,k

!$acc kernels    
      do k=1,kn-1
      do j=1,jn-1
      do i=1,mgn
           d(i,j,k) =  d(ie-mgn+i,j,k)
          ei(i,j,k) = ei(ie-mgn+i,j,k)
          v1(i,j,k) = v1(ie-mgn+i,j,k)
          v2(i,j,k) = v2(ie-mgn+i,j,k)
          v3(i,j,k) = v3(ie-mgn+i,j,k)
      enddo
      enddo
      enddo

      do k=1,kn-1
      do j=1,jn-1
      do i=1,mgn
           d(ie+i,j,k) =  d(is+i-1,j,k)
          ei(ie+i,j,k) = ei(is+i-1,j,k)
          v1(ie+i,j,k) = v1(is+i-1,j,k)
          v2(ie+i,j,k) = v2(is+i-1,j,k)
          v3(ie+i,j,k) = v3(is+i-1,j,k)
      enddo
      enddo
      enddo

      do k=1,kn-1
      do i=1,in-1
      do j=1,mgn
           d(i,j,k) =  d(i,je-mgn+j,k)
          ei(i,j,k) = ei(i,je-mgn+j,k)
          v1(i,j,k) = v1(i,je-mgn+j,k)
          v2(i,j,k) = v2(i,je-mgn+j,k)
          v3(i,j,k) = v3(i,je-mgn+j,k)
      enddo
      enddo
      enddo

      do k=1,kn-1
      do i=1,in-1
      do j=1,mgn
           d(i,je+j,k) =  d(i,js+j-1,k)
          ei(i,je+j,k) = ei(i,js+j-1,k)
          v1(i,je+j,k) = v1(i,js+j-1,k)
          v2(i,je+j,k) = v2(i,js+j-1,k)
          v3(i,je+j,k) = v3(i,js+j-1,k)
      enddo
      enddo
      enddo

      do j=1,jn-1
      do i=1,in-1
      do k=1,mgn
           d(i,j,k) =  d(i,j,ke-mgn+k)
          ei(i,j,k) = ei(i,j,ke-mgn+k)
          v1(i,j,k) = v1(i,j,ke-mgn+k)
          v2(i,j,k) = v2(i,j,ke-mgn+k)
          v3(i,j,k) = v3(i,j,ke-mgn+k)
      enddo
      enddo
      enddo

      do j=1,jn-1
      do i=1,in-1
      do k=1,mgn
           d(i,j,ke+k) =  d(i,j,ks+k-1)
          ei(i,j,ke+k) = ei(i,j,ks+k-1)
          v1(i,j,ke+k) = v1(i,j,ks+k-1)
          v2(i,j,ke+k) = v2(i,j,ks+k-1)
          v3(i,j,ke+k) = v3(i,j,ks+k-1)
      enddo
      enddo
      enddo
!$acc end kernels

      return
      end subroutine BoundaryCondition

      subroutine ConsvVariable
      use commons
      implicit none
      integer::i,j,k

!$acc kernels      
!$acc loop independent

      do k=ks,ke
      do j=js,je
      do i=is,ie
          et(i,j,k) = 0.5d0*d(i,j,k)*(   &
     &                    +v1(i,j,k)**2  &
     &                    +v2(i,j,k)**2  &
     &                    +v3(i,j,k)**2) &
     &                    +ei(i,j,k)
          mv1(i,j,k) =d(i,j,k)*v1(i,j,k)
          mv2(i,j,k) =d(i,j,k)*v2(i,j,k)
          mv3(i,j,k) =d(i,j,k)*v3(i,j,k)
      enddo
      enddo
      enddo
!$acc end kernels
      
      return
      end subroutine Consvvariable

      subroutine PrimVariable
      use commons
      implicit none
      integer::i,j,k

!$acc kernels      
!$acc loop independent
      do k=ks,ke
      do j=js,je
      do i=is,ie
          v1(i,j,k) = mv1(i,j,k)/d(i,j,k)
          v2(i,j,k) = mv2(i,j,k)/d(i,j,k)
          v3(i,j,k) = mv3(i,j,k)/d(i,j,k)

          ei(i,j,k) =  et(i,j,k)          &
     &          -0.5d0*d(i,j,k)*(         &
     &                    +v1(i,j,k)**2   &
     &                    +v2(i,j,k)**2   &
     &                    +v3(i,j,k)**2)

           p(i,j,k) =  ei(i,j,k)*(gam-1.0d0)
          cs(i,j,k) =  sqrt(gam*p(i,j,k)/d(i,j,k))
      enddo
      enddo
      enddo
!$acc end kernels
      return
      end subroutine PrimVariable

      subroutine TimestepControl
      use commons
      implicit none
      real(8)::dtl1
      real(8)::dtl2
      real(8)::dtl3
      real(8)::dtlocal
      real(8)::dtmin
      integer::i,j,k
!$acc kernels    
      dtmin=1.0d90
!$acc loop reduction(min:dtmin)
      do k=ks,ke
      do j=js,je
      do i=is,ie
         dtl1 =(x1a(i+1)-x1a(i))/(abs(v1(i,j,k)) +cs(i,j,k))
         dtl2 =(x2a(j+1)-x2a(j))/(abs(v2(i,j,k)) +cs(i,j,k))
         dtl3 =(x3a(k+1)-x3a(k))/(abs(v3(i,j,k)) +cs(i,j,k))
         dtlocal = min (dtl1,dtl2,dtl3)
         if(dtlocal .lt. dtmin) dtmin = dtlocal
      enddo
      enddo
      enddo

      dt = 0.05d0 * dtmin
!$acc end kernels
!$acc update host (dt)

!      write(6,*)"dt",dt
      return
      end subroutine TimestepControl

      subroutine StateVevtor
      use commons
      use fluxmod
      implicit none
      integer::i,j,k

!$acc kernels
!$acc loop independent
      do k=1,kn-1
      do j=1,jn-1
      do i=1,in-1
         svc(nden,i,j,k) =  d(i,j,k)
         svc(nve1,i,j,k) = v1(i,j,k)
         svc(nve2,i,j,k) = v2(i,j,k)
         svc(nve3,i,j,k) = v3(i,j,k)
         svc(nene,i,j,k) = ei(i,j,k)/d(i,j,k)
         svc(npre,i,j,k) = ei(i,j,k)*(gam-1.0d0)
      enddo
      enddo
      enddo
!$acc end kernels

      return
      end subroutine StateVevtor

      subroutine minmod(a,b,d)
!$acc routine seq
      use fluxmod, only : nhyd
      implicit none
      real(8),dimension(nhyd),intent(in)::a,b
      real(8),dimension(nhyd),intent(out)::d
      integer:: n

      do n=1,nhyd
         d(n) = sign(1.0d0,a(n))*max(0.0d0,min(abs(a(n))                &
     &                                        ,sign(1.0d0,a(n))*b(n)))
      enddo

      return
      end subroutine minmod


      subroutine vanLeer(dvp,dvm,dv)
!$acc routine seq
      use fluxmod, only : nhyd
      implicit none
      real(8),dimension(nhyd),intent(in)::dvp,dvm
      real(8),dimension(nhyd),intent(out)::dv
      integer:: n

      do n=1,nhyd
         if(dvp(n)*dvm(n) .gt. 0.0d0)then
            dv(n) =2.0d0*dvp(n)*dvm(n)/(dvp(n)+dvm(n))
         else
            dv(n) = 0.0d0
         endif

      enddo

      return
      end subroutine vanLeer



      subroutine MClimiter(a,b,c,d)
!$acc routine seq
      use fluxmod, only : nhyd
      implicit none
      real(8),dimension(nhyd),intent(in)::a,b,c
      real(8),dimension(nhyd),intent(out)::d
      integer:: n

      do n=1,nhyd
         d(n) = sign(1.0d0,a(n))*max(0.0d0,min(abs(a(n))         &
     &                                  ,sign(1.0d0,a(n))*b(n)   &
     &                                  ,sign(1.0d0,a(n))*c(n))) 
      enddo

      return
      end subroutine MClimiter

      subroutine NumericalFlux1
      use commons, only: is,ie,in,js,je,jn,ks,ke,kn,gam
      use fluxmod
      implicit none
      integer::i,j,k
      real(8),dimension(nhyd):: dsvp,dsvm,dsvc,dsv
      real(8),dimension(nhyd,in,jn,kn):: leftpr,rigtpr
      real(8),dimension(2*mflx+madd,in,jn,kn):: leftco,rigtco
      real(8),dimension(2*mflx+madd):: leftst,rigtst
      real(8),dimension(mflx):: nflux

!$acc declare create(leftco,rigtco)
!$acc declare create(leftpr,rigtpr)

!$acc data present(leftco,rigtco,leftpr,rigtpr)

!$acc kernels      
!$acc loop independent
      do k=ks,ke
      do j=js,je
!$acc loop independent private(dsv,dsvp,dsvm)
      do i=is-1,ie+1
         dsvp(:) = (svc(:,i+1,j,k) -svc(:,i,j,k)                 )
         dsvm(:) = (                svc(:,i,j,k) - svc(:,i-1,j,k))

         call vanLeer(dsvp,dsvm,dsv)
!         call minmod(dsvp,dsvm,dsv)
         leftpr(:,i+1,j,k) = svc(:,i,j,k) + 0.5d0*dsv(:)
         rigtpr(:,i  ,j,k) = svc(:,i,j,k) - 0.5d0*dsv(:)
      enddo
      enddo
      enddo
!$acc end kernels

!$acc kernels
!$acc loop independent
      do k=ks,ke
      do j=js,je
      do i=is,ie+1
         leftco(mudn,i,j,k)=leftpr(nden,i,j,k) ! rho
         leftco(muvu,i,j,k)=leftpr(nve1,i,j,k)*leftpr(nden,i,j,k)  ! rho v_x
         leftco(muvv,i,j,k)=leftpr(nve2,i,j,k)*leftpr(nden,i,j,k)  ! rho v_y
         leftco(muvw,i,j,k)=leftpr(nve3,i,j,k)*leftpr(nden,i,j,k)  ! rho v_z
         leftco(muet,i,j,k)=leftpr(nene,i,j,k)*leftpr(nden,i,j,k) &! e_i+ rho v^2/2
     &               +0.5d0*leftpr(nden,i,j,k)*(    &
     &                     +leftpr(nve1,i,j,k)**2   &
     &                     +leftpr(nve2,i,j,k)**2   &
     &                     +leftpr(nve3,i,j,k)**2)

         leftco(mfdn,i,j,k)=leftpr(nden,i,j,k)                   *leftpr(nve1,i,j,k)
         leftco(mfvu,i,j,k)=leftpr(nden,i,j,k)*leftpr(nve1,i,j,k)*leftpr(nve1,i,j,k) &
     &                     +leftpr(npre,i,j,k)
         leftco(mfvv,i,j,k)=leftpr(nden,i,j,k)*leftpr(nve2,i,j,k)*leftpr(nve1,i,j,k)
         leftco(mfvw,i,j,k)=leftpr(nden,i,j,k)*leftpr(nve3,i,j,k)*leftpr(nve1,i,j,k)
         leftco(mfet,i,j,k)=(leftpr(nene,i,j,k)*leftpr(nden,i,j,k)  &
     &               +0.5d0*leftpr(nden,i,j,k)*(   &
     &                     +leftpr(nve1,i,j,k)**2  &
     &                     +leftpr(nve2,i,j,k)**2  &
     &                     +leftpr(nve3,i,j,k)**2) &
     &                     +leftpr(npre,i,j,k)     &
     &                       )                                  *leftpr(nve1,i,j,k) 

         leftco(mcsp,i,j,k)= sqrt(gam*(gam-1.0d0)*leftpr(nene,i,j,k))
         leftco(mvel,i,j,k)= leftpr(nve1,i,j,k)
         leftco(mpre,i,j,k)= leftpr(npre,i,j,k)


         rigtco(mudn,i,j,k)=rigtpr(nden,i,j,k)
         rigtco(muvu,i,j,k)=rigtpr(nve1,i,j,k)*rigtpr(nden,i,j,k)
         rigtco(muvv,i,j,k)=rigtpr(nve2,i,j,k)*rigtpr(nden,i,j,k)
         rigtco(muvw,i,j,k)=rigtpr(nve3,i,j,k)*rigtpr(nden,i,j,k)
         rigtco(muet,i,j,k)=rigtpr(nene,i,j,k)*rigtpr(nden,i,j,k) &
     &               +0.5d0*rigtpr(nden,i,j,k)*(  &
     &                     +rigtpr(nve1,i,j,k)**2 &
     &                     +rigtpr(nve2,i,j,k)**2 &
     &                     +rigtpr(nve3,i,j,k)**2)

         rigtco(mfdn,i,j,k)=rigtpr(nden,i,j,k)                   *rigtpr(nve1,i,j,k)
         rigtco(mfvu,i,j,k)=rigtpr(nden,i,j,k)*rigtpr(nve1,i,j,k)*rigtpr(nve1,i,j,k) &
     &                     +rigtpr(npre,i,j,k)
         rigtco(mfvv,i,j,k)=rigtpr(nden,i,j,k)*rigtpr(nve2,i,j,k)*rigtpr(nve1,i,j,k)
         rigtco(mfvw,i,j,k)=rigtpr(nden,i,j,k)*rigtpr(nve3,i,j,k)*rigtpr(nve1,i,j,k)
         rigtco(mfet,i,j,k)=(rigtpr(nene,i,j,k)*rigtpr(nden,i,j,k) &
     &               +0.5d0*rigtpr(nden,i,j,k)*(   &
     &                     +rigtpr(nve1,i,j,k)**2  &
     &                     +rigtpr(nve2,i,j,k)**2  &
     &                     +rigtpr(nve3,i,j,k)**2) &
     &                     +rigtpr(npre,i,j,k)     &
     &                      )                                    *rigtpr(nve1,i,j,k)

         rigtco(mcsp,i,j,k)= sqrt(gam*(gam-1.0d0)*rigtpr(nene,i,j,k))
         rigtco(mvel,i,j,k)= rigtpr(nve1,i,j,k)
         rigtco(mpre,i,j,k)= rigtpr(npre,i,j,k)

      enddo
      enddo
      enddo
!$acc end kernels

!$acc kernels
!$acc loop independent
      do k=ks,ke
      do j=js,je
!$acc loop independent private(leftst,rigtst,nflux)
      do i=is,ie+1
         leftst(:)=leftco(:,i,j,k)
         rigtst(:)=rigtco(:,i,j,k)
!         call HLLE(leftst,rigtst,nflux)
         call HLLC(leftst,rigtst,nflux)
         nflux1(mden,i,j,k)=nflux(mden)
         nflux1(mrv1,i,j,k)=nflux(mrvu)
         nflux1(mrv2,i,j,k)=nflux(mrvv)
         nflux1(mrv3,i,j,k)=nflux(mrvw)
         nflux1(meto,i,j,k)=nflux(meto)
      enddo
      enddo
      enddo
!$acc end kernels
!$acc end data

      return
      end subroutine Numericalflux1

      subroutine NumericalFlux2
      use commons, only: is,ie,in,js,je,jn,ks,ke,kn,gam
      use fluxmod
      implicit none
      integer::i,j,k
      real(8),dimension(nhyd):: dsvp,dsvm,dsvc,dsv
      real(8),dimension(nhyd,in,jn,kn):: leftpr,rigtpr
      real(8),dimension(2*mflx+madd,in,jn,kn):: leftco,rigtco
      real(8),dimension(2*mflx+madd):: leftst,rigtst
      real(8),dimension(mflx):: nflux

!$acc declare create(leftco,rigtco)
!$acc declare create(leftpr,rigtpr)

!$acc data present(leftco,rigtco,leftpr,rigtpr)

!$acc kernels
!$acc loop independent
      do k=ks,ke
      do i=is,ie
!$acc loop independent private(dsv,dsvp,dsvm)
      do j=js-1,je+1
         dsvp(:) = (svc(:,i,j+1,k) -svc(:,i,j,k)                 )
         dsvm(:) = (                svc(:,i,j,k) - svc(:,i,j-1,k))

         call vanLeer(dsvp,dsvm,dsv)
!         call minmod(dsvp,dsvm,dsv)
         leftpr(:,i,j+1,k) = svc(:,i,j,k) + 0.5d0*dsv(:)
         rigtpr(:,i,j  ,k) = svc(:,i,j,k) - 0.5d0*dsv(:)

!         leftpr(:,i,j,k) = svc(:,i,j-1,k)
!         rigtpr(:,i,j,k) = svc(:,i,j  ,k)

       enddo
       enddo
       enddo
!$acc end kernels

!$acc kernels
!$acc loop independent
      do k=ks,ke
      do i=is,ie
      do j=js,je+1
         leftco(mudn,i,j,k)=leftpr(nden,i,j,k)
         leftco(muvw,i,j,k)=leftpr(nve1,i,j,k)*leftpr(nden,i,j,k)
         leftco(muvu,i,j,k)=leftpr(nve2,i,j,k)*leftpr(nden,i,j,k) ! rho v
         leftco(muvv,i,j,k)=leftpr(nve3,i,j,k)*leftpr(nden,i,j,k)
         leftco(muet,i,j,k)=leftpr(nene,i,j,k)*leftpr(nden,i,j,k) &
     &               +0.5d0*leftpr(nden,i,j,k)*(  &
     &                     +leftpr(nve1,i,j,k)**2 &
     &                     +leftpr(nve2,i,j,k)**2 &
     &                     +leftpr(nve3,i,j,k)**2)

         leftco(mfdn,i,j,k)=leftpr(nden,i,j,k)                   *leftpr(nve2,i,j,k) ! rho v
         leftco(mfvw,i,j,k)=leftpr(nden,i,j,k)*leftpr(nve1,i,j,k)*leftpr(nve2,i,j,k)
         leftco(mfvu,i,j,k)=leftpr(nden,i,j,k)*leftpr(nve2,i,j,k)*leftpr(nve2,i,j,k) &
     &                     +leftpr(npre,i,j,k)
         leftco(mfvv,i,j,k)=leftpr(nden,i,j,k)*leftpr(nve3,i,j,k)*leftpr(nve2,i,j,k)
         leftco(mfet,i,j,k)=(leftpr(nene,i,j,k)*leftpr(nden,i,j,k) &
     &               +0.5d0*leftpr(nden,i,j,k)*(   &
     &                     +leftpr(nve1,i,j,k)**2  &
     &                     +leftpr(nve2,i,j,k)**2  &
     &                     +leftpr(nve3,i,j,k)**2) &
     &                     +leftpr(npre,i,j,k)     &
     &                                       )*leftpr(nve2,i,j,k)

         leftco(mcsp,i,j,k)= sqrt(gam*(gam-1.0d0)*leftpr(nene,i,j,k))
         leftco(mvel,i,j,k)= leftpr(nve2,i,j,k)
         leftco(mpre,i,j,k)= leftpr(npre,i,j,k)


         rigtco(mudn,i,j,k)=rigtpr(nden,i,j,k)
         rigtco(muvw,i,j,k)=rigtpr(nve1,i,j,k)*rigtpr(nden,i,j,k)
         rigtco(muvu,i,j,k)=rigtpr(nve2,i,j,k)*rigtpr(nden,i,j,k)
         rigtco(muvv,i,j,k)=rigtpr(nve3,i,j,k)*rigtpr(nden,i,j,k)
         rigtco(muet,i,j,k)=rigtpr(nene,i,j,k)*rigtpr(nden,i,j,k) &
     &               +0.5d0*rigtpr(nden,i,j,k)*(   &
     &                     +rigtpr(nve1,i,j,k)**2  &
     &                     +rigtpr(nve2,i,j,k)**2  &
     &                     +rigtpr(nve3,i,j,k)**2)

         rigtco(mfdn,i,j,k)=rigtpr(nden,i,j,k)                   *rigtpr(nve2,i,j,k)
         rigtco(mfvw,i,j,k)=rigtpr(nden,i,j,k)*rigtpr(nve1,i,j,k)*rigtpr(nve2,i,j,k)
         rigtco(mfvu,i,j,k)=rigtpr(nden,i,j,k)*rigtpr(nve2,i,j,k)*rigtpr(nve2,i,j,k) &
     &                     +rigtpr(npre,i,j,k)
         rigtco(mfvv,i,j,k)=rigtpr(nden,i,j,k)*rigtpr(nve3,i,j,k)*rigtpr(nve2,i,j,k)
         rigtco(mfet,i,j,k)=(rigtpr(nene,i,j,k)*rigtpr(nden,i,j,k) &
     &               +0.5d0*rigtpr(nden,i,j,k)*(   &
     &                     +rigtpr(nve1,i,j,k)**2  &
     &                     +rigtpr(nve2,i,j,k)**2  &
     &                     +rigtpr(nve3,i,j,k)**2) &
     &                     +rigtpr(npre,i,j,k)     &
     &                                       )*rigtpr(nve2,i,j,k)

         rigtco(mcsp,i,j,k)= sqrt(gam*(gam-1.0d0)*rigtpr(nene,i,j,k))
         rigtco(mvel,i,j,k)= rigtpr(nve2,i,j,k)
         rigtco(mpre,i,j,k)= rigtpr(npre,i,j,k)

      enddo
      enddo
      enddo
!$acc end kernels

!$acc kernels
!$acc loop independent
      do k=ks,ke
      do i=is,ie
!$acc loop independent private(dsv,dsvp,dsvm)
      do j=js,je+1
         leftst(:)=leftco(:,i,j,k)
         rigtst(:)=rigtco(:,i,j,k)
!         call HLLE(leftst,rigtst,nflux)
         call HLLC(leftst,rigtst,nflux)
         nflux2(mden,i,j,k)=nflux(mden)
         nflux2(mrv1,i,j,k)=nflux(mrvw)
         nflux2(mrv2,i,j,k)=nflux(mrvu) ! mrv2=3, mrvu=2
         nflux2(mrv3,i,j,k)=nflux(mrvv)
         nflux2(meto,i,j,k)=nflux(meto)
      enddo
      enddo
      enddo
!$acc end kernels
!$acc end data

      return
      end subroutine Numericalflux2


      subroutine NumericalFlux3
      use commons, only: is,ie,in,js,je,jn,ks,ke,kn,gam
      use fluxmod
      implicit none
      integer::i,j,k
      real(8),dimension(nhyd):: dsvp,dsvm,dsvc,dsv
      real(8),dimension(nhyd,in,jn,kn):: leftpr,rigtpr
      real(8),dimension(2*mflx+madd,in,jn,kn):: leftco,rigtco
      real(8),dimension(2*mflx+madd):: leftst,rigtst
      real(8),dimension(mflx):: nflux

!$acc declare create(leftco,rigtco)
!$acc declare create(leftpr,rigtpr)

!$acc data present(leftco,rigtco,leftpr,rigtpr)

!$acc kernels
!$acc loop independent
      do j=js,je
      do i=is,ie
!$acc loop independent private(dsv,dsvp,dsvm)
      do k=ks-1,ke+1
         dsvp(:) = (svc(:,i,j,k+1) -svc(:,i,j,k)                 )
         dsvm(:) = (                svc(:,i,j,k) - svc(:,i,j,k-1))

         call vanLeer(dsvp,dsvm,dsv)
!         call minmod(dsvp,dsvm,dsv)
         leftpr(:,i,j,k+1) = svc(:,i,j,k) + 0.5d0*dsv(:)
         rigtpr(:,i,j,k  ) = svc(:,i,j,k) - 0.5d0*dsv(:)

       enddo
       enddo
       enddo
!$acc end kernels

!$acc kernels
!$acc loop independent
      do j=js,je
      do i=is,ie
      do k=ks,ke+1
         leftco(mudn,i,j,k)=leftpr(nden,i,j,k)
         leftco(muvv,i,j,k)=leftpr(nve1,i,j,k)*leftpr(nden,i,j,k)
         leftco(muvw,i,j,k)=leftpr(nve2,i,j,k)*leftpr(nden,i,j,k) ! rho v
         leftco(muvu,i,j,k)=leftpr(nve3,i,j,k)*leftpr(nden,i,j,k)
         leftco(muet,i,j,k)=leftpr(nene,i,j,k)*leftpr(nden,i,j,k) &
     &               +0.5d0*leftpr(nden,i,j,k)*(                  &
     &                     +leftpr(nve1,i,j,k)**2                 &
     &                     +leftpr(nve2,i,j,k)**2                 &
     &                     +leftpr(nve3,i,j,k)**2)

         leftco(mfdn,i,j,k)=leftpr(nden,i,j,k)                   *leftpr(nve3,i,j,k) ! rho v
         leftco(mfvv,i,j,k)=leftpr(nden,i,j,k)*leftpr(nve1,i,j,k)*leftpr(nve3,i,j,k)
         leftco(mfvw,i,j,k)=leftpr(nden,i,j,k)*leftpr(nve2,i,j,k)*leftpr(nve3,i,j,k)
         leftco(mfvu,i,j,k)=leftpr(nden,i,j,k)*leftpr(nve3,i,j,k)*leftpr(nve3,i,j,k) &
     &                     +leftpr(npre,i,j,k)
         leftco(mfet,i,j,k)=(leftpr(nene,i,j,k)*leftpr(nden,i,j,k) &
     &               +0.5d0*leftpr(nden,i,j,k)*(                   &
     &                     +leftpr(nve1,i,j,k)**2                  &
     &                     +leftpr(nve2,i,j,k)**2                  &
     &                     +leftpr(nve3,i,j,k)**2)                 &
     &                     +leftpr(npre,i,j,k)                     &
     &                                       )*leftpr(nve3,i,j,k)

         leftco(mcsp,i,j,k)= sqrt(gam*(gam-1.0d0)*leftpr(nene,i,j,k))
         leftco(mvel,i,j,k)= leftpr(nve3,i,j,k)
         leftco(mpre,i,j,k)= leftpr(npre,i,j,k)


         rigtco(mudn,i,j,k)=rigtpr(nden,i,j,k)
         rigtco(muvv,i,j,k)=rigtpr(nve1,i,j,k)*rigtpr(nden,i,j,k)
         rigtco(muvw,i,j,k)=rigtpr(nve2,i,j,k)*rigtpr(nden,i,j,k)
         rigtco(muvu,i,j,k)=rigtpr(nve3,i,j,k)*rigtpr(nden,i,j,k)
         rigtco(muet,i,j,k)=rigtpr(nene,i,j,k)*rigtpr(nden,i,j,k) &
     &               +0.5d0*rigtpr(nden,i,j,k)*(                  &
     &                     +rigtpr(nve1,i,j,k)**2                 &
     &                     +rigtpr(nve2,i,j,k)**2                 &
     &                     +rigtpr(nve3,i,j,k)**2)

         rigtco(mfdn,i,j,k)=rigtpr(nden,i,j,k)                   *rigtpr(nve3,i,j,k)
         rigtco(mfvv,i,j,k)=rigtpr(nden,i,j,k)*rigtpr(nve1,i,j,k)*rigtpr(nve3,i,j,k)
         rigtco(mfvw,i,j,k)=rigtpr(nden,i,j,k)*rigtpr(nve2,i,j,k)*rigtpr(nve3,i,j,k)
         rigtco(mfvu,i,j,k)=rigtpr(nden,i,j,k)*rigtpr(nve3,i,j,k)*rigtpr(nve3,i,j,k) &
     &                     +rigtpr(npre,i,j,k)
         rigtco(mfet,i,j,k)=(rigtpr(nene,i,j,k)*rigtpr(nden,i,j,k) &
     &               +0.5d0*rigtpr(nden,i,j,k)*(                   &
     &                     +rigtpr(nve1,i,j,k)**2                  &
     &                     +rigtpr(nve2,i,j,k)**2                  &
     &                     +rigtpr(nve3,i,j,k)**2)                 &
     &                     +rigtpr(npre,i,j,k)                     & 
     &                                       )*rigtpr(nve3,i,j,k)

         rigtco(mcsp,i,j,k)= sqrt(gam*(gam-1.0d0)*rigtpr(nene,i,j,k))
         rigtco(mvel,i,j,k)= rigtpr(nve3,i,j,k)
         rigtco(mpre,i,j,k)= rigtpr(npre,i,j,k)

      enddo
      enddo
      enddo
!$acc end kernels

!$acc kernels
!$acc loop independent
      do j=js,je
      do i=is,ie
!$acc loop independent private(leftst,rigtst,nflux)
      do k=ks,ke+1
         leftst(:)=leftco(:,i,j,k)
         rigtst(:)=rigtco(:,i,j,k)
!         call HLLE(leftst,rigtst,nflux)
         call HLLC(leftst,rigtst,nflux)
         nflux3(mden,i,j,k)=nflux(mden)
         nflux3(mrv1,i,j,k)=nflux(mrvv)
         nflux3(mrv2,i,j,k)=nflux(mrvw)
         nflux3(mrv3,i,j,k)=nflux(mrvu)
         nflux3(meto,i,j,k)=nflux(meto)
      enddo
      enddo
      enddo
!$acc end kernels
!$acc end data

      return
      end subroutine Numericalflux3

      subroutine HLLE(leftst,rigtst,nflux)
!$acc routine seq
      use fluxmod
      implicit none
      real(8),dimension(2*mflx+madd),intent(in)::leftst,rigtst
      real(8),dimension(mflx),intent(out)::nflux
      real(8),dimension(mflx)::ul,ur,fl,fr
      real(8)::csl,csr
      real(8):: vl, vr
      real(8):: sl, sr

      ul(1:mflx) = leftst(1:mflx)
      fl(1:mflx) = leftst(mflx+1:2*mflx)
      ur(1:mflx) = rigtst(1:mflx)
      fr(1:mflx) = rigtst(mflx+1:2*mflx)
      csl=leftst(mcsp)
      csr=rigtst(mcsp)
       vl=leftst(mvel)
       vr=rigtst(mvel)

       sl = min(vl,vr) - max(csl,csr)
       sl = min(0.0d0,sl)
       sr = max(vl,vr) + max(csl,csr)
       sr = max(0.0d0,sr)

       nflux(:) = (sr*fl(:)-sl*fr(:) +sl*sr*(ur(:)-ul(:)))/(sr-sl)

      return
      end subroutine HLLE

      subroutine HLLC(leftst,rigtst,nflux)
!$acc routine seq
!=====================================================================
!
! HLLC Scheme
!
! Purpose
! Calculation of Numerical Flux by HLLC method
!
! Reference
!  Toro EF, Spruce M, Speares W. (1992,1994)
!
! Input
! Output
!=====================================================================
      use fluxmod, only: mflx,madd                 &
     &                 , mudn,muvu,muvv,muvw,muet  &
     &                 , mfdn,mfvu,mfvv,mfvw,mfet  &
     &                 , mcsp,mvel,mpre            &
     &                 , mden,mrvu,mrvv,mrvw,meto

      implicit none
      real(8),dimension(2*mflx+madd),intent(in)::leftst,rigtst
      real(8),dimension(mflx),intent(out)::nflux

!----- U -----
! qql :: left state
! qqr :: right state
      real(8) :: rol,vxl,vyl,vzl,ptl,eel
      real(8) :: ror,vxr,vyr,vzr,ptr,eer
      real(8) :: rxl,ryl,rzl
      real(8) :: rxr,ryr,rzr
      real(8) :: ptst

!----- U* ----
! qqlst ::  left state
! qqrst :: right state
      real(8) :: rolst,vxlst,vylst,vzlst,eelst
      real(8) :: rorst,vxrst,vyrst,vzrst,eerst
      real(8) :: rxlst,rylst,rzlst
      real(8) :: rxrst,ryrst,rzrst

!----- flux ---
! fqql ::  left physical flux
! fqqr :: right physical flux
      real(8) :: frol,frxl,fryl,frzl,feel
      real(8) :: fror,frxr,fryr,frzr,feer

!----- wave speed ---
! sl ::  left-going fastest signal velocity
! sr :: right-going fastest signal velocity
! sm :: contact discontinuity velocity
! slst ::  left-going alfven velocity
! srst :: right-going alfven velocity
      real(8) :: sm,sl,sr

! cfl :: left-state Fast wave velocity
! cfr :: right-sate Fast wave velocity
      real(8) :: cfl,cfr

!--------------------
! temporary variables
      real(8) :: sdl,sdr,sdml,sdmr,isdml,isdmr,rosdl,rosdr
      real(8) :: temp
  
! no if
      real(8) :: sign1,maxs1,mins1
      real(8) :: msl,msr

!----- Step 0. ----------------------------------------------------------|

!---- Left state
        
        rol = leftst(mudn)
        eel = leftst(muet)
        rxl = leftst(muvu)
        ryl = leftst(muvv)
        rzl = leftst(muvw)
        vxl = leftst(muvu)/leftst(mudn)
        vyl = leftst(muvv)/leftst(mudn)
        vzl = leftst(muvw)/leftst(mudn)
        ptl = leftst(mpre)

!---- Right state
        
        ror = rigtst(mudn)
        eer = rigtst(muet)
        rxr = rigtst(muvu)
        ryr = rigtst(muvv)
        rzr = rigtst(muvw)
        vxr = rigtst(muvu)/rigtst(mudn)
        vyr = rigtst(muvv)/rigtst(mudn)
        vzr = rigtst(muvw)/rigtst(mudn)
        ptr = rigtst(mpre)
!----- Step 1. ----------------------------------------------------------|
! Compute wave left & right wave speed
!
         
        cfl = leftst(mcsp)
        cfr = rigtst(mcsp)

        sl = min(vxl,vxr)-max(cfl,cfr) ! note sl is negative
        sr = max(vxl,vxr)+max(cfl,cfr)
!----- Step 2. ----------------------------------------------------------|
! compute L/R fluxs
!
! Left value
        frol = leftst(mfdn)
        feel = leftst(mfet)
        frxl = leftst(mfvu)
        fryl = leftst(mfvv)
        frzl = leftst(mfvw)

! Right value
! Left value
        fror = rigtst(mfdn)
        feer = rigtst(mfet)
        frxr = rigtst(mfvu)
        fryr = rigtst(mfvv) 
        frzr = rigtst(mfvw)

!----- Step 4. ----------------------------------------------------------|
! compute middle and alfven wave
!
        sdl = sl - vxl
        sdr = sr - vxr
        rosdl = rol*sdl
        rosdr = ror*sdr

        temp = 1.0d0/(rosdr - rosdl)
! Eq. 45
        sm = (rosdr*vxr - rosdl*vxl - ptr + ptl)*temp
           
        sdml = sl - sm; isdml = 1.0d0/sdml
        sdmr = sr - sm; isdmr = 1.0d0/sdmr
        
!----- Step 5. ----------------------------------------------------------|
! compute intermediate states
!
! Eq. 49
        ptst = (rosdr*ptl-rosdl*ptr+rosdl*rosdr*(vxr-vxl))*temp

!----- Step 5A. ----------------------------------------------------------|
! compute Ul*
!

        rolst = rol*sdl   *isdml
        vxlst = sm
        rxlst = rolst*vxlst
           
        vylst = vyl
        rylst = rolst*vylst
        vzlst = vzl
        rzlst = rolst*vzlst

        eelst =(sdl*eel - ptl*vxl + ptst*sm  )*isdml

!----- Step 5B. ----------------------------------------------------------|
! compute Ur*
!

        rorst   = rosdr   *isdmr
        vxrst = sm
        rxrst = rorst*vxrst
        vyrst = vyr
        ryrst = rorst*vyrst
        vzrst = vzr
        rzrst = rorst*vzrst
           
        eerst = (sdr*eer - ptr*vxr  + ptst*sm  )*isdmr
              
!----- Step 6. ----------------------------------------------------------|
! compute flux
        sign1 = sign(1.0d0,sm)    ! 1 for sm>0, -1 for sm<0
        maxs1 =  max(0.0d0,sign1) ! 1 sm>0, 0 for sm<0
        mins1 = -min(0.0d0,sign1) ! 0 sm>0,-1 for sm<0

        msl   = min(sl  ,0.0d0)   ! 0 for sl > 0, sl for sl < 0
        msr   = max(sr  ,0.0d0)   ! S_R > 0

        nflux(mden) = (frol+msl*(rolst-rol))*maxs1 &
     &               +(fror+msr*(rorst-ror))*mins1
        nflux(meto) = (feel+msl*(eelst-eel))*maxs1 &
     &               +(feer+msr*(eerst-eer))*mins1
        nflux(mrvu) = (frxl+msl*(rxlst-rxl))*maxs1 &
     &               +(frxr+msr*(rxrst-rxr))*mins1
        nflux(mrvv) = (fryl+msl*(rylst-ryl))*maxs1 &
     &               +(fryr+msr*(ryrst-ryr))*mins1
        nflux(mrvw) = (frzl+msl*(rzlst-rzl))*maxs1 &
     &               +(frzr+msr*(rzrst-rzr))*mins1

      return
      end subroutine HLLC

      subroutine UpdateConsv
      use commons
      use fluxmod
      implicit none
      integer::i,j,k

!$acc kernels
!$acc loop independent
      do k=ks,ke
      do j=js,je
      do i=is,ie
         
         d(i,j,k) = d(i,j,k)                       &
     & +dt*(                                       &
     &  (- nflux1(mden,i+1,j,k)                    &
     &   + nflux1(mden,i  ,j,k))/(x1a(i+1)-x1a(i)) &
     & +(- nflux2(mden,i,j+1,k)                    &
     &   + nflux2(mden,i,j  ,k))/(x2a(j+1)-x2a(j)) &
     & +(- nflux3(mden,i,j,k+1)                    &
     &   + nflux3(mden,i,j,k  ))/(x3a(k+1)-x3a(k)) &
     &      )

         mv1(i,j,k) = mv1(i,j,k)                   &
     & +dt*(                                       &
     &  (- nflux1(mrv1,i+1,j,k)                    &
     &   + nflux1(mrv1,i  ,j,k))/(x1a(i+1)-x1a(i)) &
     & +(- nflux2(mrv1,i,j+1,k)                    &
     &   + nflux2(mrv1,i,j  ,k))/(x2a(j+1)-x2a(j)) &
     & +(- nflux3(mrv1,i,j,k+1)                    &
     &   + nflux3(mrv1,i,j,k  ))/(x3a(k+1)-x3a(k)) &
     &      )

         mv2(i,j,k) = mv2(i,j,k)                   &    
     & +dt*(                                       &
     &  (- nflux1(mrv2,i+1,j,k)                    &
     &   + nflux1(mrv2,i  ,j,k))/(x1a(i+1)-x1a(i)) &
     & +(- nflux2(mrv2,i,j+1,k)                    &
     &   + nflux2(mrv2,i,j  ,k))/(x2a(j+1)-x2a(j)) &
     & +(- nflux3(mrv2,i,j,k+1)                    &
     &   + nflux3(mrv2,i,j,k  ))/(x3a(k+1)-x3a(k)) &
     &      )

         mv3(i,j,k) = mv3(i,j,k)                   & 
     & +dt*(                                       &
     &  (- nflux1(mrv3,i+1,j,k)                    &
     &   + nflux1(mrv3,i  ,j,k))/(x1a(i+1)-x1a(i)) &
     & +(- nflux2(mrv3,i,j+1,k)                    &
     &   + nflux2(mrv3,i,j  ,k))/(x2a(j+1)-x2a(j)) &
     & +(- nflux3(mrv3,i,j,k+1)                    &
     &   + nflux3(mrv3,i,j,k  ))/(x3a(k+1)-x3a(k)) &
     &      )

          et(i,j,k) = et(i,j,k)                    &
     & +dt*(                                       &
     &  (- nflux1(meto,i+1,j,k)                    &
     &   + nflux1(meto,i  ,j,k))/(x1a(i+1)-x1a(i)) &
     & +(- nflux2(meto,i,j+1,k)                    &
     &   + nflux2(meto,i,j  ,k))/(x2a(j+1)-x2a(j)) &
     & +(- nflux3(meto,i,j,k+1)                    &
     &   + nflux3(meto,i,j,k  ))/(x3a(k+1)-x3a(k)) &
     &      )
      enddo
      enddo
      enddo
!$acc end kernels

      return
      end subroutine UpdateConsv

      subroutine Output
      use commons
      implicit none
      integer::i,j,k
      character(20),parameter::dirname="bindata/"
      character(40)::filename
      real(8),save::tout
      data tout / 0.0d0 / 
      integer::nout
      data nout / 1 /
      integer,parameter::unitout=17
      integer,parameter::unitbin=13
      integer,parameter:: gs=1
      integer,parameter:: nvar=5
      real(8)::x1out(is-gs:ie+gs,2)
      real(8)::x2out(js-gs:je+gs,2)
      real(8)::x3out(js-gs:je+gs,2)
      real(8)::hydout(is-gs:ie+gs,js-gs:je+gs,ks-gs:ke+gs,nvar)

      logical, save:: is_inited
      data is_inited /.false./

      if (.not. is_inited) then
         call makedirs("bindata")
         is_inited =.true.
      endif


      if(time .lt. tout+dtout) return
!$acc update host (d,v1,v2,v3,p)

      x1out(is-gs:ie+gs,1) = x1b(is-gs:ie+gs)
      x1out(is-gs:ie+gs,2) = x1a(is-gs:ie+gs)

      x2out(is-gs:ie+gs,1) = x2b(is-gs:ie+gs)
      x2out(is-gs:ie+gs,2) = x2a(is-gs:ie+gs)

      x3out(ks-gs:ke+gs,1) = x3b(ks-gs:ke+gs)
      x3out(ks-gs:ke+gs,2) = x3a(ks-gs:ke+gs)

      hydout(is-gs:ie+gs,js-gs:je+gs,ks-gs:ke+gs,1) =  d(is-gs:ie+gs,js-gs:je+gs,ks-gs:ke+gs)
      hydout(is-gs:ie+gs,js-gs:je+gs,ks-gs:ke+gs,2) = v1(is-gs:ie+gs,js-gs:je+gs,ks-gs:ke+gs)
      hydout(is-gs:ie+gs,js-gs:je+gs,ks-gs:ke+gs,3) = v2(is-gs:ie+gs,js-gs:je+gs,ks-gs:ke+gs)
      hydout(is-gs:ie+gs,js-gs:je+gs,ks-gs:ke+gs,4) = v3(is-gs:ie+gs,js-gs:je+gs,ks-gs:ke+gs)
      hydout(is-gs:ie+gs,js-gs:je+gs,ks-gs:ke+gs,5) =  p(is-gs:ie+gs,js-gs:je+gs,ks-gs:ke+gs)


      write(filename,'(a3,i5.5,a4)')"unf",nout,".dat"
      filename = trim(dirname)//filename

      open(unitout,file=filename,status='replace',form='formatted')
      write(unitout,*) "# ",time,dt
      write(unitout,*) "# ",ngrid,gs
      write(unitout,*) "# ",ngrid,gs
      write(unitout,*) "# ",ngrid,gs
      close(unitout)

      write(filename,'(a3,i5.5,a4)')"bin",nout,".dat"
      filename = trim(dirname)//filename
      open(unitbin,file=filename,status='replace',form='binary') 
      write(unitbin) x1out(:,:)
      write(unitbin) x2out(:,:)
      write(unitbin) x3out(:,:)
      write(unitbin) hydout(:,:,:,:)
      close(unitbin)
      
      write(6,*) "output:",nout,time

      nout=nout+1
      tout=time

      return
      end subroutine Output


      subroutine makedirs(outdir)
      character(len=*), intent(in) :: outdir
      character(len=256) command
      write(command, *) 'if [ ! -d ', trim(outdir), ' ]; then mkdir -p ', trim(outdir), '; fi'
      write(*, *) trim(command)
      call system(command)
      end subroutine makedirs
