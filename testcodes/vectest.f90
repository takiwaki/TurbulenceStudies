      module commons
      implicit none

      integer,parameter::ngrid=128
      integer,parameter::mgn=2
      integer,parameter::in=ngrid+2*mgn+1 &
     &                  ,jn=ngrid+2*mgn+1 &
     &                  ,kn=1
      integer,parameter::is=mgn+1 &
     &                  ,js=mgn+1 &
     &                  ,ks=1 
      integer,parameter::ie=ngrid+mgn &
     &                  ,je=ngrid+mgn &
     &                  ,ke=1

      real(8),parameter:: x1min=-0.5d0,x1max=0.5d0
      real(8),parameter:: x2min=-0.5d0,x2max=0.5d0
      real(8),dimension(in)::x1a,x1b
      real(8),dimension(jn)::x2a,x2b
      real(8),dimension(kn)::x3a,x3b

      real(8)::dx,dy

      real(8),dimension(in,jn,kn)::b1,b2,b3
      real(8),dimension(in,jn,kn)::psi,jcd,mpt,mptpo

      end module commons

      program vectest
      use commons
      implicit none
      logical::is_final
      data is_final /.false./

      write(6,*) "setup grids and fiels"
      call GenerateGrid
      call GenerateProblem
      call CurrentDensity
      call Potential
      call Snap2D


      write(6,*) "program has been finished"
      end program vectest

      subroutine GenerateGrid
      use commons
      implicit none
      integer::i,j,k
! x coordinates
      dx=(x1max-x1min)/dble(ngrid)
      do i=1,in
         x1a(i) = dx*(i-(mgn+1))+x1min
      enddo
      do i=1,in-1
         x1b(i) = 0.5d0*(x1a(i+1)+x1a(i))
      enddo
 
! y coordinates
      dy=(x2max-x2min)/dble(ngrid)
      do j=1,jn
         x2a(j) = dy*(j-(mgn+1))+x2min
      enddo

      do j=1,jn-1
         x2b(j) = 0.5d0*(x2a(j+1)+x2a(j))
      enddo

      return
      end subroutine GenerateGrid

      subroutine GenerateProblem
      use commons
      implicit none
      integer::i,j,k
      real(8)::pi
      real(8),parameter::k_ini=2.0d0

      real(8),dimension(in,jn,kn)::vpsi1b,vpsi2b
      real(8),dimension(in,jn,kn)::mpsi1b,mpsi2b
      real(8):: psinorm

      integer::seedsize
      integer,allocatable:: seed(:)
      real(8)::x

      real(8),parameter:: b0 = 1.0d0

      pi=acos(-1.0d0)
      psinorm = 1.0d0/(2.0d0*pi*k_ini)


      do k=ks,ke
      do j=js,je+1
      do i=is,ie+1
            psi(i,j,k) = psinorm * sin(k_ini*x1b(i)*2.0d0*pi/(x1max-x1min)) &
     &                           * sin(k_ini*x2b(j)*2.0d0*pi/(x2max-x2min))
         mpsi1b(i,j,k) = psinorm * sin(k_ini*x1b(i)*2.0d0*pi/(x1max-x1min)) &
     &                           * sin(k_ini*x2a(j)*2.0d0*pi/(x2max-x2min))
         mpsi2b(i,j,k) = psinorm * sin(k_ini*x1a(i)*2.0d0*pi/(x1max-x1min)) &
     &                           * sin(k_ini*x2b(j)*2.0d0*pi/(x2max-x2min))
      enddo
      enddo
      enddo

      do k=ks,ke
      do j=js,je
      do i=is,ie
         b1(i,j,k) =  b0*(mpsi1b(i,j+1,k)-mpsi1b(i,j,k))/(x2a(j+1)-x2a(j))
         b2(i,j,k) = -b0*(mpsi2b(i+1,j,k)-mpsi2b(i,j,k))/(x1a(i+1)-x1a(i))
      enddo
      enddo
      enddo

      call BoundaryCondition

      return
      end subroutine GenerateProblem

      subroutine BoundaryCondition
      use commons
      implicit none
      integer::i,j,k

      k=ks
      do j=1,jn-1
      do i=1,mgn
          b1(i,j,k) = b1(ie-mgn+i,j,k)
          b2(i,j,k) = b2(ie-mgn+i,j,k)
          b3(i,j,k) = b3(ie-mgn+i,j,k)
      enddo
      enddo

      k=ks
      do j=1,jn-1
      do i=1,mgn
          b1(ie+i,j,k) = b1(is+i-1,j,k)
          b2(ie+i,j,k) = b2(is+i-1,j,k)
          b3(ie+i,j,k) = b3(is+i-1,j,k)
      enddo
      enddo

      k=ks
      do i=1,in-1
      do j=1,mgn
          b1(i,j,k) = b1(i,je-mgn+j,k)
          b2(i,j,k) = b2(i,je-mgn+j,k)
          b3(i,j,k) = b3(i,je-mgn+j,k)
      enddo
      enddo

      k=ks
      do i=1,in-1
      do j=1,mgn
          b1(i,je+j,k) = b1(i,js+j-1,k)
          b2(i,je+j,k) = b2(i,js+j-1,k)
          b3(i,je+j,k) = b3(i,js+j-1,k)
      enddo
      enddo

      return
      end subroutine BoundaryCondition


subroutine CurrentDensity
  use commons
  implicit none
  integer::i,j,k

  character(20),parameter::dirname="output/"
  character(40)::filename
  integer,parameter::unitvor=123

  logical,save:: is_inited
  data is_inited / .false. /

  jcd(:,:,:) = 0.0d0
  k=ks
  do j=js,je
  do i=is,ie
     jcd(i,j,k)= (b2(i  ,j,k)-b2(i-1,j,k))/dx*0.5 &
               &+(b2(i+1,j,k)-b2(i  ,j,k))/dx*0.5 &
               &-(b1(i,j  ,k)-b1(i,j-1,k))/dy*0.5 &
               &-(b1(i,j+1,k)-b1(i,j  ,k))/dy*0.5 
  enddo
  enddo

  return
end subroutine CurrentDensity

module potmod
  use commons
  implicit none
  real(8),dimension(in,jn,kn):: p,rp,r,Ap,phi
  real(8)::a1,a2,a3,a4
  real(8)::tmp1,tmp2,tmp3
  real(8)::alpha,beta
end module potmod

subroutine Potential
  use commons
  use potmod
  implicit none
  integer::i,j,k
  integer::iter
  integer,parameter::itermax=1000
  real(8),parameter::eps=1.0d-8

  character(20),parameter::dirname="output/"
  character(40)::filename
  integer,parameter::unitvor=231

  logical,save:: is_inited
  data is_inited / .false. /

  mpt(:,:,:)= 0.0d0

  k=1
  do j=js,je
  do i=is+1,ie
     mpt(i,j,k) = mpt(i-1,j,k)  - (b2(i,j,k) + b2(i-1,j,k))/2.0d0*dx
  enddo
  enddo

  do j=js+1,je
  do i=is  ,ie
     mpt(i,j,k) = mpt(i,j-1,k) + (b1(i,j,k) + b1(i,j-1,k))/2.0d0*dy
  enddo
  enddo

  a1=1.0d0/dx**2
  a2=1.0d0/dy**2
  a4=2.0d0*(a1+a2)
  
  do j=js,je; do i=is,ie
     r(i,j,k) = -jcd(i,j,k)
  enddo; enddo

  do j=js,je; do i=is,ie
       p(i,j,k) = r(i,j,k)
     phi(i,j,k) = 0.0d0
  enddo; enddo

  do j=js,je
     p(is-1,j,k) = p(ie,j,k)
     p(ie+1,j,k) = p(is,j,k)
  enddo

  do i=is,ie
     p(i,js-1,k) = p(i,je,k)
     p(i,je+1,k) = p(i,js,k)    
  enddo

  itloop: do iter=1,itermax

  do j=js,je; do i=is,ie
     Ap(i,j,k) =  a1*(p(i+1,j,k)+p(i-1,j,k)) &
               & +a2*(p(i,j+1,k)+p(i,j-1,k)) &
               & -a4*(p(i,j,k))
  enddo; enddo

  tmp1=0.0d0
  do j=js,je; do i=is,ie
     tmp1= tmp1 +  r(i,j,k)**2
  enddo; enddo
  tmp2=0.0d0
  do j=js,je; do i=is,ie
     tmp2= tmp2 +  Ap(i,j,k)*p(i,j,k)
  enddo; enddo

  alpha= tmp1/tmp2
  do j=js,je; do i=is,ie
     phi(i,j,k) = phi(i,j,k) +alpha * p(i,j,k)
     rp(i,j,k) =    r(i,j,k) -alpha *Ap(i,j,k)
  enddo; enddo

  tmp3 =0.0d0
  do j=js,je; do i=is,ie
     tmp3= tmp3 + rp(i,j,k)**2      
  enddo; enddo
  if(tmp3 .lt. eps) then
     write(6,*) iter,eps
     exit itloop
  endif
  beta =  tmp3/tmp1
  do j=js,je; do i=is,ie
     p(i,j,k) = rp(i,j,k) +beta * p(i,j,k)
     r(i,j,k) = rp(i,j,k)
  enddo; enddo

  do j=js,je
     p(is-1,j,k) = p(ie,j,k)
     p(ie+1,j,k) = p(is,j,k)
  enddo

  do i=is,ie
     p(i,js-1,k) = p(i,je,k)
     p(i,je+1,k) = p(i,js,k)    
  enddo
  enddo itloop

  do j=js,je; do i=is,ie
     mptpo(i,j,k) = phi(i,j,k)
  enddo; enddo

  return
end subroutine Potential


subroutine Snap2D
  use commons
  implicit none
  integer::i,j,k

  character(20),parameter::dirname="./"
  character(40)::filename
  integer,parameter::unitvor=123

  logical,save:: is_inited
  data is_inited / .false. /
  
  write(filename,*)"test.dat"
  filename = filename
  write(6,*) filename
  open(unitvor,file=filename,status='replace',form='formatted')
  k=1
  write(unitvor,*) "# x y omega_z"
  do j=js,je
  do i=is,ie
     write(unitvor,'(10(1x,E12.3))') x1b(i),x2b(j) &
                                 & ,psi(i,j,k),jcd(i,j,k),mpt(i,j,k),mptpo(i,j,k)

  enddo
     write(unitvor,*)
  enddo

  close(unitvor)


  return
end subroutine Snap2D
 
