module fieldmod
    implicit none
    integer:: incr
    real(8):: time,dt
    integer:: in,jn,kn
    integer:: izone,jzone,kzone
    integer:: igs,jgs,kgs
    integer:: is,js,ks,ie,je,ke
    real(8),dimension(:),allocatable:: x1b,x2b
    real(8),dimension(:),allocatable:: x1a,x2a
    real(8),dimension(:,:,:),allocatable:: d,v1,v2,v3,p
    real(8),dimension(:,:,:),allocatable:: b1,b2,b3,bp
    real(8),dimension(:,:,:),allocatable:: vor, kin ! vorticity
    real(8),dimension(:,:,:),allocatable:: jcd, mag ! current density
    real(8),dimension(:,:,:),allocatable:: mpt ! magnetic potential
    real(8),dimension(:,:,:),allocatable:: mptpo ! magnetic potential
    real(8),dimension(:,:,:),allocatable:: hcr ! cross helicity
    real(8):: dx,dy
end module fieldmod

program data_analysis
  use fieldmod
  implicit none
  integer:: fbeg, fend
  logical:: flag
  integer,parameter:: unitcon=100

  INQUIRE(FILE ="control.dat",EXIST = flag)
  if(flag) then
     open (unitcon,file="control.dat" &
     &        ,status='old',form='formatted')
     read (unitcon,*) fbeg,fend
     close(unitcon)
  endif

  FILENUMBER: do incr  = fbeg,fend
     write(6,*) "file number",incr
     call ReadData
     call Vorticity
     call Potential
     call Snap2D
     call Fourier
     call Probability
  enddo FILENUMBER

end program data_analysis

subroutine ReadData
  use fieldmod
  implicit none   
  character(20),parameter::dirname="bindata/"
  character(40)::filename
  integer,parameter::unitinp=13
  integer,parameter::unitbin=14
  character(8)::dummy
  logical,save:: is_inited
  data is_inited / .false. /

  write(filename,'(a3,i5.5,a4)')"unf",incr,".dat"
  filename = trim(dirname)//filename
  open(unitinp,file=filename,status='old',form='formatted')
  read(unitinp,*) dummy,time,dt
  read(unitinp,*) dummy,izone,igs
  read(unitinp,*) dummy,jzone,jgs
  close(unitinp)
  in=izone+2*igs
  jn=jzone+2*jgs
  kn=1
  is=1+igs
  js=1+jgs
  ks=1
  ie=in-igs
  je=jn-jgs
  ke=1

  if(.not. is_inited)then
     allocate( x1b(in),x1a(in))
     allocate( x2b(jn),x2a(jn))
     allocate( d(in,jn,kn))
     allocate(v1(in,jn,kn))
     allocate(v2(in,jn,kn))
     allocate(v3(in,jn,kn))
     allocate(b1(in,jn,kn))
     allocate(b2(in,jn,kn))
     allocate(b3(in,jn,kn))
     allocate(bp(in,jn,kn))
     allocate( p(in,jn,kn))
     is_inited = .true.
  endif

  write(filename,'(a3,i5.5,a4)')"bin",incr,".dat"
  filename = trim(dirname)//filename
  open(unitbin,file=filename,status='old',form='binary')
  read(unitbin)x1b(:),x1a(:)
  read(unitbin)x2b(:),x2a(:)
  read(unitbin)  d(:,:,:)
  read(unitbin) v1(:,:,:)
  read(unitbin) v2(:,:,:)
  read(unitbin) v3(:,:,:)
  read(unitbin) b1(:,:,:)
  read(unitbin) b2(:,:,:)
  read(unitbin) b3(:,:,:)
  read(unitbin) bp(:,:,:)
  read(unitbin)  p(:,:,:)
  close(unitbin)
  
  dx = x1b(2)-x1b(1)
  dy = x2b(2)-x2b(1)

  return
end subroutine ReadData

subroutine Vorticity
  use fieldmod
  implicit none
  integer::i,j,k

  character(20),parameter::dirname="output/"
  character(40)::filename
  integer,parameter::unitvor=123

  logical,save:: is_inited
  data is_inited / .false. /

  if(.not. is_inited)then
     allocate( vor(in,jn,kn))
     allocate( jcd(in,jn,kn))
     allocate( kin(in,jn,kn))
     allocate( mag(in,jn,kn))
     allocate( Hcr(in,jn,kn))
     is_inited = .true.
  endif

  k=ks
  do j=js,je
  do i=is,ie
     vor(i,j,k)= (v2(i  ,j,k)-v2(i-1,j,k))/dx*0.5 &
               &+(v2(i+1,j,k)-v2(i  ,j,k))/dx*0.5 &
               &-(v1(i,j  ,k)-v1(i,j-1,k))/dy*0.5 &
               &-(v1(i,j+1,k)-v1(i,j  ,k))/dy*0.5
     jcd(i,j,k)= (b2(i  ,j,k)-b2(i-1,j,k))/dx*0.5 &
               &+(b2(i+1,j,k)-b2(i  ,j,k))/dx*0.5 &
               &-(b1(i,j  ,k)-b1(i,j-1,k))/dy*0.5 &
               &-(b1(i,j+1,k)-b1(i,j  ,k))/dy*0.5 
     kin(i,j,k)= 0.5d0*d(i,j,k)*( &
               & +v1(i,j,k)*v1(i,j,k) &
               & +v2(i,j,k)*v2(i,j,k) &
               & )
     mag(i,j,k)= ( &
               & +b1(i,j,k)*b1(i,j,k) &
               & +b2(i,j,k)*b2(i,j,k) &
               & )
     Hcr(i,j,k) =     v1(i,j,k) * b1(i,j,k) &
                &  +  v2(i,j,k) * b2(i,j,k) &
                &  +  v3(i,j,k) * b3(i,j,k)
  enddo
  enddo

  return
end subroutine Vorticity

module potmod
  implicit none
  real(8),dimension(:,:,:),allocatable:: p,rp,r,Ap,phi
  real(8)::a1,a2,a3,a4
  real(8)::tmp1,tmp2,tmp3
  real(8)::alpha,beta
end module potmod

subroutine Potential
  use fieldmod, press=>p
  use potmod
  implicit none
  integer::i,j,k
  integer::iter
  integer,parameter::itermax=100
  real(8),parameter::eps=3.0d-2

  character(20),parameter::dirname="output/"
  character(40)::filename
  integer,parameter::unitvor=231

  logical,save:: is_inited
  data is_inited / .false. /


  if(.not. is_inited)then
     allocate(  p(in,jn,kn))
     allocate( rp(in,jn,kn))
     allocate(  r(in,jn,kn))
     allocate( Ap(in,jn,kn))
     allocate(phi(in,jn,kn))

     allocate( mpt(in,jn,kn))
     allocate( mptpo(in,jn,kn))
     is_inited = .true.
  endif
 

  mpt(:,:,:)= 0.0d0

  k=ks
  do j=js+1,je
  do i=is  ,ie
!     mpt(i,j,k) = mpt(i,j-1,k) + (b1(i,j,k) + b1(i,j-1,k))/2.0d0*dy
  enddo
  enddo

  do j=js,je
  do i=is+1,ie
     write(6,*)(b2(i,j,k) + b2(i-1,j,k))/2.0d0*dx
     mpt(i,j,k) = mpt(i-1,j,k) - (b2(i,j,k) + b2(i-1,j,k))/2.0d0*dx
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
!  write(6,*) tmp3
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


subroutine Fourier
  use fieldmod
  implicit none
  integer::i,j,k
  integer::ik,jk,kk,rk
  integer,parameter:: nk=128
  integer,parameter:: nvar=5
  real(8),dimension(nvar):: X
  real(8),dimension(nvar):: Xtot
  real(8),dimension(nk,nk,nvar):: Xhat2Dc,Xhat2Ds
  real(8),dimension(nk):: kx,ky
  real(8),dimension(nk,nvar):: Xhat1D
  real(8):: kr
  real(8):: dkx,dky,dkr
  character(20),parameter::dirname="output/"
  character(40)::filename
  integer,parameter::unitspc=21
  integer,parameter::unittot=22
  real(8):: pi

  pi=acos(-1.0d0)

  k=1
  
  Xtot(:)=0.0d0
  do j=js,je
  do i=is,ie
     Xtot(1) = Xtot(1) + kin(i,j,k)    *dx*dy
     Xtot(2) = Xtot(2) + vor(i,j,k)**2 *dx*dy
     Xtot(3) = Xtot(3) + mpt(i,j,k)**2 *dx*dy
     Xtot(4) = Xtot(4) + Hcr(i,j,k)    *dx*dy
     Xtot(5) = Xtot(5) + mag(i,j,k)    *dx*dy
  enddo
  enddo

  do j=js,je
  do i=is,ie


  enddo
  enddo

  dkx = 1.0d0/(dx*in)
  dky = 1.0d0/(dy*jn)
  
  do ik=1,nk
     kx(ik) = ik *dkx
  enddo
  do jk=1,nk
     ky(jk) = jk *dky
  enddo

  Xhat2Dc(:,:,:) = 0.0d0
  Xhat2Ds(:,:,:) = 0.0d0

  do ik=1,nk
  do jk=1,nk

  do j=js,je
  do i=is,ie
     X(1) = kin(i,j,k) 
     X(2) = vor(i,j,k)**2
     X(3) = mpt(i,j,k)**2
     X(4) = hcr(i,j,k)
     X(5) = mag(i,j,k)

     Xhat2Dc(ik,jk,1:nvar) = Xhat2Dc(ik,jk,1:nvar)  &
 &    + X(1:nvar)                                   &
 &    * cos(2.0d0*pi*(kx(ik)*x1b(i)+ky(jk)*x2b(j))) & 
 &    *dx*dy
     Xhat2Ds(ik,jk,1:nvar) = Xhat2Ds(ik,jk,1:nvar)  &
 &    + X(1:nvar)                                   &
 &    * sin(2.0d0*pi*(kx(ik)*x1b(i)+ky(jk)*x2b(j))) & 
 &    *dx*dy

  enddo
  enddo

  enddo
  enddo

  Xhat1D(:,1:nvar) = 0.0d0
  dkr = dkx/sqrt(2.0d0) ! minimum k
  do ik=1,nk
  do jk=1,nk
     kr = sqrt(kx(ik)**2+ky(jk)**2)
     rk = min(nk,int(kr/dkr))
     Xhat1D(rk,1:nvar) = Xhat1D(rk,1:nvar) + sqrt(Xhat2Dc(ik,jk,1:nvar)**2 + Xhat2Ds(ik,jk,1:nvar)**2)*dkx*dky/dkr
  enddo
  enddo

  write(filename,'(a3,i5.5,a4)')"spc",incr,".dat"
  filename = trim(dirname)//filename
  open(unitspc,file=filename,status='replace',form='formatted')
  write(unitspc,'(1a,1(1x,E12.3))') "#",time
!                                               12345678   12345678   12345678   12345678   12345678
  write(unitspc,'(1a,6(1x,a8))') "#","1:k    ","2:E_kin ","3:ens   ","4:mpot  ","5:cross ","6:E_mag "
  do rk=1,nk
     write(unitspc,'(1x,6(1x,E12.3))') rk*dkr,Xhat1D(rk,1)/Xtot(1) & ! kinetic energy
                                     &       ,Xhat1D(rk,2)/Xtot(2) & ! enstrophy
                                     &       ,Xhat1D(rk,3)/Xtot(3) & ! magpot
                                     &       ,Xhat1D(rk,4)         & ! cross helicity
                                     &       ,Xhat1D(rk,5)/Xtot(5)   ! magnetic energy

  enddo
  close(unitspc)

  write(filename,'(a3,i5.5,a4)')"tot",incr,".dat"
  filename = trim(dirname)//filename
  open(unittot,file=filename,status='replace',form='formatted')
  write(unittot,'(1x,6(1x,E12.3))') time,Xtot(1),Xtot(2),Xtot(3),Xtot(4),Xtot(5)
  close(unittot)

  return
end subroutine Fourier
  
subroutine Probability
  use fieldmod
  implicit none
  integer::i,j,k,n
  integer,parameter:: np=100
  real(8),dimension(-np:np):: vxpro,vypro
  character(20),parameter::dirname="output/"
  character(40)::filename
  integer,parameter::unitpro=331
  real(8):: vxmax, vymax

  k=1
  vxmax=0.0d0
  vymax=0.0d0
  do j=js,je
  do i=is,ie
     vxmax = max(vxmax, abs(v1(i,j,k)))
     vymax = max(vymax, abs(v2(i,j,k)))
  enddo
  enddo

  vxpro(:)= 0.0d0
  vypro(:)= 0.0d0
  do j=js,je
  do i=is,ie 
     n= min(np,max(-np,int(v1(i,j,k)*np/vxmax)))
     vxpro(n) = vxpro(n) + dx*dy
     n= min(np,max(-np,int(v2(i,j,k)*np/vymax)))
     vypro(n) = vypro(n) + dx*dy
  enddo
  enddo

  write(filename,'(a3,i5.5,a4)')"pro",incr,".dat"
  filename = trim(dirname)//filename
  open(unitpro,file=filename,status='replace',form='formatted')
  write(unitpro,*) "# ",time
  do n=-np,np
     write(unitpro,'(4(1x,E12.3))') vxmax/np*n, vxpro(n), vymax/np*n, vypro(n)
  enddo
  close(unitpro)

  return
end subroutine Probability

subroutine Snap2D
  use fieldmod
  implicit none
  integer::i,j,k

  character(20),parameter::dirname="output/"
  character(40)::filename
  integer,parameter::unitvor=123

  logical,save:: is_inited
  data is_inited / .false. /

  write(filename,'(a3,i5.5,a4)')"vor",incr,".dat"
  filename = trim(dirname)//filename
  open(unitvor,file=filename,status='replace',form='formatted')
  k=1
  write(unitvor,*) "# ",time
  write(unitvor,*) "# x y omega_z"
  do j=js,je
  do i=is,ie
     write(unitvor,'(10(1x,E12.3))') x1b(i),x2b(j) &
                                 & ,vor(i,j,k),jcd(i,j,k),kin(i,j,k),mag(i,j,k) &
                                 & ,mpt(i,j,k),mptpo(i,j,k)

  enddo
     write(unitvor,*)
  enddo

  close(unitvor)


  return
end subroutine Snap2D
 
