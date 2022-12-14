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

!$acc declare create(incr)
!$acc declare create(time,dt)
!$acc declare create(in,jn,kn)
!$acc declare create(izone,jzone,kzone)
!$acc declare create(igs,jgs,kgs)
!$acc declare create(is,js,ks)
!$acc declare create(ie,je,ke)
      
!$acc declare create(x1a,x1b)
!$acc declare create(x2a,x2b)
      
!$acc declare create(d,v1,v2,v3,p)
!$acc declare create(b1,b2,b3,bp)
!$acc declare create(vor,kin)
!$acc declare create(jcd,mag)
!$acc declare create(mpt,hcr)
!$acc declare create(mptpo)

!$acc declare create(dx,dy)
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

  write(6,*) "data index span from",fbeg," - ",fend

  FILENUMBER: do incr  = fbeg,fend
     write(6,*) "file number",incr
     call ReadData
     call Vorticity
     call Potential
     call Snap2D
!     call Fourier
  enddo FILENUMBER

  stop
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
  write(6,*) "data size",in,jn,kn

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
  write(6,*) "open ",filename
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
  
  write(6,*) "data read end"
  
  dx = x1b(2)-x1b(1)
  dy = x2b(2)-x2b(1)

!$acc update device (in,jn,kn)
!$acc update device (is,js,ks)
!$acc update device (ie,je,ke)
!$acc update device (x1a,x1b)
!$acc update device (x2a,x2b)
!$acc update device (d,v1,v2,v3,p)
!$acc update device (b1,b2,b3,bp)
!$acc update device (dx,dy)
  
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
!$acc update device (vor)
!$acc update device (jcd)
!$acc update device (kin)
!$acc update device (mag)
!$acc update device (Hcr)
     is_inited = .true.
  endif

!$acc kernels
  k=1
!$acc loop collapse(2) independent
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

!$acc end kernels
!$acc update host (vor)
!$acc update host (jcd)
!$acc update host (kin)
!$acc update host (mag)
!$acc update host (Hcr)

  return
end subroutine Vorticity

module potmod
  implicit none
  real(8),dimension(:,:,:),allocatable:: p,rp,r,Ap,phi
  real(8)::a1,a2,a3,a4
  real(8)::tmp1,tmp2,tmp3
  real(8)::alpha,beta
!$acc declare create(p,rp,r,Ap,phi)
!$acc declare create(a1,a2,a3,a4)
!$acc declare create(tmp1,tmp2,tmp3)
!$acc declare create(alpha,beta)

  real(8)::sumtest 
!$acc declare create(sumtest)

end module potmod

subroutine Potential
  use fieldmod, press=>p
  use potmod
  implicit none
  integer::i,j,k
  integer::iter
  integer,parameter::itermax=1000
  real(8),parameter::eps=3.0d-2

  character(20),parameter::dirname="output/"
  character(40)::filename
  integer,parameter::unitvor=231

  logical,save:: is_inited
  data is_inited / .false. /

  write(6,*) "calculate vector potential"

  if(.not. is_inited)then
     allocate(  p(in,jn,kn));   p(:,:,:)=0.0d0
     allocate( rp(in,jn,kn));  rp(:,:,:)=0.0d0
     allocate(  r(in,jn,kn));   r(:,:,:)=0.0d0
     allocate( Ap(in,jn,kn));  Ap(:,:,:)=0.0d0
     allocate(phi(in,jn,kn)); phi(:,:,:)=0.0d0
!$acc update device (p,rp,r,Ap,phi)

     allocate( mpt(in,jn,kn))
     allocate( mptpo(in,jn,kn))
!$acc update device (mpt,mptpo)
     is_inited = .true.
  endif

!$acc kernels
  k=ks
!$acc loop collapse(2) independent
  do j=js,je
  do i=is,ie
  mpt(i,j,k)= 0.0d0
  enddo
  enddo
!$acc end kernels

!$acc kernels
  k=ks
!$acc loop collapse(2) independent
  do j=js+1,je
  do i=is  ,ie
     mpt(i,j,k) = mpt(i,j-1,k) + (b1(i,j,k) + b1(i,j-1,k))/2.0d0*dy
  enddo
  enddo
!$acc end kernels

!$acc kernels
  k=ks
!$acc loop collapse(2) independent
  do j=js,je
  do i=is+1,ie
     mpt(i,j,k) = mpt(i-1,j,k) - (b2(i,j,k) + b2(i-1,j,k))/2.0d0*dx
  enddo
  enddo
!$acc end kernels

  a1=1.0d0/dx**2
  a2=1.0d0/dy**2
  a4=2.0d0*(a1+a2)
!$acc update device (a1,a2,a4)

!$acc kernels
  k=ks
!$acc loop collapse(2) independent
  do j=js,je; do i=is,ie
     r(i,j,k) = -jcd(i,j,k)
  enddo; enddo
!$acc end kernels

!$acc kernels
  k=ks
!$acc loop collapse(2) independent
  do j=js,je; do i=is,ie
       p(i,j,k) = r(i,j,k)
     phi(i,j,k) = 0.0d0
  enddo; enddo
!$acc end kernels

!$acc kernels
  k=ks
!$acc loop independent
  do j=js,je
     p(is-1,j,k) = p(ie,j,k)
     p(ie+1,j,k) = p(is,j,k)
  enddo

!$acc loop independent
  do i=is,ie
     p(i,js-1,k) = p(i,je,k)
     p(i,je+1,k) = p(i,js,k)    
  enddo
!$acc end kernels

!$acc loop seq
  itloop: do iter=1,itermax
     if(iter .eq. itermax) then
        write(6,*) "iteration reaches max"
     endif

!$acc kernels
  k=ks           
!$acc loop collapse(2) independent
  do j=js,je; do i=is,ie
     Ap(i,j,k) =  a1*(p(i+1,j,k)+p(i-1,j,k)) &
               & +a2*(p(i,j+1,k)+p(i,j-1,k)) &
               & -a4*(p(i,j,k))
  enddo; enddo
!$acc end kernels
!$acc update host (a1,a2,a4)

  tmp1=0.0d0
!$acc update device (tmp1)
!$acc kernels
  k=ks
!$acc loop collapse(2) reduction(+:tmp1)
  do j=js,je; do i=is,ie
     tmp1 = tmp1 + r(i,j,k)**2
  enddo; enddo
!$acc end kernels
!$acc update host (tmp1)
  
  tmp2=0.0d0
!$acc update device (tmp2)
!$acc kernels
  k=ks
!$acc loop collapse(2) reduction(+:tmp2)
  do j=js,je; do i=is,ie
     tmp2 = tmp2 +  Ap(i,j,k)*p(i,j,k)
  enddo; enddo
!$acc end kernels
!$acc update host (tmp2)

  alpha= tmp1/tmp2
!$acc update device (alpha)

!$acc kernels
  k=ks
!$acc loop collapse(2) independent
  do j=js,je; do i=is,ie
     phi(i,j,k) = phi(i,j,k) + alpha * p(i,j,k)
      rp(i,j,k) =   r(i,j,k) - alpha *Ap(i,j,k)
  enddo; enddo
!$acc end kernels

  tmp3=0.0d0
!$acc update device (tmp3)
!$acc kernels
  k=ks
!$acc loop collapse(2) reduction(+:tmp3)
  do j=js,je; do i=is,ie
     tmp3= tmp3 + rp(i,j,k)**2      
  enddo; enddo
!$acc end kernels
!$acc update host (tmp3)
  
  if(tmp3 .lt. eps) then
     write(6,*) "iteration end at step=",iter,"eps=",eps
     exit itloop
  endif
  
  beta =  tmp3/tmp1
!$acc update device (beta)
  
!$acc kernels
  k=ks 
!$acc loop collapse(2) independent
  do j=js,je; do i=is,ie
     p(i,j,k) = rp(i,j,k) + beta * p(i,j,k)
     r(i,j,k) = rp(i,j,k)
  enddo; enddo
!$acc end kernels

!$acc kernels
  k=ks
!$acc loop independent
  do j=js,je
     p(is-1,j,k) = p(ie,j,k)
     p(ie+1,j,k) = p(is,j,k)
  enddo
!$acc end kernels

!$acc kernels
  k=ks
!$acc loop independent
  do i=is,ie
     p(i,js-1,k) = p(i,je,k)
     p(i,je+1,k) = p(i,js,k)    
  enddo
!$acc end kernels
  
!  write(6,*) "debug2",iter, tmp1, tmp2, tmp3, alpha,beta

  enddo itloop

!$acc kernels
  k=ks
!$acc loop collapse(2) independent
  do j=js,je; do i=is,ie
     mptpo(i,j,k) = phi(i,j,k)
  enddo; enddo
!$acc end kernels
!$acc update host (mpt,phi,mptpo)
  write(6,*) "debug", mpt(is+1,js+1,ks),mptpo(is+1,js+1,ks),phi(is+1,js+1,ks)
  return
end subroutine Potential

module spctrmod
implicit none
  real(8),dimension(:,:,:),allocatable:: X2D
!$acc declare create(X2D)
  integer,parameter:: nk=128
  integer,parameter:: nvar=5
  real(8),dimension(nvar):: Xtot
  real(8),dimension(nk,nk,nvar):: Xhat2DC,Xhat2DS
  real(8),dimension(nk):: kx,ky,kz
  real(8),dimension(nk,nvar):: Xhat1D
  real(8):: kr
  real(8):: dkx,dky,dkz,dkr
  
!$acc declare create(Xtot)
!$acc declare create(dkx,dky,dkz)
!$acc declare create(kx,ky,kz)
!$acc declare create(Xhat2DC,Xhat2DS,Xhat1D)
  
  real(8) :: pi
!$acc declare create(pi)
end module spctrmod

subroutine Fourier
  use fieldmod
  use spctrmod
  implicit none
  integer::i,j,k,n
  integer::ik,jk,kk,rk
  real(8):: Xtotloc
  real(8):: Xhat2DCloc,Xhat2DSloc
  character(20),parameter::dirname="output/"
  character(40)::filename
  integer,parameter::unitspc=21
  integer,parameter::unittot=22
  logical,save:: is_inited
  data is_inited / .false. /

  if(.not. is_inited)then
     allocate( X2D(is:ie,js:je,nvar))
     pi=acos(-1.0d0)
!$acc update device (X2D)
!$acc update device (pi)
     is_inited = .true.
  endif
!$acc kernels
  k=ks
!$acc loop collapse(2) independent
  do j=js,je
  do i=is,ie
     X2D(i,j,1) = kin(i,j,k)
     X2D(i,j,2) = vor(i,j,k)**2
     X2D(i,j,3) = mpt(i,j,k)**2
     X2D(i,j,4) = Hcr(i,j,k)**2
     X2D(i,j,5) = mag(i,j,k)**2
  enddo
  enddo
!$acc end kernels

  
!$acc kernels
!$acc loop independent private(Xtotloc)
  do n=1,nvar
     Xtotloc = 0.0d0
  k=ks
!$acc loop collapse(2) reduction(+:Xtotloc)
  do j=js,je
  do i=is,ie
     Xtotloc = Xtotloc + X2D(i,j,n) * dx*dy
  enddo
  enddo
     Xtot(n) = Xtotloc
  enddo
!$acc end kernels
!$acc update host (Xtot)

  dkx = 1.0d0/(dx*(in-2*igs))
  dky = 1.0d0/(dy*(jn-2*jgs))
!$acc update device (dkx,dky)
  
  do ik=1,nk
     kx(ik) = ik *dkx
  enddo
  do jk=1,nk
     ky(jk) = jk *dky
  enddo
!$acc update device (kx,ky)


!$acc kernels
!$acc loop collapse(3) independent private(Xhat2DCloc,Xhat2DSloc)
  nloop: do n=1,nvar
  do jk=1,nk
  do ik=1,nk 
     Xhat2DCloc = 0.0d0    
!$acc loop collapse(2) reduction(+:Xhat2DCloc)
  do j=js,je
  do i=is,ie
     Xhat2DCloc = Xhat2DCloc &
 &    + X2D(i,j,n) &
 &    * cos(2.0d0*pi*(kx(ik)*x1b(i)+ky(jk)*x2b(j))) & 
 &    * dx*dy

  enddo
  enddo
     Xhat2DC(ik,jk,n) = Xhat2DCloc
  
     Xhat2DSloc = 0.0d0
!$acc loop collapse(2) reduction(+:Xhat2DSloc)
  do j=js,je
  do i=is,ie
     Xhat2DSloc = Xhat2DSloc &
 &    + X2D(i,j,n) &
 &    * sin(2.0d0*pi*(kx(ik)*x1b(i)+ky(jk)*x2b(j))) & 
 &    * dx*dy

  enddo
  enddo
     Xhat2DS(ik,jk,n) = Xhat2DSloc 
  enddo
  enddo

  enddo nloop
!$acc end kernels
!$acc update host (Xhat2DC)
!$acc update host (Xhat2DS)

  

  Xhat1D(:,1:nvar) = 0.0d0
  dkr = dkx/sqrt(2.0d0) ! minimum k
  do ik=1,nk
  do jk=1,nk
     kr = sqrt(kx(ik)**2+ky(jk)**2)
     rk = min(nk,int(kr/dkr))
     Xhat1D(rk,1:nvar) = Xhat1D(rk,1:nvar) + sqrt(Xhat2DC(ik,jk,1:nvar)**2 + Xhat2DS(ik,jk,1:nvar)**2)*dkx*dky/dkr
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

subroutine Snap2D
  use fieldmod
  implicit none
  integer::i,j,k

  character(20),parameter::dirname="output/"
  character(40)::filename
  integer,parameter::unitvor=123

  logical,save:: is_inited
  data is_inited / .false. /

!$acc update host (vor,jcd,kin,mag,mpt,mptpo)
    
  write(filename,'(a3,i5.5,a4)')"vor",incr,".dat"
  filename = trim(dirname)//filename
  open(unitvor,file=filename,status='replace',form='formatted')
  k=ks
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
 
