module fieldmod
    implicit none
    integer:: incr
    real(8):: time,dt
    integer:: in,jn,kn
    integer:: izone,jzone,kzone
    integer:: igs,jgs,kgs
    integer:: is,js,ks,ie,je,ke
    real(8),dimension(:),allocatable:: x1b,x2b,x3b
    real(8),dimension(:),allocatable:: x1a,x2a,x3a
    real(8),dimension(:,:,:),allocatable:: d,v1,v2,v3,p
    real(8),dimension(:,:,:),allocatable:: vor1,vor2,vor3
    real(8),dimension(:,:,:),allocatable:: kin,hk
    real(8):: dx,dy,dz

!$acc declare create(incr)
!$acc declare create(time,dt)
!$acc declare create(in,jn,kn)
!$acc declare create(izone,jzone,kzone)
!$acc declare create(igs,jgs,kgs)
!$acc declare create(is,js,ks)
!$acc declare create(ie,je,ke)
      
!$acc declare create(x1a,x1b)
!$acc declare create(x2a,x2b)
!$acc declare create(x3a,x3b)
      
!$acc declare create(d,v1,v2,v3,p)
!$acc declare create(vor1,vor2,vor3)
!$acc declare create(kin,hk)

!$acc declare create(dx,dy,dz)

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

  FILENUMBER: do incr  = fbeg,fend,10
     write(6,*) "file number",incr
     call ReadData
     call Vorticity
     call Fourier
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
  read(unitinp,*) dummy,kzone,kgs
  close(unitinp)
  in= izone+2*igs
  jn= jzone+2*jgs
  kn= kzone+2*kgs

  is=1+igs
  js=1+jgs
  ks=1+kgs
  ie=in-igs
  je=jn-jgs
  ke=kn-kgs

  if(.not. is_inited)then
     allocate( x1b(in),x1a(in))
     allocate( x2b(jn),x2a(jn))
     allocate( x3b(kn),x3a(kn))
     allocate( d(in,jn,kn))
     allocate(v1(in,jn,kn))
     allocate(v2(in,jn,kn))
     allocate(v3(in,jn,kn))
     allocate( p(in,jn,kn))
     is_inited = .true.
  endif

  write(filename,'(a3,i5.5,a4)')"bin",incr,".dat"
  filename = trim(dirname)//filename
  open(unitbin,file=filename,status='old',form='binary')
  read(unitbin)x1b(:),x1a(:)
  read(unitbin)x2b(:),x2a(:)
  read(unitbin)x3b(:),x3a(:)
  read(unitbin)  d(:,:,:)
  read(unitbin) v1(:,:,:)
  read(unitbin) v2(:,:,:)
  read(unitbin) v3(:,:,:)
  read(unitbin)  p(:,:,:)
  close(unitbin)
  
  dx = x1b(2)-x1b(1)
  dy = x2b(2)-x2b(1)
  dz = x3b(2)-x3b(1)

!$acc update device (x1a,x1b)
!$acc update device (x2a,x2b)
!$acc update device (x3a,x3b)
!$acc update device (d,v1,v2,v3,p)
!$acc update device (dx,dy,dz)

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
     allocate( vor1(in,jn,kn))
     allocate( vor2(in,jn,kn))
     allocate( vor3(in,jn,kn))
     allocate(  kin(in,jn,kn))
     allocate(   hk(in,jn,kn))
     is_inited = .true.
  endif

!$acc kernels
!$acc loop collapse(3) independent
  do k=ks,ke
  do j=js,je
  do i=is,ie
     vor1(i,j,k)= (v2(i  ,j,k)-v2(i-1,j,k))/dx*0.5 &
                &+(v2(i+1,j,k)-v2(i  ,j,k))/dx*0.5 &
                &-(v3(i,j  ,k)-v3(i,j-1,k))/dy*0.5 &
                &-(v3(i,j+1,k)-v3(i,j  ,k))/dy*0.5 
     vor2(i,j,k)= (v3(i  ,j,k)-v3(i-1,j,k))/dx*0.5 &
                &+(v3(i+1,j,k)-v3(i  ,j,k))/dx*0.5 &
                &-(v1(i,j  ,k)-v1(i,j-1,k))/dy*0.5 &
                &-(v1(i,j+1,k)-v1(i,j  ,k))/dy*0.5 
     vor3(i,j,k)= (v2(i  ,j,k)-v2(i-1,j,k))/dx*0.5 &
                &+(v2(i+1,j,k)-v2(i  ,j,k))/dx*0.5 &
                &-(v1(i,j  ,k)-v1(i,j-1,k))/dy*0.5 &
                &-(v1(i,j+1,k)-v1(i,j  ,k))/dy*0.5 
      kin(i,j,k)= 0.5d0*d(i,j,k)*( v1(i,j,k)*v1(i,j,k) &
                &                 +v2(i,j,k)*v2(i,j,k) &
                &                 +v3(i,j,k)*v3(i,j,k))

       hk(i,j,k)=  v1(i,j,k)*vor1(i,j,k) &
                & +v2(i,j,k)*vor2(i,j,k) &
                & +v3(i,j,k)*vor3(i,j,k)
  enddo
  enddo
  enddo
!$acc end kernels

!!$acc update host (vor1,vor2,vor3)
!!$acc update host (kin,hk)

!  write(6,*) "debug 1",kin(is,js,ks),hk(is,js,ks)
  
  return
end subroutine Vorticity
  
subroutine Fourier
  use fieldmod
  implicit none
  integer::i,j,k,n
  integer::ik,jk,kk,rk
  integer,parameter:: nk=128
  integer,parameter:: nvar=2
  real(8),dimension(:,:,:,:),allocatable:: X3D
  real(8)                :: Xtotloc
  real(8),dimension(nvar):: Xtot
  real(8),dimension(nk,nk,nk,nvar):: Xhat3DC,Xhat3DS
  real(8)                         :: Xhat3DCloc,Xhat3DSloc
  real(8),dimension(nk):: kx,ky,kz
  real(8),dimension(nk,nvar):: Xhat1D
  real(8):: kr
  real(8):: dkx,dky,dkz,dkr
  character(20),parameter::dirname="output/"
  character(40)::filename
  integer,parameter::unitspc=21
  integer,parameter::unittot=22
  real(8),save:: pi
  logical,save:: is_inited
  data is_inited / .false. /
!$acc declare create(X3D)
!$acc declare create(Xtot)
!$acc declare create(dkx,dky,dkz)
!$acc declare create(kx,ky,kz)
!$acc declare create(Xhat3DC,Xhat3DS,Xhat1D)
!$acc declare create(pi)

  if(.not. is_inited)then
     allocate( X3D(is:ie,js:je,ks:ke,nvar))
     pi=acos(-1.0d0)
!$acc update device (X3D)
!$acc update device (pi)
     is_inited = .true.
  endif

!$acc kernels
!$acc loop collapse(3) independent
  do k=ks,ke
  do j=js,je
  do i=is,ie
     X3D(i,j,k,1) = kin(i,j,k)
     X3D(i,j,k,2) =  hk(i,j,k)
  enddo
  enddo
  enddo
!$acc end kernels
  
!!$acc update host (X3D) 
! write(6,*) "debug 2",X3D(is,js,ks,1),X3D(is,js,ks,2)
  

!$acc kernels
!$acc loop independent private(Xtotloc)
  do n=1,nvar
     Xtotloc = 0.0d0
!$acc loop collapse(3)reduction(+:Xtotloc)
  do k=ks,ke
  do j=js,je
  do i=is,ie
     Xtotloc = Xtotloc + X3D(i,j,k,n)*dx*dy*dz
  enddo
  enddo
  enddo
     Xtot(n) = Xtotloc
  enddo
!$acc end kernels
!$acc update host (Xtot)
  
!  write(6,*) "debug 3",Xtot(1),Xtot(2)
  
  dkx = 1.0d0/(dx*(in-2*igs))
  dky = 1.0d0/(dy*(jn-2*jgs))
  dkz = 1.0d0/(dz*(kn-2*kgs))
!$acc update device (dkx,dky,dkz)
  
  do ik=1,nk
     kx(ik) = ik *dkx
  enddo
  do jk=1,nk
     ky(jk) = jk *dky
  enddo
  do kk=1,nk
     ky(kk) = kk *dkz
  enddo
!$acc update device (kx,ky,kz)

!$acc kernels
!$acc loop collapse(4) independent private(Xhat3DCloc)
  nloop: do n=1,nvar
  do kk=1,nk
  do jk=1,nk
  do ik=1,nk 
     Xhat3DCloc = 0.0d0    
!$acc loop collapse(3) reduction(+:Xhat3DCloc)
  do k=ks,ke
  do j=js,je
  do i=is,ie
     Xhat3DCloc = Xhat3DCloc &
 &    + X3D(i,j,k,n) &
 &    * cos(2.0d0*pi*(kx(ik)*x1b(i)+ky(jk)*x2b(j)+kz(kk)*x3b(k) )) & 
 &    *dx*dy*dz

  enddo
  enddo
  enddo
     Xhat3DC(ik,jk,kk,n) = Xhat3DCloc
  enddo
  enddo
  enddo

  enddo nloop
!$acc end kernels
!$acc update host (Xhat3DC)

!$acc kernels
!$acc loop collapse(4) independent private(Xhat3DSloc)
  do n=1,nvar  
  do kk=1,nk
  do jk=1,nk
  do ik=1,nk     
     Xhat3DSloc = 0.0d0
!$acc loop collapse(3) reduction(+:Xhat3DSloc)
  do k=ks,ke
  do j=js,je
  do i=is,ie
     Xhat3DSloc = Xhat3DSloc &
 &    + X3D(i,j,k,n) &
 &    * sin(2.0d0*pi*(kx(ik)*x1b(i)+ky(jk)*x2b(j)+kz(kk)*x3b(k) )) & 
 &    *dx*dy*dz

  enddo
  enddo
  enddo
     Xhat3DS(ik,jk,kk,n) = Xhat3DSloc 
  enddo
  enddo
  enddo
  enddo
!$acc end kernels
!$acc update host (Xhat3DS)

  Xhat1D(:,:) = 0.0d0
  dkr = dkx/sqrt(3.0d0) ! minimum k
  do kk=1,nk
  do jk=1,nk
  do ik=1,nk
     kr = sqrt(kx(ik)**2 +ky(jk)**2 +kz(kk)**2)
     rk = min(nk,int(kr/dkr))
     Xhat1D(rk,1:nvar) = Xhat1D(rk,1:nvar) + sqrt(   Xhat3DS(ik,jk,kk,1:nvar)**2 & 
                                           &       + Xhat3DC(ik,jk,kk,1:nvar)**2 &
                                           &      )*dkx*dky*dkz
  enddo
  enddo
  enddo

  write(filename,'(a3,i5.5,a4)')"spc",incr,".dat"
  filename = trim(dirname)//filename
  open(unitspc,file=filename,status='replace',form='formatted')
  write(unitspc,*) "# ",time
  do rk=1,nk
     write(unitspc,'(3(1x,E12.3))') rk*dkr,Xhat1D(rk,1)/Xtot(1),Xhat1D(rk,2)/Xtot(2)
  enddo
  close(unitspc)

  write(filename,'(a3,i5.5,a4)')"tot",incr,".dat"
  filename = trim(dirname)//filename
  open(unittot,file=filename,status='replace',form='formatted')
  write(unittot,'(6(1x,E12.3))') time,Xtot(1),Xtot(2)
  close(unittot)

  return
end subroutine Fourier
