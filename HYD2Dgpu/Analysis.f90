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
    real(8),dimension(:,:,:),allocatable:: vor
    real(8),dimension(:,:,:),allocatable:: kin
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
!$acc declare create(vor)
!$acc declare create(kin)

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

  FILENUMBER: do incr  = fbeg,fend
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
  close(unitinp)
  in=izone+2*igs
  jn=jzone+2*jgs
  kn=1

  is=1+igs
  js=1+jgs
  ie=in-igs
  je=jn-jgs

  if(.not. is_inited)then
     allocate( x1b(in),x1a(in))
     allocate( x2b(jn),x2a(jn))
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
  read(unitbin)  d(:,:,:)
  read(unitbin) v1(:,:,:)
  read(unitbin) v2(:,:,:)
  read(unitbin) v3(:,:,:)
  read(unitbin)  p(:,:,:)
  close(unitbin)
  
  dx = x1b(2)-x1b(1)
  dy = x2b(2)-x2b(1)

!$acc update device (x1a,x1b)
!$acc update device (x2a,x2b)
!$acc update device (d,v1,v2,v3,p)
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
     allocate( kin(in,jn,kn))
!!$acc update device (vor1)
!!$acc update device (kin)
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
     kin(i,j,k)= 0.5d0*d(i,j,k)*( &
               & +v1(i,j,k)*v1(i,j,k) &
               & +v2(i,j,k)*v2(i,j,k) &
               & )
  enddo
  enddo
!$acc end kernels
!$acc update host (vor)
!$acc update host (kin)

  write(filename,'(a3,i5.5,a4)')"vor",incr,".dat"
  filename = trim(dirname)//filename
  open(unitvor,file=filename,status='replace',form='formatted')

  write(unitvor,'(1a,4(1x,E12.3))') "#",time
  write(unitvor,'(1a,4(1x,a8))') "#","1:x    ","2:y     ","3:omg_z ","4:E_kin "
  do j=js,je
  do i=is,ie
     write(unitvor,'(4(1x,E12.3))') x1b(i),x2b(j),vor(i,j,k),kin(i,j,k)
  enddo
     write(unitvor,*)
  enddo

  close(unitvor)


  return
end subroutine Vorticity
  
module spctrmod
implicit none
  real(8),dimension(:,:,:),allocatable:: X2D
!$acc declare create(X2D)
  integer,parameter:: nk=128
  integer,parameter:: nvar=2
  real(8),dimension(nvar):: Xtot
  real(8),dimension(nk,nk,nvar):: Xhat2DC,Xhat2DS
  real(8),dimension(nk):: kx,ky,kz
  real(8),dimension(nk,nvar):: Xhat1D
  real(8):: kr
  real(8):: dkx,dky,dkz,dkr
  
  real(8) :: pi
!$acc declare create(Xtot)
!$acc declare create(dkx,dky,dkz)
!$acc declare create(kx,ky,kz)
!$acc declare create(Xhat2DC,Xhat2DS,Xhat1D)
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
  enddo
  enddo
!$acc end kernels

!$acc kernels
!$acc loop independent private(Xtotloc)
  do n=1,nvar
     Xtotloc = 0.0d0
  k=ks
!$acc loop collapse(2)reduction(+:Xtotloc)
  do j=js,je
  do i=is,ie
     Xtotloc = Xtotloc + X2D(i,j,n)*dx*dy
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
 &    *dx*dy

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
 &    *dx*dy

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
  do rk=1,nk
     write(unitspc,'(1x,6(1x,E12.3))') rk*dkr,Xhat1D(rk,1)/Xtot(1) & ! kinetic energy
                                           & ,Xhat1D(rk,2)/Xtot(2)   ! enstrophy
  enddo
  close(unitspc)

  write(filename,'(a3,i5.5,a4)')"tot",incr,".dat"
  filename = trim(dirname)//filename
  open(unittot,file=filename,status='replace',form='formatted')
  write(unittot,'(6(1x,E12.3))') time,Xtot(1),Xtot(2)
  close(unittot)

  return
end subroutine Fourier
