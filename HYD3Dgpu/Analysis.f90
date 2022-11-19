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
    real(8):: dx,dy,dz
    real(8):: Etot,Vtot
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
     call Probability
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
     is_inited = .true.
  endif

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
  enddo
  enddo
  enddo

  write(filename,'(a3,i5.5,a4)')"vor",incr,".dat"
  filename = trim(dirname)//filename
  open(unitvor,file=filename,status='replace',form='formatted')

  k=1
  write(unitvor,*) "# ",time
  write(unitvor,*) "# x y omega_z"
  do j=js,je
  do i=is,ie
     write(unitvor,'(3(1x,E12.3))') x1b(i),x2b(j),vor3(i,j,k)
  enddo
     write(unitvor,*)
  enddo

  close(unitvor)


  return
end subroutine Vorticity
  
subroutine Fourier
  use fieldmod
  implicit none
  integer::i,j,k
  integer::ik,jk,kk,rk
  integer,parameter:: nk=100
  real(8),dimension(nk,nk,nk):: Ehat3Dc,Vhat3Dc,Ehat3Ds,Vhat3Ds
  real(8),dimension(nk):: kx,ky,kz
  real(8),dimension(nk):: Ehat1D,Vhat1D
  real(8):: kr
  real(8):: dkx,dky,dkz,dkr
  character(20),parameter::dirname="output/"
  character(40)::filename
  integer,parameter::unitspc=21
  real(8):: pi

  pi=acos(-1.0d0)

  Etot=0.0d0
  Vtot=0.0d0
  do k=ks,ke
  do j=js,je
  do i=is,ie
     Etot = Etot &
 &    + 0.5d0*d(i,j,k)                              &
 &    *(v1(i,j,k)*v1(i,j,k) + v2(i,j,k)*v2(i,j,k)+ v3(i,j,k)*v3(i,j,k))  & 
 &    *dx*dy*dz

     Vtot = Vtot &
 &    + vor3(i,j,k)**2                               & 
 &    *dx*dy*dz

  enddo
  enddo
  enddo

  dkx = 1.0d0/(dx*in)
  dky = 1.0d0/(dy*jn)
  dkz = 1.0d0/(dz*jn)
  
  do ik=1,nk
     kx(ik) = ik *dkx
  enddo
  do jk=1,nk
     ky(jk) = jk *dky
  enddo
  do kk=1,nk
     ky(kk) = kk *dkz
  enddo

  Ehat3Dc(:,:,:) = 0.0d0
  Ehat3Ds(:,:,:) = 0.0d0
  Vhat3Dc(:,:,:) = 0.0d0
  Vhat3Ds(:,:,:) = 0.0d0

  do kk=1,nk
  do jk=1,nk
  do ik=1,nk

  do k=ks,ke
  do j=js,je
  do i=is,ie
     Ehat3Dc(ik,jk,kk) = Ehat3Dc(ik,jk,kk) &
 &    + 0.5d0*d(i,j,k)                              &
 &    *(v1(i,j,k)*v1(i,j,k) + v2(i,j,k)*v2(i,j,k) + v3(i,j,k)*v3(i,j,k))  &
 &    * cos(2.0d0*pi*(kx(ik)*x1b(i)+ky(jk)*x2b(j))) & 
 &    *dx*dy*dz
     Ehat3Ds(ik,jk,kk) = Ehat3Ds(ik,jk,kk) &
 &    + 0.5d0*d(i,j,k)                              &
 &    *(v1(i,j,k)*v1(i,j,k) + v2(i,j,k)*v2(i,j,k) + v3(i,j,k)*v3(i,j,k))  &
 &    * sin(2.0d0*pi*(kx(ik)*x1b(i)+ky(jk)*x2b(j))) & 
 &    *dx*dy*dz

     Vhat3Dc(ik,jk,kk) = Vhat3Dc(ik,jk,kk) &
 &    + vor3(i,j,k)**2                               &
 &    * cos(2.0d0*pi*(kx(ik)*x1b(i)+ky(jk)*x2b(j))) & 
 &    *dx*dy*dz

     Vhat3Ds(ik,jk,kk) = Vhat3Ds(ik,jk,kk) &
 &    + vor3(i,j,k)**2                               &
 &    * sin(2.0d0*pi*(kx(ik)*x1b(i)+ky(jk)*x2b(j))) & 
 &    *dx*dy*dz

  enddo
  enddo
  enddo

  enddo
  enddo
  enddo

  Ehat1D(:) = 0.0d0
  dkr = dkx/sqrt(2.0d0) ! minimum k
  do ik=1,nk
  do jk=1,nk
     kr = sqrt(kx(ik)**2+ky(jk)**2)
     rk = min(nk,int(kr/dkr))
     Ehat1D(rk) = Ehat1D(rk) + sqrt(Ehat3Dc(ik,jk,kk)**2 + Ehat3Ds(ik,jk,kk)**2)*dkx*dky
     Vhat1D(rk) = Vhat1D(rk) + sqrt(Vhat3Dc(ik,jk,kk)**2 + Vhat3Ds(ik,jk,kk)**2)*dkx*dky 
  enddo
  enddo

  write(filename,'(a3,i5.5,a4)')"spc",incr,".dat"
  filename = trim(dirname)//filename
  open(unitspc,file=filename,status='replace',form='formatted')
  write(unitspc,*) "# ",time
  do rk=1,nk
     write(unitspc,'(3(1x,E12.3))') rk*dkr,Ehat1D(rk)/Etot,Vhat1D(rk)/Vtot
  enddo
  close(unitspc)

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
