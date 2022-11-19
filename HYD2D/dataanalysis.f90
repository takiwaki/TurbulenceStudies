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
    real(8):: dx,dy
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
     is_inited = .true.
  endif

  k=1
  do j=js,je
  do i=is,ie
     vor(i,j,k)= (v2(i  ,j,k)-v2(i-1,j,k))/dx*0.5 &
               &+(v2(i+1,j,k)-v2(i  ,j,k))/dx*0.5 &
               &-(v1(i,j  ,k)-v1(i,j-1,k))/dy*0.5 &
               &-(v1(i,j+1,k)-v1(i,j  ,k))/dy*0.5 
  enddo
  enddo

  write(filename,'(a3,i5.5,a4)')"vor",incr,".dat"
  filename = trim(dirname)//filename
  open(unitvor,file=filename,status='replace',form='formatted')

  write(unitvor,*) "# ",time
  write(unitvor,*) "# x y omega_z"
  do j=js,je
  do i=is,ie
     write(unitvor,'(3(1x,E12.3))') x1b(i),x2b(j),vor(i,j,k)
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
  real(8),dimension(nk,nk):: Ehat2Dc,Vhat2Dc,Ehat2Ds,Vhat2Ds
  real(8),dimension(nk):: kx,ky
  real(8),dimension(nk):: Ehat1D,Vhat1D
  real(8):: kr
  real(8):: dkx,dky,dkr
  character(20),parameter::dirname="output/"
  character(40)::filename
  integer,parameter::unitspc=21
  real(8):: pi

  pi=acos(-1.0d0)

  k=1
  
  Etot=0.0d0
  Vtot=0.0d0
  do j=js,je
  do i=is,ie
     Etot = Etot &
 &    + 0.5d0*d(i,j,k)                              &
 &    *(v1(i,j,k)*v1(i,j,k) + v2(i,j,k)*v2(i,j,k))  & 
 &    *dx*dy

     Vtot = Vtot &
 &    + vor(i,j,k)**2                               & 
 &    *dx*dy

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

  Ehat2Dc(:,:) = 0.0d0
  Ehat2Ds(:,:) = 0.0d0
  Vhat2Dc(:,:) = 0.0d0
  Vhat2Ds(:,:) = 0.0d0

  do ik=1,nk
  do jk=1,nk

  do j=js,je
  do i=is,ie
     Ehat2Dc(ik,jk) = Ehat2Dc(ik,jk) &
 &    + 0.5d0*d(i,j,k)                              &
 &    *(v1(i,j,k)*v1(i,j,k) + v2(i,j,k)*v2(i,j,k))  &
 &    * cos(2.0d0*pi*(kx(ik)*x1b(i)+ky(jk)*x2b(j))) & 
 &    *dx*dy
     Ehat2Ds(ik,jk) = Ehat2Ds(ik,jk) &
 &    + 0.5d0*d(i,j,k)                              &
 &    *(v1(i,j,k)*v1(i,j,k) + v2(i,j,k)*v2(i,j,k))  &
 &    * sin(2.0d0*pi*(kx(ik)*x1b(i)+ky(jk)*x2b(j))) & 
 &    *dx*dy

     Vhat2Dc(ik,jk) = Vhat2Dc(ik,jk) &
 &    + vor(i,j,k)**2                               &
 &    * cos(2.0d0*pi*(kx(ik)*x1b(i)+ky(jk)*x2b(j))) & 
 &    *dx*dy

     Vhat2Ds(ik,jk) = Vhat2Ds(ik,jk) &
 &    + vor(i,j,k)**2                               &
 &    * sin(2.0d0*pi*(kx(ik)*x1b(i)+ky(jk)*x2b(j))) & 
 &    *dx*dy

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
     Ehat1D(rk) = Ehat1D(rk) + sqrt(Ehat2Dc(ik,jk)**2 + Ehat2Ds(ik,jk)**2)*dkx*dky
     Vhat1D(rk) = Vhat1D(rk) + sqrt(Vhat2Dc(ik,jk)**2 + Vhat2Ds(ik,jk)**2)*dkx*dky 
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
