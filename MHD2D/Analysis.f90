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
    real(8),dimension(:,:,:),allocatable:: vor ! vorticity
    real(8),dimension(:,:,:),allocatable:: jcd ! current density
    real(8),dimension(:,:,:),allocatable:: mpt ! magnetic potential
    real(8),dimension(:,:,:),allocatable:: hcr ! cross helicity
    real(8):: dx,dy
    real(8):: Etot,Vtot,Mtot,Ctot
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
     call CrossHelicity
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
     is_inited = .true.
  endif

  k=1
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
  enddo
  enddo

  write(filename,'(a3,i5.5,a4)')"vor",incr,".dat"
  filename = trim(dirname)//filename
  open(unitvor,file=filename,status='replace',form='formatted')

  write(unitvor,*) "# ",time
  write(unitvor,*) "# x y omega_z"
  do j=js,je
  do i=is,ie
     write(unitvor,'(4(1x,E12.3))') x1b(i),x2b(j),vor(i,j,k),jcd(i,j,k)
  enddo
     write(unitvor,*)
  enddo

  close(unitvor)


  return
end subroutine Vorticity

subroutine Potential
  use fieldmod
  implicit none
  integer::i,j,k

  character(20),parameter::dirname="output/"
  character(40)::filename
  integer,parameter::unitvor=231

  logical,save:: is_inited
  data is_inited / .false. /


  if(.not. is_inited)then
     allocate( mpt(in,jn,kn))
     is_inited = .true.
  endif
 

  mpt(:,:,:)= 0.0d0

  k=1
  do j=js,je
  do i=is+1,ie
     mpt(i,j,k) = mpt(i-1,j,k)  + (b2(i,j,k) + b2(i-1,j,k))/2.0d0*dx
  enddo
  enddo

  do j=js+1,je
  do i=is  ,ie
     mpt(i,j,k) = mpt(i,j-1,k) + (b1(i,j,k) + b1(i,j-1,k))/2.0d0*dy
  enddo
  enddo

  return
end subroutine Potential

subroutine CrossHelicity
  use fieldmod
  implicit none
  integer::i,j,k

  character(20),parameter::dirname="output/"
  character(40)::filename
  integer,parameter::unitvor=231

  logical,save:: is_inited
  data is_inited / .false. /


  if(.not. is_inited)then
     allocate( Hcr(in,jn,kn))
     is_inited = .true.
  endif
 

  Hcr(:,:,:)= 0.0d0

  k=1
  do j=js,je
  do i=is,ie
     Hcr(i,j,k) =     v1(i,j,k) * b1(i,j,k) &
                &  +  v2(i,j,k) * b2(i,j,k) &
                &  +  v3(i,j,k) * b3(i,j,k)
  enddo
  enddo

  return
end subroutine CrossHelicity


subroutine Fourier
  use fieldmod
  implicit none
  integer::i,j,k
  integer::ik,jk,kk,rk
  integer,parameter:: nk=100
  integer,parameter:: nvar=4
  real(8),dimension(nvar):: X
  real(8),dimension(nk,nk,nvar):: Xhat2Dc,Xhat2Ds
  real(8),dimension(nk):: kx,ky
  real(8),dimension(nk,nvar):: Xhat1D
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

     Mtot = Mtot &
 &    + mpt(i,j,k)**2                               & 
 &    *dx*dy

     Ctot = Ctot &
 &    + Hcr(i,j,k)                                  & 
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

  Xhat2Dc(:,:,:) = 0.0d0
  Xhat2Ds(:,:,:) = 0.0d0

  do ik=1,nk
  do jk=1,nk

  do j=js,je
  do i=is,ie
     X(1) = 0.5d0*d(i,j,k)                          &
 &    *(v1(i,j,k)*v1(i,j,k) + v2(i,j,k)*v2(i,j,k)) 
     X(2) = vor(i,j,k)**2
     X(3) = mpt(i,j,k)**2
     X(4) = hcr(i,j,k)

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
     Xhat1D(rk,1:nvar) = Xhat1D(rk,1:nvar) + sqrt(Xhat2Dc(ik,jk,1:nvar)**2 + Xhat2Ds(ik,jk,1:nvar)**2)*dkx*dky 
  enddo
  enddo

  write(filename,'(a3,i5.5,a4)')"spc",incr,".dat"
  filename = trim(dirname)//filename
  open(unitspc,file=filename,status='replace',form='formatted')
  write(unitspc,*) "# ",time
  do rk=1,nk
     write(unitspc,'(6(1x,E12.3))') rk*dkr,Xhat1D(rk,1)/Etot &
                                  &       ,Xhat1D(rk,2)/Vtot &
                                  &       ,Xhat1D(rk,3)/Mtot &
                                  &       ,Xhat1D(rk,4)/Ctot
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
