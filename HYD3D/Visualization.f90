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
     call Snap2D
     call VISITOUT3D
!     call Fourier
!     call Probability
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
     allocate(  kin(in,jn,kn))
     allocate(   hk(in,jn,kn))
     is_inited = .true.
  endif

  do k=ks,ke
  do j=js,je
  do i=is,ie
     vor1(i,j,k)= (v3(i,j  ,k)-v3(i,j-1,k))/dy*0.5 &
                &+(v3(i,j+1,k)-v3(i,j  ,k))/dy*0.5 &
                &-(v2(i,j,k  )-v2(i,j,k-1))/dz*0.5 &
                &-(v2(i,j,k+1)-v2(i,j,k  ))/dz*0.5 
     vor2(i,j,k)= (v1(i,j,k  )-v1(i,j,k-1))/dz*0.5 &
                &+(v1(i,j,k+1)-v1(i,j,k  ))/dz*0.5 &
                &-(v3(i  ,j,k)-v3(i-1,j,k))/dx*0.5 &
                &-(v3(i+1,j,k)-v3(i  ,j,k))/dx*0.5 
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


      subroutine VISITOUT3D
      use fieldmod
      use hdf5
      implicit none 
      integer,parameter::ndim=3                   ! max size of the dimension for output file
      integer(hid_t) :: file_id                   ! file identifier
      integer ::   rank                           ! dataset rank
      integer(hsize_t), dimension(ndim) :: dims   ! dataset dimensions
      character(80) :: dsetname                   ! dataset name
      integer :: error
      character*(80) :: fname,fpath
      character*(80) :: fnxmf,fpxmf
      integer unitxmf
      integer,save:: imax,jmax,kmax
      character*(80) :: DIM3D,DIM1D,TIMEEXP

      character(40),parameter:: vstdir = "./vstdata/" 

      character(2):: id

      integer,parameter:: nfld=5
      integer,save:: vnum
      data vnum /0/

! Variable to esatimate NODE value

      real*8, dimension(:,:,:,:), allocatable,save:: V3DCEND
      real*4, dimension(:,:,:,:), allocatable,save:: V3DCENF
      real*4,dimension(:),allocatable,save::XF,YF,ZF

      integer:: i,j,k

      logical,save::is_inited
      data is_inited / .false. /

      id="tb"

! Initialize
      if(.not. is_inited)then
         imax = in-2*igs
         jmax = jn-2*jgs
         kmax = kn-2*kgs

         allocate(XF(imax+1))
         allocate(YF(jmax+1))
         allocate(ZF(kmax+1))
         allocate(V3DCENF   (1:imax  ,1:jmax  ,1:kmax  ,nfld))
         
         XF(1:imax+1)= real(x1a(1:imax+1))
         YF(1:jmax+1)= real(x2a(1:jmax+1))
         ZF(1:kmax+1)= real(x3a(1:kmax+1))
  
         call makedirs(vstdir)

         is_inited = .true.
      endif
      
!      if(.not. mod(incr,10) ==0 )return
      vnum=vnum+1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Writing HDF5 file
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Determine the output filename
            write(fname,"(a2,i5.5,a3)") &
     &  id,vnum,'.h5'
      fpath =  trim(vstdir)//fname
! Open HDF5
      CALL h5open_f (error)
! Note:
! H5F_ACC_TRUNC_F: delete the file if that is exist
! H5F_ACC_EXCL_F : fail writing to the file if that  is exist
      CALL h5fcreate_f(fpath, H5F_ACC_TRUNC_F, file_id, error)
      dims(:) = 0
      rank = 1
      dims(1) = imax+1
      dsetname="/x"
      call  hdf5OUT(file_id,rank,dims,dsetname,XF)
      rank=1
      dims(1) = jmax+1
      dsetname="/y"
      call  hdf5OUT(file_id,rank,dims,dsetname,YF)
      rank=1
      dims(1) = kmax+1
      dsetname="/z"
      call  hdf5OUT(file_id,rank,dims,dsetname,ZF)

      do k=1,kmax
      do j=1,jmax
      do i=1,imax
         V3DCENF(i,j,k,1) = real( kin(igs+i,jgs+j,kgs+k) )
         V3DCENF(i,j,k,2) = real(  hk(igs+i,jgs+j,kgs+k) )
         V3DCENF(i,j,k,3) = real(  v1(igs+i,jgs+j,kgs+k) )
         V3DCENF(i,j,k,4) = real(  v2(igs+i,jgs+j,kgs+k) )
         V3DCENF(i,j,k,5) = real(  v3(igs+i,jgs+j,kgs+k) )
      enddo
      enddo
      enddo

      rank=3
      dims(1) = imax
      dims(2) = jmax
      dims(3) = kmax
      dsetname="/kinene"
      call  hdf5OUT(file_id,rank,dims,dsetname,V3DCENF(1,1,1,1))
      dsetname="/kinhel"
      call  hdf5OUT(file_id,rank,dims,dsetname,V3DCENF(1,1,1,2))

      dsetname="/vx"
      call  hdf5OUT(file_id,rank,dims,dsetname,V3DCENF(1,1,1,3))
      dsetname="/vy"
      call  hdf5OUT(file_id,rank,dims,dsetname,V3DCENF(1,1,1,4))
      dsetname="/vz"
      call  hdf5OUT(file_id,rank,dims,dsetname,V3DCENF(1,1,1,5))

! Close the File
      call  h5Fclose_f(file_id,error)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Writing XDMF file
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Determine the output filename
            write(fnxmf,"(a2,i5.5,a4)") & 
     &  id,vnum,'.xmf'
      fpxmf =  trim(vstdir)//fnxmf
! FileOpen
      unitxmf=1221
      open(unitxmf,file=fpxmf,status='unknown',form='formatted')
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Header
      write(unitxmf,'(a)')'<?xml version="1.0" ?>'
      write(unitxmf,'(a)')'<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" []>'
      write(unitxmf,'(a)')'<Xdmf Version="2.0">'
      write(unitxmf,'(a)')'  <Domain>'
      write(unitxmf,'(a)')'    <Grid Name="mesh" GridType="Uniform">'
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Time
      write(TIMEEXP,"(I0)") int((time))
      write(unitxmf,'(a)')'      <Time Value="'//trim(TIMEEXP)//'"/>'
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Coordinate
      write(DIM3D,"(I0,1x,I0,1x,I0)") kmax+1,jmax+1,imax+1
      write(unitxmf,'(a)')'      <Topology' &
     &               //' TopologyType="3DRectMesh" NumberOfElements="' &
     &               //trim(DIM3D)//'"/>'
      write(unitxmf,'(a)')'      <Geometry GeometryType="VXVYVZ">'
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! x coordinate
      write(DIM1D,"(I0)") imax+1
      dsetname = "x"
      write(unitxmf,'(a)')&
     &                '        <DataItem Dimensions="' // trim(DIM1D)&
     &              //'" Name="' // trim(dsetname) // '"'
      write(unitxmf,'(a)') '          ' // &
     &               'NumberType="Float" Precision="4" Format="HDF">'
      write(unitxmf,'(a)')'	  ' // trim(fname) //":/x"
      write(unitxmf,'(a)')'        </DataItem>'
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! y coordinate
      write(DIM1D,"(I0)") jmax+1
      dsetname="y"
      write(unitxmf,*)'        <DataItem Dimensions="'//trim(DIM1D)&
     &              //'" Name="' // trim(dsetname) // '"'
      write(unitxmf,*) '          ' // &
     &              'NumberType="Float" Precision="4" Format="HDF">'
      write(unitxmf,'(a)')'	  '//trim(fname)//":/y"
      write(unitxmf,'(a)')'        </DataItem>'
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! z coordinate
      write(DIM1D,"(I0)") kmax+1
      dsetname="z"
      write(unitxmf,'(a)')'        <DataItem Dimensions="'//trim(DIM1D) &
     &              //'" Name="'//trim(dsetname) // '"'
      write(unitxmf,'(a)') '          ' // &
     &              'NumberType="Float" Precision="4" Format="HDF">'
      write(unitxmf,'(a)')'	  '//trim(fname)//":/z"
      write(unitxmf,'(a)')'        </DataItem>'
      write(unitxmf,'(a)')'      </Geometry>'
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! kinetic energy
      dsetname="kinene"
      write(unitxmf,'(a)')'      <Attribute Name="'//trim(dsetname) &
     &              //'" AttributeType="Scalar" Center="Cell">'
      write(DIM3D,"(I0,1x,I0,1x,I0)") kmax,jmax,imax
      write(unitxmf,'(a)')'        ' &
     &                    //'<DataItem Dimensions="'//trim(DIM3D)//'"'
      write(unitxmf,'(a)')'          ' // &
     &                'NumberType="Float" Precision="4" Format="HDF">'
      write(unitxmf,'(a)')'	  '//trim(fname)//":/kinene"
      write(unitxmf,'(a)')'        </DataItem>'
      write(unitxmf,'(a)')'      </Attribute>'
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! kinematic helicity
      dsetname="kinhel"
      write(unitxmf,'(a)')'      <Attribute Name="'//trim(dsetname) &
     &              //'" AttributeType="Scalar" Center="Cell">'
      write(DIM3D,"(I0,1x,I0,1x,I0)") kmax,jmax,imax
      write(unitxmf,'(a)')'        ' &
     &                    //'<DataItem Dimensions="'//trim(DIM3D)//'"'
      write(unitxmf,'(a)')'          ' // &
     &                'NumberType="Float" Precision="4" Format="HDF">'
      write(unitxmf,'(a)')'	  '//trim(fname)//":/kinhel"
      write(unitxmf,'(a)')'        </DataItem>'
      write(unitxmf,'(a)')'      </Attribute>'

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! vx
      dsetname="vx"
      write(unitxmf,'(a)')'      <Attribute Name="'//trim(dsetname) &
     &              //'" AttributeType="Scalar" Center="Cell">'
      write(DIM3D,"(I0,1x,I0,1x,I0)") kmax,jmax,imax
      write(unitxmf,'(a)')'        ' &
     &                    //'<DataItem Dimensions="'//trim(DIM3D)//'"'
      write(unitxmf,'(a)')'          ' // &
     &                'NumberType="Float" Precision="4" Format="HDF">'
      write(unitxmf,'(a)')'	  '//trim(fname)//":/vx"
      write(unitxmf,'(a)')'        </DataItem>'
      write(unitxmf,'(a)')'      </Attribute>'

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! vy
      dsetname="vy"
      write(unitxmf,'(a)')'      <Attribute Name="'//trim(dsetname) &
     &              //'" AttributeType="Scalar" Center="Cell">'
      write(DIM3D,"(I0,1x,I0,1x,I0)") kmax,jmax,imax
      write(unitxmf,'(a)')'        ' &
     &                    //'<DataItem Dimensions="'//trim(DIM3D)//'"'
      write(unitxmf,'(a)')'          ' // &
     &                'NumberType="Float" Precision="4" Format="HDF">'
      write(unitxmf,'(a)')'	  '//trim(fname)//":/vy"
      write(unitxmf,'(a)')'        </DataItem>'
      write(unitxmf,'(a)')'      </Attribute>'

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! vz
      dsetname="vz"
      write(unitxmf,'(a)')'      <Attribute Name="'//trim(dsetname) &
     &              //'" AttributeType="Scalar" Center="Cell">'
      write(DIM3D,"(I0,1x,I0,1x,I0)") kmax,jmax,imax
      write(unitxmf,'(a)')'        ' &
     &                    //'<DataItem Dimensions="'//trim(DIM3D)//'"'
      write(unitxmf,'(a)')'          ' // &
     &                'NumberType="Float" Precision="4" Format="HDF">'
      write(unitxmf,'(a)')'	  '//trim(fname)//":/vz"
      write(unitxmf,'(a)')'        </DataItem>'
      write(unitxmf,'(a)')'      </Attribute>'


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Footer
      write(unitxmf,'(a)')'    </Grid>'
      write(unitxmf,'(a)')'  </Domain>'

      write(unitxmf,'(a)')'</Xdmf>'
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      close(unitxmf)
      return
      end subroutine VISITOUT3D

      subroutine hdf5OUT(file_id,rank,dims,dsetname,dset)
      use hdf5
      implicit none
      integer(hid_t),intent(in) :: file_id                    ! file identifier
      integer,intent(in)  :: rank                             ! dataset rank
      integer(hsize_t), dimension(rank),intent(in)  :: dims   ! dataset dimensions
      character(80),intent(in):: dsetname                     ! dataset name
      real(4),intent(in) :: dset
      integer(hid_t) :: dset_id       ! dataset identifier
      integer(hid_t) :: dspace_id     ! dataspace identifier
      integer :: error
! Create the data space for the data set
      call h5Screate_simple_f(rank, dims, dspace_id, error)
! Get dset_id for data set
      call h5Dcreate_f(file_id,dsetname,h5t_native_real,dspace_id, &
     &                 dset_id,error)
! Write the data
      call h5Dwrite_f(dset_id, h5t_native_real, dset, dims, error)
! Close the data id
      call h5Dclose_f(dset_id, error)
! Close the data space
      call h5Sclose_f(dspace_id, error)

      return
      end subroutine hdf5out

      subroutine makedirs(outdir)
      character(len=*), intent(in) :: outdir
      character(len=256) command
      write(command, *) 'if [ ! -d ', trim(outdir), ' ]; then mkdir -p ', trim(outdir), '; fi'
      write(*, *) trim(command)
      call system(command)
      end subroutine makedirs

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
  k=ks
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
end subroutine Snap2D
 
