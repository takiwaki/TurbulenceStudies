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
    real(8),dimension(:,:,:),allocatable:: b1,b2,b3,bp
    real(8),dimension(:,:,:),allocatable:: vor1,vor2,vor3
    real(8),dimension(:,:,:),allocatable:: jcd1,jcd2,jcd3
    real(8),dimension(:,:,:),allocatable:: mpt1,mpt2,mpt3
    real(8),dimension(:,:,:),allocatable:: kin,hk
    real(8),dimension(:,:,:),allocatable:: mag,hm,hmm,hcr
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
!$acc declare create(b1,b2,b3,bp)
!$acc declare create(vor1,vor2,vor3)
!$acc declare create(jcd1,jcd2,jcd3)
!$acc declare create(mpt1,mpt2,mpt3)
!$acc declare create(kin,hk)
!$acc declare create(mag,hmm,hm,hcr)

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
     write(6,*) "file number count",incr
     call ReadData
!     call GenerateProblem
     call Vorticity
     call Potential
!     call Snap2D
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
  
  in=izone+2*igs
  jn=jzone+2*jgs
  kn=kzone+2*kgs
  write(6,*) "data size",in,jn,kn
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
  read(unitbin)x3b(:),x3a(:)
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
  dz = x3b(2)-x3b(1)

!  write(6,*) "dx",dx,dy,dz

!$acc update device (in,jn,kn)
!$acc update device (is,js,ks)
!$acc update device (ie,je,ke)
  
!$acc update device (x1a,x1b)
!$acc update device (x2a,x2b)
!$acc update device (x3a,x3b)
!$acc update device (d,v1,v2,v3,p)
!$acc update device (b1,b2,b3,bp)
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
     allocate( jcd1(in,jn,kn))
     allocate( jcd2(in,jn,kn))
     allocate( jcd3(in,jn,kn))
     allocate(  kin(in,jn,kn))
     allocate(   hk(in,jn,kn))
     allocate(  mag(in,jn,kn))
     allocate(  hmm(in,jn,kn))
     allocate(  hcr(in,jn,kn))
!$acc update device (vor1,vor2,vor3)
!$acc update device (jcd1,jcd2,jcd3)
!$acc update device (kin,hk)
!$acc update device (mag,hmm,hcr)
     is_inited = .true.
  endif

!$acc kernels
!$acc loop collapse(3) independent
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

     jcd1(i,j,k)= (b3(i,j  ,k)-b3(i,j-1,k))/dy*0.5 &
                &+(b3(i,j+1,k)-b3(i,j  ,k))/dy*0.5 &
                &-(b2(i,j,k  )-b2(i,j,k-1))/dz*0.5 &
                &-(b2(i,j,k+1)-b2(i,j,k  ))/dz*0.5 
     jcd2(i,j,k)= (b1(i,j,k  )-b1(i,j,k-1))/dz*0.5 &
                &+(b1(i,j,k+1)-b1(i,j,k  ))/dz*0.5 &
                &-(b3(i  ,j,k)-b3(i-1,j,k))/dx*0.5 &
                &-(b3(i+1,j,k)-b3(i  ,j,k))/dx*0.5 
     jcd3(i,j,k)= (b2(i  ,j,k)-b2(i-1,j,k))/dx*0.5 &
                &+(b2(i+1,j,k)-b2(i  ,j,k))/dx*0.5 &
                &-(b1(i,j  ,k)-b1(i,j-1,k))/dy*0.5 &
                &-(b1(i,j+1,k)-b1(i,j  ,k))/dy*0.5

      mag(i,j,k)= 0.5d0*( b1(i,j,k)*b1(i,j,k) &
                &        +b2(i,j,k)*b2(i,j,k) &
                &        +b3(i,j,k)*b3(i,j,k))

       hmm(i,j,k)=  b1(i,j,k)*jcd1(i,j,k) &
                 & +b2(i,j,k)*jcd2(i,j,k) &
                 & +b3(i,j,k)*jcd3(i,j,k)
      hcr(i,j,k)=    ( b1(i,j,k)*v1(i,j,k) &
                &     +b2(i,j,k)*v2(i,j,k) &
                &     +b3(i,j,k)*v3(i,j,k))

  enddo
  enddo
  enddo
!$acc end kernels

!!$acc update host (kin,hk)
!!$acc update host (mag,hmm,hcr)
!  write(6,*)"debug1",kin(is,js,ks),hk(is,js,ks),mag(is,js,ks),hmm(is,js,ks),hcr(is,js,ks)
  return
end subroutine Vorticity

module potmod
  implicit none
  real(8),dimension(:,:,:),allocatable:: p,rp,r,Ap,phi
  real(8):: a1,a2,a3,a4
  real(8):: tmp1,tmp2,tmp3
  real(8):: alpha,beta
!$acc declare create(p,rp,r,Ap,phi)
!$acc declare create(a1,a2,a3,a4)
!$acc declare create(tmp1,tmp2,tmp3)
!$acc declare create(alpha,beta)
end module potmod

subroutine Potential
  use fieldmod, pressure => p
  use potmod
  implicit none
  integer::i,j,k,n
  integer:: iter
  integer,parameter:: itermax=1000
  real(8),parameter::eps=1.0d-1
  logical,save:: is_inited
  data is_inited / .false. /

  if(.not. is_inited)then
     allocate(  p(in,jn,kn))
     allocate( rp(in,jn,kn))
     allocate(  r(in,jn,kn))
     allocate( Ap(in,jn,kn))
     allocate(phi(in,jn,kn))
!$acc update device (p,rp,r,Ap,phi)
     
     allocate( mpt1(in,jn,kn))
     allocate( mpt2(in,jn,kn))
     allocate( mpt3(in,jn,kn))
!$acc update device (mpt1,mpt2,mpt3)
     allocate(   hm(in,jn,kn))
!$acc update device (hm)
     is_inited = .true.
  endif
  
  a1=1.0d0/dx**2
  a2=1.0d0/dy**2
  a3=1.0d0/dz**2
  a4=2.0d0*(a1+a2+a3)
!$acc update device (a1,a2,a3,a4)
  
  nloop: do n=1,3
!$acc kernels
     select case(n)
     case(1)
!$acc loop collapse(3) independent
        do k=ks,ke;  do j=js,je; do i=is,ie
           r(i,j,k) = -jcd1(i,j,k)
        enddo; enddo; enddo
     case(2)
!$acc loop collapse(3) independent
        do k=ks,ke;  do j=js,je; do i=is,ie
           r(i,j,k) = -jcd2(i,j,k)
        enddo; enddo; enddo    
     case(3)
!$acc loop collapse(3) independent
        do k=ks,ke;  do j=js,je; do i=is,ie
           r(i,j,k) = -jcd3(i,j,k)
        enddo; enddo; enddo   
     end select
!$acc end kernels
     
!$acc kernels
!$acc loop collapse(3) independent
     do k=ks,ke;  do j=js,je; do i=is,ie
          p(i,j,k) = r(i,j,k)
        phi(i,j,k) = 0.0d0
     enddo; enddo; enddo
!$acc end kernels
     
     itloop: do iter=1,itermax

!$acc kernels
     !==========================
     ! boundary condition
     !==========================
!$acc loop collapse(2) independent
     do k=ks,ke; do j=js,je
        p(is-1,j,k) = p(ie,j,k)
        p(ie+1,j,k) = p(is,j,k)    
     enddo;  enddo

!$acc loop collapse(2) independent
     do k=ks,ke; do i=is,ie
        p(i,js-1,k) = p(i,je,k)
        p(i,je+1,k) = p(i,js,k)    
     enddo;  enddo
     
!$acc loop collapse(2) independent
     do j=js,je; do i=is,ie
        p(i,j,ks-1) = p(i,j,ke)
        p(i,j,ke+1) = p(i,j,ks)    
     enddo;  enddo
!$acc end kernels
   
!$acc kernels           
!$acc loop collapse(3) independent
        do k=ks,ke; do j=js,je; do i=is,ie
           Ap(i,j,k) =  a1*(p(i+1,j,k)+p(i-1,j,k)) &
                     & +a2*(p(i,j+1,k)+p(i,j-1,k)) &
                     & +a3*(p(i,j,k+1)+p(i,j,k-1)) &
                     & -a4*(p(i,j,k))
        enddo; enddo;  enddo
!$acc end kernels

  tmp1=0.0d0
!$acc update device (tmp1)
!$acc kernels
!$acc loop collapse(3) reduction(+:tmp1)
        do k=ks,ke;  do j=js,je; do i=is,ie
           tmp1= tmp1 +  r(i,j,k)**2
        enddo; enddo; enddo
!$acc end kernels
!$acc update host (tmp1)
        
  tmp2=0.0d0
!$acc update device (tmp2)
!$acc kernels
!$acc loop collapse(3) reduction(+:tmp2)
        do k=ks,ke;  do j=js,je; do i=is,ie
           tmp2= tmp2 + Ap(i,j,k)*p(i,j,k)
        enddo; enddo; enddo
!$acc end kernels
!$acc update host (tmp2)

  alpha= tmp1/tmp2
!$acc update device (alpha)
        
!$acc kernels
!$acc loop collapse(3) independent
        do k=ks,ke;  do j=js,je; do i=is,ie
           phi(i,j,k) = phi(i,j,k) + alpha * p(i,j,k)
            rp(i,j,k) =   r(i,j,k) - alpha *Ap(i,j,k)
        enddo; enddo; enddo
!$acc end kernels
        
  tmp3=0.0d0
!$acc update device (tmp3)
!$acc kernels
!$acc loop collapse(3) reduction(+:tmp3)
        do k=ks,ke;  do j=js,je; do i=is,ie
           tmp3 = tmp3 + rp(i,j,k)**2      
        enddo; enddo; enddo
!$acc end kernels
!$acc update host (tmp3)
        
        if(tmp3 .lt. eps) then
           write(6,*) "jcd",n," found in",iter,"iterations"
           exit itloop
        endif
        
        if(iter .eq. itermax) then
           write(6,*) "iteration reaches max",itermax
           write(6,*) "error=",tmp3
        endif
  beta =  tmp3/tmp1
!$acc update device (beta)
!$acc kernels
!$acc loop collapse(3) independent
        do k=ks,ke;  do j=js,je; do i=is,ie
            p(i,j,k) = rp(i,j,k) + beta * p(i,j,k)
            r(i,j,k) = rp(i,j,k)
        enddo; enddo; enddo
!$acc end kernels
!        write(6,*) tmp1,tmp2,tmp3,alpha,beta
!        stop
     enddo itloop
     
!$acc kernels
     select case(n)
     case(1)
!$acc loop collapse(3) independent
        do k=ks,ke;  do j=js,je; do i=is,ie
           mpt1(i,j,k) = phi(i,j,k)
        enddo; enddo; enddo
     case(2)
!$acc loop collapse(3) independent
        do k=ks,ke;  do j=js,je; do i=is,ie
           mpt2(i,j,k) = phi(i,j,k)
        enddo; enddo; enddo
     case(3)
!$acc loop collapse(3) independent
        do k=ks,ke;  do j=js,je; do i=is,ie
           mpt3(i,j,k) = phi(i,j,k)
        enddo; enddo; enddo
     end select
!$acc end kernels
  enddo nloop

!$acc kernels
!$acc loop collapse(3) independent
  do k=ks,ke
  do j=js,je
  do i=is,ie
       hm(i,j,k)=  b1(i,j,k)*mpt1(i,j,k) &
                & +b2(i,j,k)*mpt2(i,j,k) &
                & +b3(i,j,k)*mpt3(i,j,k)
  enddo
  enddo
  enddo
!$acc end kernels

  return
end subroutine Potential

module spctrmod
implicit none
  real(8),dimension(:,:,:,:),allocatable:: X3D
!$acc declare create(X3D)
  integer,parameter:: nk=128
  integer,parameter:: nvar=5
  real(8),dimension(nvar):: Xtot
  real(8),dimension(nk,nk,nk,nvar):: Xhat3DC,Xhat3DS
  real(8),dimension(nk):: kx,ky,kz
  real(8),dimension(nk,nvar):: Xhat1D
  real(8):: kr
  real(8):: dkx,dky,dkz,dkr
  
  real(8) :: pi
!$acc declare create(Xtot)
!$acc declare create(dkx,dky,dkz)
!$acc declare create(kx,ky,kz)
!$acc declare create(Xhat3DC,Xhat3DS,Xhat1D)
!$acc declare create(pi)
end module spctrmod
subroutine Fourier
  use fieldmod
  use spctrmod
  implicit none
  integer::i,j,k,n
  integer::ik,jk,kk,rk
  real(8):: Xtotloc
  real(8):: Xhat3DCloc,Xhat3DSloc
  character(20),parameter::dirname="output/"
  character(40)::filename
  integer,parameter::unitspc=21
  integer,parameter::unittot=22
  integer,save::nout
  data nout / 1 /
  logical,save:: is_inited
  data is_inited / .false. /
  
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
     X3D(i,j,k,3) =  hm(i,j,k)
!     X3D(i,j,k,3) = hmm(i,j,k)
     X3D(i,j,k,4) = Hcr(i,j,k)
     X3D(i,j,k,5) = mag(i,j,k)
  enddo
  enddo
  enddo
!$acc end kernels  

!!$acc update host (X3D)
!  write(6,*)"debug2",X3D(is,js,ks,1:nvar)
 
!$acc kernels
!$acc loop independent
  do n=1,nvar
     Xtotloc = 0.0d0
!$acc loop collapse(3) reduction(+:Xtotloc)
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
! 
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
  do n=1,nvar
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
  enddo
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

  
  write(filename,'(a3,i5.5,a4)')"spc",nout,".dat"
  filename = trim(dirname)//filename
  open(unitspc,file=filename,status='replace',form='formatted')
  write(unitspc,*) "# ",time
  do rk=1,nk
     write(unitspc,'(7(1x,E12.3))') rk*dkr,Xhat1D(rk,1)/Xtot(1) &
                                         &,Xhat1D(rk,2) &
                                         &,Xhat1D(rk,3) &
                                         &,Xhat1D(rk,4) &
                                         &,Xhat1D(rk,5)/Xtot(5)
  enddo
  close(unitspc)

  write(filename,'(a3,i5.5,a4)')"tot",nout,".dat"
  filename = trim(dirname)//filename
  open(unittot,file=filename,status='replace',form='formatted')
  write(unittot,'(6(1x,E12.3))') time,Xtot(1),Xtot(2),Xtot(3),Xtot(4),Xtot(5)
  close(unittot)

  nout=nout+1
  
  return
end subroutine Fourier

subroutine Snap2D
  use fieldmod
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
  k=ks
  write(unitvor,*) "# x y omega_z"
  do j=js,je
  do i=is,ie
     write(unitvor,'(10(1x,E12.3))') x1b(i),x2b(j) &
                                 & ,jcd3(i,j,k),mpt3(i,j,k)

  enddo
     write(unitvor,*)
  enddo

  close(unitvor)


  return
end subroutine Snap2D

      subroutine GenerateProblem
      use fieldmod
      implicit none
      integer::i,j,k
      real(8)::pi
      real(8)::Ahl,Bhl,Chl
      real(8),parameter::k_ini=2.0d0

      real(8),dimension(in,jn,kn)::vpsi1b,vpsi2b
      real(8),dimension(in,jn,kn)::mpsi1b,mpsi2b
      real(8):: psinorm

      integer::seedsize
      integer,allocatable:: seed(:)
      real(8)::x

      real(8),parameter:: b0 = 2.0d0
      real(8),parameter:: x1max = 0.5d0, x2max = 0.5d0,x3max = 0.5d0
      real(8),parameter:: x1min =-0.5d0, x2min =-0.5d0,x3min =-0.5d0

      pi=acos(-1.0d0)

      Ahl = 0.5d0
      Bhl = 0.5d0
      Chl = 0.5d0

      do k=ks,ke
      do j=js,je
      do i=is,ie
         b1(i,j,k) = b0*(  Ahl*sin(2.0d0*pi*(k_ini*x3b(k)/(x3max-x3min))) &
   &                     + Chl*cos(2.0d0*pi*(k_ini*x2b(j)/(x2max-x2min))))
         b2(i,j,k) = b0*(  Bhl*sin(2.0d0*pi*(k_ini*x1b(i)/(x1max-x1min))) &
   &                     + Ahl*cos(2.0d0*pi*(k_ini*x3b(k)/(x3max-x3min))))
         b3(i,j,k) = b0*(  Chl*sin(2.0d0*pi*(k_ini*x2b(j)/(x2max-x2min))) &
   &                     + Bhl*cos(2.0d0*pi*(k_ini*x1b(i)/(x1max-x1min))))
!         write(6,*)b1(i,j,k),b2(i,j,k),b3(i,j,k)
      enddo
      enddo
      enddo

      call BoundaryCondition

      return
      end subroutine GenerateProblem


      subroutine BoundaryCondition
      use fieldmod
      implicit none
      integer::i,j,k

! x-boundary
      do k=1,kn-1
      do j=1,jn-1
      do i=1,igs
          b1(i,j,k) = b1(ie-igs+i,j,k)
          b2(i,j,k) = b2(ie-igs+i,j,k)
          b3(i,j,k) = b3(ie-igs+i,j,k)
      enddo
      enddo
      enddo

      
      do k=1,kn-1
      do j=1,jn-1
      do i=1,igs
          b1(ie+i,j,k) = b1(is+i-1,j,k)
          b2(ie+i,j,k) = b2(is+i-1,j,k)
          b3(ie+i,j,k) = b3(is+i-1,j,k)
      enddo
      enddo
      enddo

! y-boundary   
      do k=1,kn-1
      do i=1,in-1
      do j=1,jgs
          b1(i,j,k) = b1(i,je-jgs+j,k)
          b2(i,j,k) = b2(i,je-jgs+j,k)
          b3(i,j,k) = b3(i,je-jgs+j,k)
      enddo
      enddo
      enddo

      
      do k=1,kn-1
      do i=1,in-1
      do j=1,jgs
          b1(i,je+j,k) = b1(i,js+j-1,k)
          b2(i,je+j,k) = b2(i,js+j-1,k)
          b3(i,je+j,k) = b3(i,js+j-1,k)
      enddo
      enddo
      enddo

! z-boundary
      do j=1,jn-1
      do i=1,in-1
      do k=1,kgs
          b1(i,j,k) = b1(i,j,ke-kgs+k)
          b2(i,j,k) = b2(i,j,ke-kgs+k)
          b3(i,j,k) = b3(i,j,ke-kgs+k)
      enddo
      enddo
      enddo

      
      do j=1,jn-1
      do i=1,in-1
      do k=1,kgs
          b1(i,j,ke+k) = b1(i,j,ks+k-1)
          b2(i,j,ke+k) = b2(i,j,ks+k-1)
          b3(i,j,ke+k) = b3(i,j,ks+k-1)
      enddo
      enddo
      enddo

      return
      end subroutine BoundaryCondition
