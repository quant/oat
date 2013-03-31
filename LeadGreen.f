! version %I% of %G%                    
!Compute transmission and reflection matrix
subroutine green(m,Nx,Vxy,be,ef,p,nume,w,vr,nume2,w2,vr2,ldvr,tt,rt)
!subroutine green(m,Nx,Vxy,be,ef,p,nume,w,vr,ldvr,tt,rt)
  use params
  implicit none

  integer, intent(in) :: m,Nx      !number of points along y and along x 
  integer, intent(in) :: nume(1:2*m)      !number of points along y and along x 
  integer, intent(in) :: nume2(1:2*m)      !number of points along y and along x 
  real(WP), intent(in) :: vxy(1:Nx,1:m)  !potential in the channel
  real(WP), intent(in) :: be    !cell size
  real(WP), intent(in) :: ef    !Fermi energy
  complex(WP), intent(out) :: tt(1:m,1:m), rt(1:m,1:m)
  complex(WP), intent(in) :: w(2*m) !eigenvalues
  complex(WP), intent(in) :: w2(2*m) !eigenvalues
  complex(WP), intent(in) :: p(m) !eigenvalues
  integer, intent(in) :: ldvr   !shall be >= 2*n
  complex(WP), intent(in) :: vr(ldvr,2*m) !eigenvectors
  complex(WP), intent(in) :: vr2(ldvr,2*m) !eigenvectors
  integer, allocatable :: IPIV(:)
  integer :: i,j,k,l
  complex(WP), allocatable :: Gd0(:,:), Fp2(:,:), Up2(:,:), UVp2(:,:)
  complex(WP), allocatable :: Gd(:,:), Goff(:,:), RR(:,:)
  complex(WP), allocatable :: Fp(:,:), FVp(:,:), Fm(:,:), FVm(:,:)
  complex(WP), allocatable :: Up(:,:), UVp(:,:), Um(:,:), UVm(:,:),WORKI(:)
  complex(WP), allocatable :: lambdap(:), lambdam(:), lambdap2(:)
  integer :: info,LDAI
  integer, save :: LworkI
  complex(WP) csump,csump2,csumm,csum
  complex(WP), parameter::zzero=(0.0,0.0)
  complex(WP), parameter::zone = (1.0,0.0)   

! definition  of matrixs for Green recursive functions
  if (lworki .eq. 0) lworki = 64*M !TODO: lworki as estimate on m

  allocate(Up(1:m,1:m))
  allocate(Up2(1:m,1:m))
  allocate(Um(1:m,1:m))
  allocate(UVp(1:m,1:m))
  allocate(UVp2(1:m,1:m))
  allocate(UVm(1:m,1:m))
  allocate(lambdap2(1:m))
  allocate(lambdap(1:m))
  allocate(lambdam(1:m))
  allocate(WORKI(max(lworki,1)))
  allocate(IPIV(1:m))

  LDAI=m
    do j=1,m
       k=nume2(j)
       Up2(1:m,j) = Vr2(1:m,k)
       lambdap2(j) = w2(k)
    end do
    do j=1,m
       k=nume(j)
       Up(1:m,j) = Vr(1:m,k)
       lambdap(j) = w(k)
    end do
    do j=m+1,2*m
       k=nume(j)
       Um(1:m,j-m) = Vr(1:m,k)
       lambdam(j-m) = w(k)
    end do
       UVp2=Up2
       UVp=Up
       UVm=Um
        call zgetrf( M, M, Uvp, LDAI, IPIV, INFO )
123     call zgetri( M, Uvp, LDAI, IPIV, WORKI,LWORKI,INFO )
        if (info .eq. 0) then 
          if (real(worki(1)) .gt. lworki) then
            print *, 'lworki adjusted: ',lworki,'->',int(real(worki(1)))
            lworki = worki(1) 
            deallocate(worki)
            allocate(worki(lworki))
          endif
        else if (info .eq. -6) then 
          print *, 'lworki adjusted: ',lworki,'->',int(real(worki(1)))
          lworki = worki(1) 
          deallocate(worki)
          allocate(worki(lworki))
          goto 123
        else
          print *, 'cannot continue, info=', info, ' worki(1)=', worki(1)
          stop
        endif

        call zgetrf( M, M, UVp2, LDAI, IPIV, INFO )
        call zgetri( M, UVp2, LDAI, IPIV, WORKI,LWORKI,INFO )

        call zgetrf( M, M, UVm, LDAI, IPIV, INFO )
        call zgetri( M, Uvm, LDAI, IPIV, WORKI,LWORKI,INFO )
  allocate(Fp2(1:m,1:m))
  allocate(Fp(1:m,1:m))
  allocate(Fm(1:m,1:m))
  allocate(FVp(1:m,1:m))
  allocate(FVm(1:m,1:m))

        do k=1,M
        do i=1,M
        csump=(0.0_WP,0.0_WP)
        csump2=(0.0_WP,0.0_WP)
        csumm=(0.0_WP,0.0_WP)
        do j=1,M
        csump2=csump2+Up2(i,j)*lambdap2(j)*UVp2(j,k)
        csump=csump+Up(i,j)*lambdap(j)*UVp(j,k)
        csumm=csumm+Um(i,j)*lambdam(j)*UVm(j,k)
        end do
        Fp2(i,k)=csump2
        Fp(i,k)=csump
        Fm(i,k)=csumm
        end do
        end do
        FVp=Fp
        FVm=Fm
        call zgetrf( M, M, FVp, LDAI, IPIV, INFO )
        call zgetri( M, FVp, LDAI, IPIV, WORKI,LWORKI,INFO )
        call zgetrf( M, M, FVm, LDAI, IPIV, INFO )
        call zgetri( M, FVm, LDAI, IPIV, WORKI,LWORKI,INFO )
  allocate(Gd(1:m,1:m))
  allocate(Gd0(1:m,1:m))
  allocate(Goff(1:m,1:m))
  allocate(RR(1:m,1:m))

        Gd = (0.0_WP,0.0_WP)
        Goff = (0.0_WP,0.0_WP)
        RR = (0.0_WP,0.0_WP)
        do i =1,m
           Goff(i,i)=1.0_WP
        end do
        do l=Nx,2,-1
          do k=1, M
          do i=1, M
            csum = conjg(P(i))*RR(i,k)
            if (l.eq.Nx) then
               csum = csum + conjg(P(i))*Fp2(i,k) !!!!!!!!!!!!!
            end if
            if(i.eq.k) csum = csum + be*(Ef-Vxy(l,i))-4.
            if(k.eq.i+1.or.k.eq.i-1) csum = csum + 1.
            RR(i,k) = csum
          end do
          end do

        call zgetrf( M, M, RR, LDAI, IPIV, INFO )
        call zgetri( M, RR, LDAI, IPIV, WORKI,LWORKI,INFO )
          do i=1, M
          do k=1, M
            RR(k,i)=-RR(k,i)*P(i)
          end do
          end do
!          Gd=matmul(Goff,RR)
!          Goff=Gd
       call zgemm('N','N',m,m,m,zone,&
       &Goff(1:m,1:m),size(Goff,dim=1),&
       &RR(1:m,1:m),size(RR,dim=1),&
       &zzero,Gd(1:m,1:m),size(Gd,dim=1))
       Goff = Gd
          end do  ! end loop l=Nx,2,-1
          do i=1, M
          do k=1, M
          csum = conjg(P(k))*RR(k,i)+P(k)*FVm(k,i)
          if(i.eq.k) csum = csum + be*(Ef-Vxy(1,i))-4. !?????????
          if(k.eq.i+1.or.k.eq.i-1) csum = csum + 1.
            RR(k,i) = csum
          end do
          end do
!!$
        call zgetrf( M, M, RR, LDAI, IPIV, INFO )
        call zgetri( M, RR, LDAI, IPIV, WORKI,LWORKI,INFO )
!        Gd=matmul(Goff,RR)
!        Goff=Gd
       call zgemm('N','N',m,m,m,zone,&
       &Goff(1:m,1:m),size(Goff,dim=1),&
       &RR(1:m,1:m),size(RR,dim=1),&
       &zzero,Gd(1:m,1:m),size(Gd,dim=1))
       Goff = Gd
        Gd0=RR
!           do l=1,M
!           do k=1,M
!                csum=(0.0_WP,0.0_WP)
!              do j=1,M
!                csum=csum+Goff(k,j)*P(j)*(FVp(j,l)-FVm(j,l))
!              end do
!                Gd(k,l)=csum
!           end do
!           end do
          do i=1, M
          do k=1, M
            Goff(k,i) = Goff(k,i)*P(i)
            FVp(k,i)  = FVp(k,i)-FVm(k,i)
          end do
          end do
       call zgemm('N','N',m,m,m,zone,&
       &Goff(1:m,1:m),size(Goff,dim=1),&
       &FVp(1:m,1:m),size(FVp,dim=1),&
       &zzero,Gd(1:m,1:m),size(Gd,dim=1))

!           Goff=matmul(Gd,Up)
       call zgemm('N','N',m,m,m,zone,&
       &Gd(1:m,1:m),size(Gd,dim=1),&
       &Up(1:m,1:m),size(Up,dim=1),&
       &zzero,Goff(1:m,1:m),size(Goff,dim=1))
!           tt=-matmul(UVp,Goff)
       call zgemm('N','N',m,m,m,-zone,&
       &UVp2(1:m,1:m),size(UVp2,dim=1),& !!!!!!!zamenila UVp na UVp2
       &Goff(1:m,1:m),size(Goff,dim=1),&
       &zzero,tt(1:m,1:m),size(tt,dim=1))
        do l=1,M
        do k=1,M
        csum=(0.0_WP,0.0_WP)
        do j=1,M
        csum=csum-Gd0(k,j)*P(j)*FVp(j,l)
!        csum=csum-Gd0(k,j)*P(j)*(FVp(j,l)-FVm(j,l))
        end do
        if(k.eq.l) csum=csum-1
        Gd(k,l)=csum
       end do
       end do
!       Goff=matmul(Gd,Up)
       call zgemm('N','N',m,m,m,zone,&
       &Gd(1:m,1:m),size(Gd,dim=1),&
       &Up(1:m,1:m),size(Up,dim=1),&
       &zzero,Goff(1:m,1:m),size(Goff,dim=1))
!       rt=matmul(UVm,Goff)
       call zgemm('N','N',m,m,m,zone,&
       &UVm(1:m,1:m),size(UVm,dim=1),&
       &Goff(1:m,1:m),size(Goff,dim=1),&
       &zzero,rt(1:m,1:m),size(rt,dim=1))
  deallocate(Up)
  deallocate(Um)
  deallocate(UVp)
  deallocate(UVm)
  deallocate(lambdap)
  deallocate(lambdam)
  deallocate(WORKI)
  deallocate(IPIV) 
  deallocate(Fp)
  deallocate(Fm)
  deallocate(FVp)
  deallocate(FVm)
  deallocate(Gd)
  deallocate(Gd0)
  deallocate(Goff)
  deallocate(RR)

end subroutine green
