      module adaptive_timestep

      ! module for adaptively setting dt
      ! --------------------

#include "cppdefs.opt"


      use dimensions, only: nx,ny,nz
      use scalars, only: dt, iic, g, rho0, pi
      use param, only: N,Lm,Mm,LLm, MMm, mynode, ocean_grid_comm
      use grid, only: pn, pm, rmask, f, h
      use ocean_vars, only: Hz, We, Wi, FlxU, FlxV, z_r
      use mpi_f08, only: mpi_double_precision, mpi_max
      use eos_vars, only: rho1,qp1,qp2
      use comm_vars, only: grdmin,grdmax
      use mixing, only: bvf

      implicit none

      private

#include "adaptive_timestep.opt"

      real :: Cu_Adv=0    ! Advective CFL number
      real,allocatable,dimension(:,:) :: fastest_bc_wavespeed
      real,allocatable,dimension(:,:) :: fastest_bt_wavespeed
      real :: fastest_bc_wavespeed_global
      real :: fastest_bt_wavespeed_global
      real :: smallest_wavenumber

      public set_adaptive_dt

      contains

! ----------------------------------------------------------------------
      subroutine set_adaptive_dt

      ! Set the adaptive timestep and number of barotropic timesteps per
      ! baroclinic timestep, depending on the maximum IGW speed,
      ! advective CFL, and barotropic GW speed

      implicit none

      ! local
      integer :: i,j,ierr
      real :: dt_adv, dt_bc, dt_bt
      logical,save :: done=.false.

      if (.not. done) then
        allocate( fastest_bc_wavespeed(GLOBAL_2D_ARRAY) );
        fastest_bc_wavespeed(:,:) = 0

        allocate( fastest_bt_wavespeed(GLOBAL_2D_ARRAY) );
        fastest_bt_wavespeed(:,:) = 0

        smallest_wavenumber = 2*pi/(grdmax*MAX(LLm,MMm))

        ! Calculate fastest barotropic wave speed
        ! Find global maximum

        do j=1,ny
          do i=1,nx        
            fastest_bt_wavespeed(i,j) = sqrt((f(i,j)/smallest_wavenumber)**2 + g*h(i,j))
          enddo
        enddo

        call MPI_Reduce(fastest_bt_wavespeed, fastest_bt_wavespeed_global, 1, mpi_double_precision,
     &    MPI_MAX, 0, ocean_grid_comm, ierr)

        ! Broadcast rank where minimum occurs to all processes
        call MPI_Bcast(fastest_bt_wavespeed_global,1,mpi_double_precision,0,ocean_grid_comm,ierr)

        done=.true.
      endif

      ! ------ Maximum timestep based on advective CFL ------

      ! Cu_Adv = dt_old * max(Volume Flux leaving cell / Cell Volume)
      ! dt_new = max_adv_CFL * max(Cell Volume / Volume Flux leaving cell) 
      ! ----> dt_new = max_adv_CFL * dt_old / Cu_Adv
      call calc_advective_CFL

      dt_adv = max_adv_CFL * dt / MAX(Cu_Adv,1e-10)   ! Don't divide by zero 

      ! ------ Maximum timestep based on baroclinic wave speed ------ 
      call calc_bc_wave_speed

      dt_bc = max_bc_CFL * grdmin / fastest_bc_wavespeed_global

      ! ------ Maximum barotropic timestep based on barotropic wave speed ------ 
      dt_bt = max_bt_CFL * grdmin / fastest_bt_wavespeed_global

      write(*,*) "YAY1", dt_adv, dt_bc, dt_bt      

      end subroutine set_adaptive_dt !]
!-------------------------------------------------------------------
      subroutine calc_advective_CFL

! Find maximum Advective Courant number: in the code segment below it
! is defined in the finite-volume sense as the sum of fluxes directed
! OUTWARD from grid box Hz divided by its control volume. Essentially
! it is the fraction of water replaced within the grid box during one
! time step.   Note that under this definition, Cu=1 is the condition
! when flux-split first-order upstream advection scheme looses its
! positive definiteness property (hence stability).

      implicit none

      ! local
      integer :: i,j,k
      integer :: ierr
      real :: ciV, cw, cx
      real :: my_Cu_Adv

      my_Cu_Adv=0
      do k=1,nz
        do j=1,ny
          do i=1,nx
#  ifdef MASKING
            ciV=dt*rmask(i,j)*pm(i,j)*pn(i,j)/Hz(i,j,k)
#  else
            ciV=dt*pm(i,j)*pn(i,j)/Hz(i,j,k)
#  endif
            cw=ciV*( max(We(i,j,k)+Wi(i,j,k), 0.)
     &                      -min(We(i,j,k-1)+Wi(i,j,k-1), 0.) )

            cx=cw+ciV*( max(FlxU(i+1,j,k), 0.) -min(FlxU(i,j,k), 0.)
     &                 +max(FlxV(i,j+1,k), 0.) -min(FlxV(i,j,k), 0.))

            if (cx > my_Cu_Adv) then
              my_Cu_Adv=cx
            endif
          enddo
        enddo
      enddo

      ! Find global maximum
      call MPI_Reduce(my_Cu_Adv, Cu_Adv, 1, mpi_double_precision,
     &  MPI_MAX, 0, ocean_grid_comm, ierr)

      ! Broadcast rank where minimum occurs to all processes
      call MPI_Bcast(Cu_Adv,1,mpi_double_precision,0,ocean_grid_comm,ierr)

      end subroutine calc_advective_CFL !]
!-------------------------------------------------------------------
      subroutine calc_bc_wave_speed

! Find maximum baroclinic wave speed, following the Lagrangian
! formulation in Robey and Grooms (2024).
! We intend to find the eigenvalues of M, given in Eq. 4.20.
! The IGW dispersion relation is given in Eq 4.21, and the wave speed = omega / k.  
! The fastest wave speed uses smallest k = 2*pi/domain size.
! We will solve for the eigenvalues of M-1, given in A.8.
! Eigenvalues of M are 1/eigenvalues of M^-1, and have units of m2/s2.

! g'_k = g * (rho_k+1 - rho_k)/rho_0 


! DSTERF computes all eigenvalues of a symmetric tridiagonal matrix
! using the Pal-Walker-Kahan variant of the QL or QR algorithm.

! SUBROUTINE DSTERF( N, D, E, INFO )

! Arguments
!  =========

!  N       (input) INTEGER
!    The order of the matrix.  N >= 0.

!  D       (input/output) DOUBLE PRECISION array, dimension (N)
!    On entry, the n diagonal elements of the tridiagonal matrix.
!    On exit, if INFO = 0, the eigenvalues in ascending order.

!  E       (input/output) DOUBLE PRECISION array, dimension (N-1)
!    On entry, the (n-1) subdiagonal elements of the tridiagonal
!    matrix.
!    On exit, E has been destroyed.
!  INFO    (output) INTEGER
!    = 0:  successful exit
!    < 0:  if INFO = -i, the i-th argument had an illegal value
!    > 0:  the algorithm failed to find all of the eigenvalues in
!          a total of 30*N iterations; if INFO = i, then i
!          elements of E have not converged to zero.




      implicit none

      ! local
      integer :: i,j,k
      integer :: info, ierr
      real :: dpth, largest_eigenvalue
      real :: smallest_positive
      real :: fastest_local
      real, dimension(N) :: D=0
      logical, dimension(N) :: is_positive
      real, dimension(N) :: rho
      real, dimension(N-1) :: E=0
      real, dimension(2:N) :: g_prime
      real :: tot, N_sum, dist    

      do j=1,ny
        do i=1,nx

          N_sum = 0
          dist = 0
          tot = 0
          do k=1,nz
            dpth = -z_r(i,j,k)
            rho(k)=rho1(i,j,k)+qp1(i,j,k)*dpth*(1.-qp2*dpth)

            N_sum = N_sum + 0.5*(bvf(i,j,k)+bvf(i,j,k-1))*Hz(i,j,k) 
            dist = dist + Hz(i,j,k)
            if (N_sum>0) then
              tot = tot + sqrt(N_sum/dist)*dist
              N_sum = 0
              dist = 0
            endif
          enddo
          tot = tot/pi

          do k=2,nz
            g_prime(k)= g * (rho(k-1)-rho(k)) / rho0
          enddo

          ! Solve for the eigenvalues of a real, symmetric matrix that 
          ! is similar to M^-1.
          ! Note that ROMS' layer indexing is backwards from the normal
          ! convention, meaning that the vertical index increases
          ! UPWARDS. Here we will follow the convention in RR24 A.8, which
          ! is to have the shallowest layers at the top of the matrix M

          ! The eigenvalues will have units of s2/m2. The output of
          ! dsterf will put the smallest eigenvalue of M-1 FIRST, which
          ! corresponds to the largest eigenvalue of M
          E(1) = -1 / (g_prime(nz)*sqrt(Hz(i,j,nz)*Hz(i,j,nz-1)))
          D(1) = E(1) + 1/(g*Hz(i,j,nz))        
          do k=2,nz-1
            E(k) = -1 / (g_prime(nz-k+1)*sqrt(Hz(i,j,nz-k+1)*Hz(i,j,nz-k)))
            D(k) = -(E(k) + E(k-1))
          enddo
          D(nz) = -E(nz-1)

          call dsterf(nz,D,E,info)


          ! The disperion relation is
          ! \omega = \pm \sqrt(f^2 + \lambda K^2),
          ! where K^2 = k^2 + l^2
          ! so the baroclinic wave speed will be
          ! \omega/K = \pm \sqrt(f^2/K^2 + \lambda)
          ! The fastest wave speed will be associated with the smallest
          ! possible value of K^2 = (2*pi/MAX(Lx,Ly))**2          

          ! Find the smallest positive
          is_positive = D > 0.0

          if (ANY(is_positive)) then
            smallest_positive = MINVAL(D, MASK=is_positive)
!            fastest_bc_wavespeed(i,j) = sqrt((f(i,j)/smallest_wavenumber)**2 + 1/smallest_positive)
            fastest_bc_wavespeed(i,j) = sqrt(1.0/smallest_positive)
          else
            fastest_bc_wavespeed(i,j) = 0
          endif


          !if ((i==5) .and. (j==5)) then
          !    do k=1,nz
          !    write(*,*) "Y", k, D(k)
          !    enddo
          !endif

        enddo
      enddo

      fastest_local = MAXVAL(fastest_bc_wavespeed)

      ! Find global maximum
      call MPI_Reduce(fastest_local, fastest_bc_wavespeed_global, 1, mpi_double_precision,
     &  MPI_MAX, 0, ocean_grid_comm, ierr)

      ! Broadcast rank where minimum occurs to all processes
      call MPI_Bcast(fastest_bc_wavespeed_global,1,mpi_double_precision,0,ocean_grid_comm,ierr)

      end subroutine calc_bc_wave_speed !]
!-------------------------------------------------------------------
      end module adaptive_timestep
