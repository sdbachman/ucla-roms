      module adaptive_timestep

      ! module for adaptively setting dt
      ! --------------------

#include "cppdefs.opt"


      use dimensions, only: nx,ny,nz
      use scalars, only: dt, iic, g, rho0
      use param, only: N, mynode, ocean_grid_comm
      use grid, only: pn, pm, rmask
      use ocean_vars, only: Hz, We, Wi, FlxU, FlxV, z_r
      use mpi_f08, only: mpi_double_precision, mpi_max
      use eos_vars, only: rho1,qp1,qp2

      implicit none

      private

#include "adaptive_timestep.opt"
      real :: Cu_Adv=0    ! Advective CFL number

      public set_adaptive_dt

      contains

! ----------------------------------------------------------------------
      subroutine set_adaptive_dt

      ! Set the adaptive timestep and number of barotropic timesteps per
      ! baroclinic timestep, depending on the maximum IGW speed,
      ! advective CFL, and barotropic GW speed

      implicit none

!      integer, intent(in)                 :: N_src

      ! local
      real :: dt_adv 

      ! ------ Maximum timestep based on advective CFL ------

      ! Cu_Adv = dt_old * max(Volume Flux leaving cell / Cell Volume)
      ! dt_new = max_adv_CFL * max(Cell Volume / Volume Flux leaving cell) 
      ! ----> dt_new = max_adv_CFL * dt_old / Cu_Adv
      call calc_advective_CFL

      dt_adv = max_adv_CFL * dt / MAX(Cu_Adv,1e-10)   ! Don't divide by zero 

      write(*,*) "YAY1", dt_adv

      ! ------ Maximum timestep based on baroclinic wave speed ------ 
      call calc_bc_wave_speed




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

!      do i = 2,N
!        D(i) = i
!        E(i) = i
!      enddo

!      call dsterf(N,D,E,info)



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
      real :: dpth
      real, dimension(N) :: D=0
      real, dimension(N) :: rho
      real, dimension(N-1) :: E=0
      real, dimension(2:N) :: g_prime

      do j=1,ny
        do i=1,nx

          do k=1,nz
            dpth = -z_r(i,j,k)
            rho(k)=rho1(i,j,k)+qp1(i,j,k)*dpth*(1.-qp2*dpth)
          enddo

          do k=2,nz
            g_prime(k)= g * (rho(k-1)-rho(k)) / rho0
          enddo

          ! Note that ROMS' layer indexing is backwards from the normal
          ! convention, meaning that the vertical index increases
          ! UPWARDS. Here we will follow the convention in RR24 A.8, which
          ! is to have the shallowest layers at the top of the matrix M

!          D(1) = ((g_prime(nz)+g) /        
!     &        (g_prime(nz)*g)) / Hz(i,j,nz)
!          E(1) = -1 / (g_prime(nz)*Hz(i,j,nz-1))
!          do k=2,nz-1          
!            D(k) = ((g_prime(nz-k+1)+g_prime(nz-k+2)) / 
!     &        (g_prime(nz-k+1)*g_prime(nz-k+2))) / Hz(i,j,nz-k+1)
!            E(k) = -1 / (g_prime(nz-k+1)*Hz(i,j,nz-k))
!          enddo
!          D(N) = 1/(g_prime(2)*Hz(i,j,2))

          D(1) = 1/(g*Hz(i,j,nz)) + 1/(g_prime(nz)*Hz(i,j,nz))
          E(1) = -1 / (g_prime(nz)*Hz(i,j,nz-1))
          do k=2,nz-1
            D(k) = (1/g_prime(nz-k+1)+1/g_prime(nz-k+2)) / Hz(i,j,nz-k+1)
            E(k) = -1 / (g_prime(nz-k+1)*Hz(i,j,nz-k))
          enddo
          D(N) = 1/(g_prime(2)*Hz(i,j,2))         

          call dsterf(nz,D,E,info)

          if ((i==5) .and. (j==5)) then
            do k=2,nz-1
              write(*,*) D(k), (E(k-1)+E(k))
            enddo
          endif

        enddo
      enddo
!
!      ! Find global maximum
!      call MPI_Reduce(my_Cu_Adv, Cu_Adv, 1, mpi_double_precision,
!     &  MPI_MAX, 0, ocean_grid_comm, ierr)
!
!      ! Broadcast rank where minimum occurs to all processes
!      call MPI_Bcast(Cu_Adv,1,mpi_double_precision,0,ocean_grid_comm,ierr)

      end subroutine calc_bc_wave_speed !]
!-------------------------------------------------------------------

      end module adaptive_timestep
