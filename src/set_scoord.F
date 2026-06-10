#include "cppdefs.opt"
#ifdef SOLVE3D
      subroutine compute_scoord_curves
      ! Stretching curves Cs_w, Cs_r from theta_s, theta_b, hc.

      use scoord, only: theta_b, theta_s, cs_w, cs_r
      use scalars, only: n

      implicit none

      integer :: k
      real :: ds, sc, CSF

      ds=1.D0/dble(N)
      Cs_w(N)=0.D0
      do k=N-1,1,-1
        sc=ds*dble(k-N)
        Cs_w(k)=CSF(sc, theta_s,theta_b)
      enddo
      Cs_w(0)=-1.D0
      do k=1,N
        sc=ds*(dble(k-N)-0.5D0)
        Cs_r(k)=CSF(sc, theta_s,theta_b)
      enddo
      end subroutine compute_scoord_curves


      subroutine print_scoord_diagnostics
      ! Sample z-levels; requires hmin/hmax from setup_grid2.

      use param, only: mynode
      use comm_vars, only: hmax, hmin
      use scoord, only: hc, cs_w
      use scalars, only: n

      implicit none

      integer :: k
      real :: ds, sc, z1, zhc, z2, z3

      ds=1.D0/dble(N)

      mpi_master_only write(*,'(/1x,A/)')
     &        'Vertical S-coordinate system (z at W-points):'

      mpi_master_only write(*,
     &     '(A7,2x,A11,2x,A14,2x,A16,2x,A14,2x,A14,2x,A14)')
     &     'level', 'S-coord', 'Cs-curve', 'Z at hmin',
     &     'at hc', 'half way', 'at hmax'

      do k=N,0,-1
        sc=ds*dble(k-N)
        z1=hmin*(hc*sc + hmin*Cs_w(k))/(hc+hmin)
        zhc=0.5*hc*(sc + Cs_w(k))
        z2=0.5*hmax*(hc*sc + 0.5*hmax*Cs_w(k))/(hc+0.5*hmax)
        z3=hmax*(hc*sc + hmax*Cs_w(k))/(hc+hmax)

        mpi_master_only write(*,
     &'(I7,2X,ES11.4,2X,ES14.7,2X,ES16.7,2X,ES14.7,2X,ES14.7,2X,ES14.7)')
     &       k, ds*(k-N), Cs_w(k), z1, zhc, z2, z3
      enddo
      end subroutine print_scoord_diagnostics


      subroutine set_scoord
      call compute_scoord_curves
      call print_scoord_diagnostics
      end subroutine set_scoord


      function CSF(sc, theta_s,theta_b)
      implicit none
      real*8 CSF, sc, theta_s,theta_b,csrf

      if (theta_s > 0.D0) then
        csrf=(1.D0-cosh(theta_s*sc))/(cosh(theta_s)-1.D0)
      else
        csrf=-sc**2
      endif
      if (theta_b > 0.D0) then
        CSF=(exp(theta_b*csrf)-1.D0)/(1.D0-exp(-theta_b))
      else
        CSF=csrf
      endif
      end function CSF

#else
      subroutine set_scoord_empty
      end
#endif /* SOLVE3D */
