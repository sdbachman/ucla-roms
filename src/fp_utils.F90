module fp_utils

  implicit none
  private

  public :: clear_fp_exceptions, restore_fp_exceptions

contains

  subroutine clear_fp_exceptions()
    interface
      subroutine clear_fp_exceptions_c() bind(c, name="clear_fp_exceptions")
      end subroutine clear_fp_exceptions_c
      subroutine restore_fp_exceptions_c() bind(c, name="restore_fp_exceptions")
      end subroutine restore_fp_exceptions_c
    end interface
    call clear_fp_exceptions_c()
  end subroutine clear_fp_exceptions

  subroutine restore_fp_exceptions()
    interface
      subroutine restore_fp_exceptions_c() bind(c, name="restore_fp_exceptions")
      end subroutine restore_fp_exceptions_c
    end interface
    call restore_fp_exceptions_c()
  end subroutine restore_fp_exceptions

end module fp_utils
