/* Mask floating-point traps before calling C libraries (NetCDF/HDF5)
   from Fortran builds using -init=snan.  Intel's Fortran runtime can
   raise forrtl error (65) when C code performs benign FP operations. */
#include <fenv.h>

#if defined(__INTEL_COMPILER) || defined(__INTEL_LLVM_COMPILER)
#include <immintrin.h>
#endif

static fenv_t saved_fenv;
static int fp_saved = 0;

void clear_fp_exceptions(void)
{
#if defined(__INTEL_COMPILER) || defined(__INTEL_LLVM_COMPILER)
  /* Mask SSE/AVX FP exceptions in MXCSR. */
  unsigned int mxcsr = _mm_getcsr();
  mxcsr &= ~0x003F; /* clear sticky status flags */
  mxcsr |= 0x1F80;  /* mask invalid, denorm, div0, overflow, underflow */
  _mm_setcsr(mxcsr);
#endif

  if (!fp_saved) {
    (void) feholdexcept(&saved_fenv);
    fp_saved = 1;
  } else {
    (void) feclearexcept(FE_ALL_EXCEPT);
  }
  (void) fedisableexcept(FE_ALL_EXCEPT);
}

void restore_fp_exceptions(void)
{
  if (fp_saved) {
    (void) fesetenv(&saved_fenv);
    fp_saved = 0;
  }
}
