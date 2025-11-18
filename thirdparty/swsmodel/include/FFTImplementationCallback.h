//
// File: FFTImplementationCallback.h
//
// MATLAB Coder version            : 25.2
// C/C++ source code generated on  : 18-Nov-2025 10:58:34
//

#ifndef FFTIMPLEMENTATIONCALLBACK_H
#define FFTIMPLEMENTATIONCALLBACK_H

// Include Files
#include "rtwtypes.h"
#include "coder_array.h"
#include <cstddef>
#include <cstdlib>

// Type Definitions
namespace coder {
namespace internal {
namespace fft {
class FFTImplementationCallback {
public:
  static int get_algo_sizes(int nfft, boolean_T useRadix2, int &nRows);
  static void r2br_r2dit_trig_impl(const array<double, 1U> &x,
                                   int unsigned_nRows,
                                   const array<double, 2U> &costab,
                                   const array<double, 2U> &sintab,
                                   array<creal_T, 1U> &y);
  static void r2br_r2dit_trig_impl(const array<creal_T, 1U> &x,
                                   int unsigned_nRows,
                                   const array<double, 2U> &costab,
                                   const array<double, 2U> &sintab,
                                   array<creal_T, 1U> &y);
  static void doHalfLengthBluestein(
      const array<double, 1U> &x, array<creal_T, 1U> &y, int nrowsx, int nRows,
      int nfft, const array<creal_T, 1U> &wwc, const array<double, 2U> &costab,
      const array<double, 2U> &sintab, const array<double, 2U> &costabinv,
      const array<double, 2U> &sintabinv);
  static void b_doHalfLengthBluestein(
      const array<double, 1U> &x, array<creal_T, 1U> &y, int nrowsx, int nRows,
      int nfft, const array<creal_T, 1U> &wwc, const array<double, 2U> &costab,
      const array<double, 2U> &sintab, const array<double, 2U> &costabinv,
      const array<double, 2U> &sintabinv);

protected:
  static void generate_twiddle_tables(int nRows, array<double, 2U> &costab,
                                      array<double, 2U> &sintab,
                                      array<double, 2U> &sintabinv);
  static void get_half_twiddle_tables(
      const array<double, 2U> &costab, const array<double, 2U> &sintab,
      const array<double, 2U> &costabinv, const array<double, 2U> &sintabinv,
      array<double, 2U> &hcostab, array<double, 2U> &hsintab,
      array<double, 2U> &hcostabinv, array<double, 2U> &hsintabinv);
  static void doHalfLengthRadix2(const array<double, 1U> &x,
                                 array<creal_T, 1U> &y, int unsigned_nRows,
                                 const array<double, 2U> &costab,
                                 const array<double, 2U> &sintab);
};

} // namespace fft
} // namespace internal
} // namespace coder

#endif
//
// File trailer for FFTImplementationCallback.h
//
// [EOF]
//
