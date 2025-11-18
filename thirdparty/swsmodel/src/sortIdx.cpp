//
// File: sortIdx.cpp
//
// MATLAB Coder version            : 25.2
// C/C++ source code generated on  : 18-Nov-2025 10:58:34
//

// Include Files
#include "sortIdx.h"
#include "rt_nonfinite.h"

// Function Definitions
//
// Arguments    : int idx_data[]
//                double x_data[]
//                int offset
//                int np
//                int nq
//                int iwork_data[]
//                double xwork_data[]
// Return Type  : void
//
namespace coder {
namespace internal {
void merge(int idx_data[], double x_data[], int offset, int np, int nq,
           int iwork_data[], double xwork_data[])
{
  if (nq != 0) {
    int iout;
    int n;
    int p;
    int q;
    n = np + nq;
    for (int j{0}; j < n; j++) {
      q = offset + j;
      iwork_data[j] = idx_data[q];
      xwork_data[j] = x_data[q];
    }
    p = 0;
    q = np;
    iout = offset - 1;
    int exitg1;
    do {
      exitg1 = 0;
      iout++;
      if (xwork_data[p] <= xwork_data[q]) {
        idx_data[iout] = iwork_data[p];
        x_data[iout] = xwork_data[p];
        if (p + 1 < np) {
          p++;
        } else {
          exitg1 = 1;
        }
      } else {
        idx_data[iout] = iwork_data[q];
        x_data[iout] = xwork_data[q];
        if (q + 1 < n) {
          q++;
        } else {
          q = iout - p;
          for (int j{p + 1}; j <= np; j++) {
            iout = q + j;
            idx_data[iout] = iwork_data[j - 1];
            x_data[iout] = xwork_data[j - 1];
          }
          exitg1 = 1;
        }
      }
    } while (exitg1 == 0);
  }
}

} // namespace internal
} // namespace coder

//
// File trailer for sortIdx.cpp
//
// [EOF]
//
