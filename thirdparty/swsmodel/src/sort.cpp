//
// File: sort.cpp
//
// MATLAB Coder version            : 25.2
// C/C++ source code generated on  : 18-Nov-2025 10:58:34
//

// Include Files
#include "sort.h"
#include "rt_nonfinite.h"
#include "sortIdx.h"
#include <cmath>

// Function Definitions
//
// Arguments    : double x_data[]
//                const int x_size[2]
//                int idx_data[]
//                int idx_size[2]
// Return Type  : void
//
namespace coder {
namespace internal {
void sort(double x_data[], const int x_size[2], int idx_data[], int idx_size[2])
{
  double vwork_data[8];
  int dim;
  int npages;
  int pagesize;
  int vlen;
  int vstride;
  dim = 1;
  if (x_size[0] != 1) {
    dim = 0;
  }
  vlen = x_size[dim];
  idx_size[0] = static_cast<signed char>(x_size[0]);
  idx_size[1] = static_cast<signed char>(x_size[1]);
  vstride = 1;
  for (int k{0}; k < dim; k++) {
    vstride *= x_size[0];
  }
  npages = 1;
  dim += 2;
  for (int k{dim}; k < 3; k++) {
    npages *= x_size[1];
  }
  pagesize = vlen * vstride;
  for (int i{0}; i < npages; i++) {
    int pageoffset;
    pageoffset = i * pagesize;
    for (int j{0}; j < vstride; j++) {
      int iidx_data[8];
      int idx0;
      idx0 = pageoffset + j;
      for (int k{0}; k < vlen; k++) {
        vwork_data[k] = x_data[idx0 + k * vstride];
        iidx_data[k] = 0;
      }
      if (vlen != 0) {
        double xwork_data[8];
        double x4[4];
        int iwork_data[8];
        int i1;
        int i2;
        int i4;
        int ib;
        int nNaNs;
        signed char idx4[4];
        x4[0] = 0.0;
        idx4[0] = 0;
        x4[1] = 0.0;
        idx4[1] = 0;
        x4[2] = 0.0;
        idx4[2] = 0;
        x4[3] = 0.0;
        idx4[3] = 0;
        nNaNs = 0;
        ib = 0;
        for (int k{0}; k < vlen; k++) {
          iwork_data[k] = 0;
          if (std::isnan(vwork_data[k])) {
            dim = (vlen - nNaNs) - 1;
            iidx_data[dim] = k + 1;
            xwork_data[dim] = vwork_data[k];
            nNaNs++;
          } else {
            ib++;
            idx4[ib - 1] = static_cast<signed char>(k + 1);
            x4[ib - 1] = vwork_data[k];
            if (ib == 4) {
              double d;
              double d1;
              int b_i;
              int b_i1;
              dim = k - nNaNs;
              if (x4[0] <= x4[1]) {
                i1 = 1;
                i2 = 2;
              } else {
                i1 = 2;
                i2 = 1;
              }
              if (x4[2] <= x4[3]) {
                ib = 3;
                i4 = 4;
              } else {
                ib = 4;
                i4 = 3;
              }
              d = x4[i1 - 1];
              d1 = x4[ib - 1];
              if (d <= d1) {
                if (x4[i2 - 1] <= d1) {
                  b_i = i1;
                  b_i1 = i2;
                  i1 = ib;
                  i2 = i4;
                } else if (x4[i2 - 1] <= x4[i4 - 1]) {
                  b_i = i1;
                  b_i1 = ib;
                  i1 = i2;
                  i2 = i4;
                } else {
                  b_i = i1;
                  b_i1 = ib;
                  i1 = i4;
                }
              } else if (d <= x4[i4 - 1]) {
                if (x4[i2 - 1] <= x4[i4 - 1]) {
                  b_i = ib;
                  b_i1 = i1;
                  i1 = i2;
                  i2 = i4;
                } else {
                  b_i = ib;
                  b_i1 = i1;
                  i1 = i4;
                }
              } else {
                b_i = ib;
                b_i1 = i4;
              }
              iidx_data[dim - 3] = idx4[b_i - 1];
              iidx_data[dim - 2] = idx4[b_i1 - 1];
              iidx_data[dim - 1] = idx4[i1 - 1];
              iidx_data[dim] = idx4[i2 - 1];
              vwork_data[dim - 3] = x4[b_i - 1];
              vwork_data[dim - 2] = x4[b_i1 - 1];
              vwork_data[dim - 1] = x4[i1 - 1];
              vwork_data[dim] = x4[i2 - 1];
              ib = 0;
            }
          }
        }
        i4 = vlen - nNaNs;
        if (ib > 0) {
          signed char perm[4];
          perm[1] = 0;
          perm[2] = 0;
          perm[3] = 0;
          if (ib == 1) {
            perm[0] = 1;
          } else if (ib == 2) {
            if (x4[0] <= x4[1]) {
              perm[0] = 1;
              perm[1] = 2;
            } else {
              perm[0] = 2;
              perm[1] = 1;
            }
          } else if (x4[0] <= x4[1]) {
            if (x4[1] <= x4[2]) {
              perm[0] = 1;
              perm[1] = 2;
              perm[2] = 3;
            } else if (x4[0] <= x4[2]) {
              perm[0] = 1;
              perm[1] = 3;
              perm[2] = 2;
            } else {
              perm[0] = 3;
              perm[1] = 1;
              perm[2] = 2;
            }
          } else if (x4[0] <= x4[2]) {
            perm[0] = 2;
            perm[1] = 1;
            perm[2] = 3;
          } else if (x4[1] <= x4[2]) {
            perm[0] = 2;
            perm[1] = 3;
            perm[2] = 1;
          } else {
            perm[0] = 3;
            perm[1] = 2;
            perm[2] = 1;
          }
          dim = static_cast<unsigned char>(ib);
          for (int k{0}; k < dim; k++) {
            i1 = (i4 - ib) + k;
            i2 = perm[k];
            iidx_data[i1] = idx4[i2 - 1];
            vwork_data[i1] = x4[i2 - 1];
          }
        }
        dim = nNaNs >> 1;
        for (int k{0}; k < dim; k++) {
          i1 = i4 + k;
          i2 = iidx_data[i1];
          ib = (vlen - k) - 1;
          iidx_data[i1] = iidx_data[ib];
          iidx_data[ib] = i2;
          vwork_data[i1] = xwork_data[ib];
          vwork_data[ib] = xwork_data[i1];
        }
        if ((static_cast<unsigned int>(nNaNs) & 1U) != 0U) {
          dim += i4;
          vwork_data[dim] = xwork_data[dim];
        }
        if (i4 > 1) {
          i2 = i4 >> 2;
          ib = 4;
          while (i2 > 1) {
            if ((static_cast<unsigned int>(i2) & 1U) != 0U) {
              i2--;
              dim = ib * i2;
              i1 = i4 - dim;
              if (i1 > ib) {
                merge(iidx_data, vwork_data, dim, ib, i1 - ib, iwork_data,
                      xwork_data);
              }
            }
            dim = ib << 1;
            i2 >>= 1;
            for (int k{0}; k < i2; k++) {
              merge(iidx_data, vwork_data, k * dim, ib, ib, iwork_data,
                    xwork_data);
            }
            ib = dim;
          }
          if (i4 > ib) {
            merge(iidx_data, vwork_data, 0, ib, i4 - ib, iwork_data,
                  xwork_data);
          }
        }
      }
      for (int k{0}; k < vlen; k++) {
        dim = idx0 + k * vstride;
        x_data[dim] = vwork_data[k];
        idx_data[dim] = iidx_data[k];
      }
    }
  }
}

} // namespace internal
} // namespace coder

//
// File trailer for sort.cpp
//
// [EOF]
//
