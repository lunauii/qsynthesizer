//
// File: upfirdnCoreImpl.cpp
//
// MATLAB Coder version            : 25.2
// C/C++ source code generated on  : 18-Nov-2025 10:58:34
//

// Include Files
#include "upfirdnCoreImpl.h"
#include "rt_nonfinite.h"
#include "swsmodel_rtwutil.h"
#include "coder_array.h"
#include "omp.h"

// Function Definitions
//
// Arguments    : double xCol
//                const array<double, 1U> &hCl
//                double LhD
//                double pD
//                double qD
//                array<double, 1U> &y
// Return Type  : void
//
namespace coder {
namespace b_signal {
namespace internal {
namespace upfirdn {
void upfirdnCoreImpl(double xCol, const array<double, 1U> &hCl, double LhD,
                     double pD, double qD, array<double, 1U> &y)
{
  array<int, 2U> hStartList;
  array<int, 2U> inpMultList;
  double yout;
  int Lh;
  int Ly;
  int b_r;
  int hClStart;
  int ihstart;
  int ixend;
  int m;
  int mm;
  int p;
  int r;
  int xdzxm;
  int xm;
  int xym;
  int ym;
  Lh = static_cast<int>(LhD);
  p = static_cast<int>(pD);
  m = static_cast<int>(LhD);
  if (static_cast<int>(LhD) > 0) {
    m = 1;
  } else if (static_cast<int>(LhD) < 0) {
    m = -1;
  }
  if (m != 0) {
    if (static_cast<int>(qD) == 0) {
      m = static_cast<int>(LhD);
    } else {
      m = static_cast<int>(LhD) -
          div_s32(static_cast<int>(LhD), static_cast<int>(qD)) *
              static_cast<int>(qD);
    }
  } else {
    if (static_cast<int>(qD) == 0) {
      m = static_cast<int>(LhD);
    } else {
      m = static_cast<int>(LhD) -
          div_s32(static_cast<int>(LhD), static_cast<int>(qD)) *
              static_cast<int>(qD);
    }
    m -= static_cast<int>(qD);
  }
  if (m != 0) {
    Ly = div_s32(static_cast<int>(LhD), static_cast<int>(qD));
  } else {
    Ly = div_s32(static_cast<int>(LhD), static_cast<int>(qD)) - 1;
  }
  y.set_size(Ly + 1);
  for (int i{0}; i <= Ly; i++) {
    y[i] = 0.0;
  }
  m = div_s32(static_cast<int>(qD), static_cast<int>(pD));
  if (static_cast<int>(pD) == 0) {
    xm = static_cast<int>(qD);
  } else {
    xm = static_cast<int>(qD) - m * static_cast<int>(pD);
  }
  xdzxm = m * static_cast<int>(pD) + xm;
  hStartList.set_size(1, static_cast<int>(pD));
  inpMultList.set_size(1, static_cast<int>(pD));
  if (static_cast<int>(static_cast<int>(pD) < 1600)) {
    for (int ii{0}; ii < p; ii++) {
      if (static_cast<int>(pD) == 0) {
        ym = ii;
      } else {
        ym = ii -
             static_cast<int>(static_cast<unsigned int>(ii) /
                              static_cast<unsigned int>(static_cast<int>(pD))) *
                 static_cast<int>(pD);
      }
      xym = xm * ym;
      if (static_cast<int>(pD) == 0) {
        r = xym;
      } else {
        r = xym - xym / static_cast<int>(pD) * static_cast<int>(pD);
      }
      hStartList[ii] = r + 1;
      inpMultList[ii] = (m * ym + div_s32(xym, static_cast<int>(pD))) + 1;
    }
  } else {
#pragma omp parallel for num_threads(omp_get_max_threads()) private(r, xym, ym)

    for (int ii = 0; ii < p; ii++) {
      if (static_cast<int>(pD) == 0) {
        ym = ii;
      } else {
        ym = ii -
             static_cast<int>(static_cast<unsigned int>(ii) /
                              static_cast<unsigned int>(static_cast<int>(pD))) *
                 static_cast<int>(pD);
      }
      xym = xm * ym;
      if (static_cast<int>(pD) == 0) {
        r = xym;
      } else {
        r = xym - xym / static_cast<int>(pD) * static_cast<int>(pD);
      }
      hStartList[ii] = r + 1;
      inpMultList[ii] = (m * ym + div_s32(xym, static_cast<int>(pD))) + 1;
    }
  }
#pragma omp parallel for num_threads(omp_get_max_threads()) private(           \
        hClStart, ixend, yout, ihstart, b_r, mm)

  for (int iout = 0; iout <= Ly; iout++) {
    hClStart = div_s32(iout, p);
    if (p == 0) {
      b_r = iout;
    } else {
      b_r = iout - iout / p * p;
    }
    ihstart = hStartList[b_r];
    ixend = xdzxm * hClStart + inpMultList[b_r];
    if (ixend > 1) {
      ihstart = hStartList[b_r] + p * (ixend - 1);
      ixend = 1;
    }
    yout = 0.0;
    if (Lh < ihstart) {
      hClStart = 0;
    } else {
      hClStart = div_s32(Lh - ihstart, p) + 1;
    }
    if (hClStart <= ixend) {
      ixend = hClStart;
    }
    hClStart = ihstart - p;
    for (mm = 0; mm < ixend; mm++) {
      yout += xCol * hCl[(hClStart + p) - 1];
    }
    y[iout] = yout;
  }
}

//
// Arguments    : const array<double, 2U> &xCol
//                const array<double, 1U> &hCl
//                double LxD
//                double LhD
//                double nChansD
//                double pD
//                double qD
//                array<double, 2U> &y
// Return Type  : void
//
void upfirdnCoreImpl(const array<double, 2U> &xCol,
                     const array<double, 1U> &hCl, double LxD, double LhD,
                     double nChansD, double pD, double qD, array<double, 2U> &y)
{
  array<int, 2U> hStartList;
  array<int, 2U> inpMultList;
  double yout;
  int Lh;
  int Lx;
  int Lxup;
  int Ly;
  int countInner1;
  int hClStart;
  int ihstart;
  int inputOffset;
  int kmax;
  int m;
  int mm;
  int nChans;
  int outOffset;
  int p;
  int r;
  int xColStart;
  int xdzxm;
  int xym;
  int ym;
  Lx = static_cast<int>(LxD);
  Lh = static_cast<int>(LhD);
  nChans = static_cast<int>(nChansD);
  p = static_cast<int>(pD);
  Lxup = static_cast<int>(pD) * (static_cast<int>(LxD) - 1) +
         static_cast<int>(LhD);
  m = Lxup;
  if (Lxup > 0) {
    m = 1;
  } else if (Lxup < 0) {
    m = -1;
  }
  if (m != 0) {
    if (static_cast<int>(qD) == 0) {
      m = Lxup;
    } else {
      m = Lxup - div_s32(Lxup, static_cast<int>(qD)) * static_cast<int>(qD);
    }
  } else {
    if (static_cast<int>(qD) == 0) {
      m = Lxup;
    } else {
      m = Lxup - div_s32(Lxup, static_cast<int>(qD)) * static_cast<int>(qD);
    }
    m -= static_cast<int>(qD);
  }
  if (m != 0) {
    Ly = div_s32(Lxup, static_cast<int>(qD)) + 1;
  } else {
    Ly = div_s32(Lxup, static_cast<int>(qD));
  }
  if (static_cast<int>(nChansD) <= 1) {
    kmax = 1;
  } else {
    kmax = static_cast<int>(nChansD);
  }
  y.set_size(Ly, kmax);
  m = Ly * kmax;
  for (int kk{0}; kk < m; kk++) {
    y[kk] = 0.0;
  }
  outOffset = 0;
  inputOffset = 0;
  m = div_s32(static_cast<int>(qD), static_cast<int>(pD));
  if (static_cast<int>(pD) == 0) {
    Lxup = static_cast<int>(qD);
  } else {
    Lxup = static_cast<int>(qD) - m * static_cast<int>(pD);
  }
  xdzxm = m * static_cast<int>(pD) + Lxup;
  hStartList.set_size(1, static_cast<int>(pD));
  inpMultList.set_size(1, static_cast<int>(pD));
  if (static_cast<int>(static_cast<int>(pD) < 1600)) {
    for (int ii{0}; ii < p; ii++) {
      if (static_cast<int>(pD) == 0) {
        ym = ii;
      } else {
        ym = ii -
             static_cast<int>(static_cast<unsigned int>(ii) /
                              static_cast<unsigned int>(static_cast<int>(pD))) *
                 static_cast<int>(pD);
      }
      xym = Lxup * ym;
      if (static_cast<int>(pD) == 0) {
        r = xym;
      } else {
        r = xym - xym / static_cast<int>(pD) * static_cast<int>(pD);
      }
      hStartList[ii] = r + 1;
      inpMultList[ii] = (m * ym + div_s32(xym, static_cast<int>(pD))) + 1;
    }
  } else {
#pragma omp parallel for num_threads(omp_get_max_threads()) private(r, xym, ym)

    for (int ii = 0; ii < p; ii++) {
      if (static_cast<int>(pD) == 0) {
        ym = ii;
      } else {
        ym = ii -
             static_cast<int>(static_cast<unsigned int>(ii) /
                              static_cast<unsigned int>(static_cast<int>(pD))) *
                 static_cast<int>(pD);
      }
      xym = Lxup * ym;
      if (static_cast<int>(pD) == 0) {
        r = xym;
      } else {
        r = xym - xym / static_cast<int>(pD) * static_cast<int>(pD);
      }
      hStartList[ii] = r + 1;
      inpMultList[ii] = (m * ym + div_s32(xym, static_cast<int>(pD))) + 1;
    }
  }
  for (int kk{0}; kk < kmax; kk++) {
#pragma omp parallel for num_threads(omp_get_max_threads()) private(           \
        hClStart, xColStart, countInner1, yout, ihstart, mm)

    for (int iout = 0; iout < Ly; iout++) {
      if (static_cast<unsigned int>(p) == 0U) {
        hClStart = MAX_int32_T;
      } else {
        hClStart = static_cast<int>(static_cast<unsigned int>(iout) /
                                    static_cast<unsigned int>(p));
      }
      if (p == 0) {
        xColStart = iout;
      } else {
        xColStart = iout - static_cast<int>(static_cast<unsigned int>(iout) /
                                            static_cast<unsigned int>(p)) *
                               p;
      }
      ihstart = hStartList[xColStart];
      hClStart = xdzxm * hClStart + inpMultList[xColStart];
      if (hClStart > Lx) {
        ihstart = hStartList[xColStart] + p * (hClStart - Lx);
        hClStart = Lx;
      }
      yout = 0.0;
      if (Lh < ihstart) {
        countInner1 = 0;
      } else {
        countInner1 = div_s32(Lh - ihstart, p) + 1;
      }
      if (countInner1 > hClStart) {
        countInner1 = hClStart;
      }
      xColStart = (inputOffset + hClStart) + 1;
      hClStart = ihstart - p;
      for (mm = 0; mm < countInner1; mm++) {
        yout += xCol[(xColStart - mm) - 2] * hCl[(hClStart + p * (mm + 1)) - 1];
      }
      y[outOffset + iout] = yout;
    }
    if (nChans != 1) {
      inputOffset += Lx;
    }
    outOffset += Ly;
  }
}

} // namespace upfirdn
} // namespace internal
} // namespace b_signal
} // namespace coder

//
// File trailer for upfirdnCoreImpl.cpp
//
// [EOF]
//
