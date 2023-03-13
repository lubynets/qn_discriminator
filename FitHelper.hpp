#ifndef FIT_HELPER_H
#define FIT_HELPER_H

#include "TF1.h"
#include "TGraphErrors.h"
#include "TMatrixDSym.h"
#include "TMatrixD.h"

inline float EvalError(double* x, std::pair<TF1*, TMatrixDSym*> f_and_cov) {// add check if npar of func is equal to dim cov
  const int Npar = f_and_cov.first->GetNpar();
  TMatrixD dfdp(Npar, 1);
  for (int i = 0; i < Npar; i++) {
    dfdp[i][0] = f_and_cov.first->GradientPar(i, x);
  }

  TMatrixD dfdp_T = dfdp;
  dfdp_T.T();

  float result = std::sqrt((dfdp_T * (*f_and_cov.second) * dfdp)[0][0]);

  if(!std::isfinite(result)) {
    result = 0.;
  }

  return result;
}

inline TGraphErrors* FuncWithErrors(std::pair<TF1*, TMatrixDSym*> f_and_cov) {
  TGraphErrors* graph = new TGraphErrors();
  const int Nsteps = 1000;
  const float left = f_and_cov.first->GetXmin();
  const float right = f_and_cov.first->GetXmax();
  const float step = (right - left) / Nsteps;

  double x = left;
  int i = 0;
  while (x <= right) {
    const float y = f_and_cov.first->Eval(x);
    float ey;
    if (f_and_cov.second != nullptr) {
      ey = EvalError(&x, f_and_cov);
    }
    else {
      ey = 0.;
    }
    graph->SetPoint(i, x, y);
    graph->SetPointError(i, 0, ey);
    x += step;
    i++;
  }

  return graph;
}

#endif//FIT_HELPER_H
