#include "Fitter.hpp"

#include "TCanvas.h"
#include "TFile.h"
#include "TFitResult.h"
#include "TMatrixD.h"

void Fitter::Fit() {
  const int Npar = 3;

  const float graphleft = graph_v_->GetPointX(0);//TODO read this info from Qn::Axis
  const float graphright = graph_v_->GetPointX(graph_v_->GetN() - 1);

  MyFunctor funct(shape_);
  funct.SetMu(mu_);
  TF1* f = new TF1("f", funct, graphleft - 0.0001, graphright + 0.0001, Npar);
  TMatrixDSym* cov = new TMatrixDSym(f->GetNpar());

  TFitResultPtr frptr = graph_v_->Fit("f", "S0");
  *cov = frptr->GetCovarianceMatrix();
  cov->Print();

  v_fit_ = {f, cov};
  graph_fit_ = FuncWithErrors(v_fit_);

  fit_chi2_ = f->GetChisquare();
  fit_ndf_ = f->GetNDF();

  v_fit_bckgr_.first = new TF1("f2", "[0] + [1]*(x-[2])", graphleft - 0.01, graphright + 0.01);// contribution from bckgr to flow
  v_fit_bckgr_.first->SetParameter(2, mu_);

  v_fit_bckgr_.second = new TMatrixDSym(f->GetNpar());

  for (int i = 0; i < Npar; i++) {
    if (i != 0) {
      v_fit_bckgr_.first->SetParameter(i - 1, f->GetParameter(i));
      for(int j=1; j<Npar; j++) {
        (*(v_fit_bckgr_.second))[i-1][j-1] = (*cov)[i][j];
      }
    }

    fit_params_.push_back(f->GetParameter(i));
    fit_params_errors_.push_back(f->GetParError(i));
  }

  graph_fit_bckgr_ = FuncWithErrors(v_fit_bckgr_);
}

float Fitter::EvalError(double* x, std::pair<TF1*, TMatrixDSym*> f_and_cov) const {// add check if npar of func is equal to dim cov
  const int Npar = f_and_cov.first->GetNpar();
  TMatrixD dfdp(Npar, 1);
  for (int i = 0; i < Npar; i++) {
    dfdp[i][0] = f_and_cov.first->GradientPar(i, x);
//     dfdp[i][0] = MyGetGradientPar(f_and_cov.first, i, x[0], 0.0001);
  }

  TMatrixD dfdp_T = dfdp;
  dfdp_T.T();

  float result = std::sqrt((dfdp_T * (*f_and_cov.second) * dfdp)[0][0]);

  if(!std::isfinite(result)) {
    result = 0.;
  }

  return result;
}

TGraphErrors* Fitter::FuncWithErrors(std::pair<TF1*, TMatrixDSym*> f_and_cov) const {
  TGraphErrors* graph = new TGraphErrors();
  const int Nsteps = 1000;
  const float left = f_and_cov.first->GetXmin();
  const float right = f_and_cov.first->GetXmax();
  const float step = (right - left) / Nsteps;

  double x = left;
  int i = 0;
  while (x <= right) {
    const float y = f_and_cov.first->Eval(x);
    const float ey = EvalError(&x, f_and_cov);
    graph->SetPoint(i, x, y);
    graph->SetPointError(i, 0, ey);
    x += step;
    i++;
  }

  return graph;
}

// double Fitter::MyGetGradientPar(TF1* f, int i, double x, double eps) const {
//   const double par_backup = f->GetParameter(i);
//
//   f->SetParameter(i, par_backup+eps/2);
//   const double value_up = f->Eval(x);
//
//   f->SetParameter(i, par_backup-eps/2);
//   const double value_low = f->Eval(x);
//
//   f->SetParameter(i, par_backup);
//
//   return (value_up - value_low)/eps;
// }
