#include "ShapeFitter.hpp"

#include "TFitResult.h"
#include "TMatrixD.h"

#include <array>

ShapeFitter::ShapeFitter(TH1F* histo) {
  histo_all_ = histo;
}

void ShapeFitter::Fit() {
  const float left_external = mu_ - 15 * sigma_;//TODO make it settable
  const float left_internal = mu_ - 8 * sigma_;
  const float right_internal = mu_ + 8 * sigma_;
  const float right_external = mu_ + 15 * sigma_;

  TH1F* histo_without_peak = ExcludeInterval(histo_all_, left_internal, right_internal);

  DefineBckgrFunc(left_external, right_external);
  FitBckgr(histo_without_peak);
  graph_bckgr_ = FuncWithErrors({bckgr_fit_, bckgr_fit_cov_});
  chi2_bckgr_fit_ = bckgr_fit_->GetChisquare() / bckgr_fit_->GetNDF();

  histo_sgnl_ = SubtractBckgr(histo_all_, {bckgr_fit_, bckgr_fit_cov_}, left_external, right_external);
  DefineSgnlFunc(histo_sgnl_, left_external, right_external);// TODO check if the same ranges as bckgr
  FitSgnl(histo_sgnl_);
  graph_sgnl_ = FuncWithErrors({sgnl_fit_, sgnl_fit_cov_});
  chi2_sgnl_fit_ = sgnl_fit_->GetChisquare() / sgnl_fit_->GetNDF();

  DefineAllFunc(left_external, right_external);
  FitAll();
  graph_all_ = FuncWithErrors({all_fit_, nullptr});
  regraph_all_ = FuncWithErrors({all_refit_, all_refit_cov_});
  chi2_all_fit_ = all_refit_->GetChisquare() / all_refit_->GetNDF();

  RedefineBckgrAndSgnl(left_external, right_external);
  regraph_bckgr_ = FuncWithErrors({bckgr_refit_, nullptr});
  regraph_sgnl_ = FuncWithErrors({sgnl_refit_, nullptr});
}

TH1F* ShapeFitter::ExcludeInterval(TH1F* histo, float left, float right) const {
  TH1F* histo_out = (TH1F*) histo->Clone();
  for (int iBin = histo_out->FindBin(left); iBin <= histo_out->FindBin(right); iBin++)
    histo_out->SetBinContent(iBin, 0);

  return histo_out;
}

TH1F* ShapeFitter::SubtractBckgr(TH1F* histo, std::pair<TF1*, TMatrixDSym*> f_and_cov, float left, float right) const {
  TH1F* histo_out = (TH1F*) histo->Clone();
  histo_out->Sumw2();
  histo_out->Add(f_and_cov.first, -1);

  for (int iBin = 1; iBin <= histo_out->GetNbinsX(); iBin++) {
    if (histo_out->GetBinCenter(iBin) < left || histo_out->GetBinCenter(iBin) > right) {
      histo_out->SetBinContent(iBin, 0);
      histo_out->SetBinError(iBin, 0);
    } else {
      const float eh = histo->GetBinError(iBin);
      double binCenter = histo_out->GetBinCenter(iBin);
      const float ef = EvalError(&binCenter, f_and_cov);// maybe create EvalErrorSq also to avoid sqrt and then **2
      const float ee = std::sqrt(eh * eh + ef * ef);
      histo_out->SetBinError(iBin, ee);
    }
  }

  return histo_out;
}

float ShapeFitter::EvalError(double* x, std::pair<TF1*, TMatrixDSym*> f_and_cov) const {// add check if npar of func is equal to dim cov
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

TGraphErrors* ShapeFitter::FuncWithErrors(std::pair<TF1*, TMatrixDSym*> f_and_cov) const {
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

void ShapeFitter::FitAll() {
  std::cout << "ShapeFitter::FitAll()\n";
  TFitResultPtr frptr = histo_all_->Fit(all_refit_, "RS0");
  *all_refit_cov_ = frptr->GetCovarianceMatrix();
  all_refit_cov_->Print();
}

void ShapeFitter::DefineAllFunc(float left, float right) {
  const int Npar_bckgr = 5;
  const int Npar_sgnl = 8;
  const int Npar = Npar_bckgr + Npar_sgnl;// TODO generalize # of parameters

  all_fit_ = new TF1("all_fit", AllShape, left, right, Npar);
  all_refit_ = new TF1("all_refit", AllShape, left, right, Npar);

  for (int i = 0; i < Npar_bckgr; i++) {
    all_fit_->SetParameter(i, bckgr_fit_->GetParameter(i));
    all_refit_->SetParameter(i, bckgr_fit_->GetParameter(i));
    double minlimit, maxlimit;
    bckgr_fit_->GetParLimits(i, minlimit, maxlimit);
    all_fit_->SetParLimits(i, minlimit, maxlimit);
    all_refit_->SetParLimits(i, minlimit, maxlimit);
  }

  for (int i = 0; i < Npar_sgnl; i++) {
    all_fit_->SetParameter(Npar_bckgr + i, sgnl_fit_->GetParameter(i));
    all_refit_->SetParameter(Npar_bckgr + i, sgnl_fit_->GetParameter(i));
    double minlimit, maxlimit;
    sgnl_fit_->GetParLimits(i, minlimit, maxlimit);
    all_fit_->SetParLimits(Npar_bckgr + i, minlimit, maxlimit);
    all_refit_->SetParLimits(Npar_bckgr + i, minlimit, maxlimit);
  }

  all_refit_cov_ = new TMatrixDSym(all_refit_->GetNpar());
}

void ShapeFitter::RedefineBckgrAndSgnl(float left, float right) {
  //   bckgr_refit_ = (TF1*) bckgr_fit_->Clone();                   // FIXME Clone is not working properly (function is cornery)
  //   sgnl_refit_ = (TF1*) sgnl_fit_->Clone();                     // Accordig to root.cern info this function should NOT be cloned

  const int Npar_bckgr = bckgr_fit_->GetNpar();
  const int Npar_sgnl = sgnl_fit_->GetNpar();

  bckgr_refit_ = new TF1("bckgr_refit", BckgrShape, left, right, Npar_bckgr);
  sgnl_refit_ = new TF1("sgnl_refit", SgnlShape, left, right, Npar_sgnl);

  for (int i = 0; i < Npar_bckgr; i++) {
    bckgr_refit_->SetParameter(i, all_refit_->GetParameter(i));
    bckgr_refit_->SetParError(i, all_refit_->GetParError(i));
  }

  for (int i = 0; i < Npar_sgnl; i++) {
    sgnl_refit_->SetParameter(i, all_refit_->GetParameter(Npar_bckgr + i));
    sgnl_refit_->SetParError(i, all_refit_->GetParError(Npar_bckgr + i));
  }
}

double AllShape(double* x, double* par) {
  return BckgrShape(x, par) + SgnlShape(x, &par[5]);// TODO generalize # of parameters smhw. Maybe with functors? They should have getters...
}

void ShapeFitter::FitBckgr(TH1F* histo) {
  std::cout << "ShapeFitter::FitBckgr()\n";
  TFitResultPtr frptr = histo->Fit(bckgr_fit_, "RS0");
  *bckgr_fit_cov_ = frptr->GetCovarianceMatrix();
  bckgr_fit_cov_->Print();
}

void ShapeFitter::DefineBckgrFunc(float left, float right) {
  const int Npar = 5;
  bckgr_fit_ = new TF1("bckgr_fit", BckgrShape, left, right, Npar);
  bckgr_fit_->FixParameter(4, mu_);
  bckgr_fit_cov_ = new TMatrixDSym(bckgr_fit_->GetNpar());
}

double BckgrShape(double* x, double* par) {
  return par[0] + par[1] * (x[0]-par[4]) + par[2] * pow((x[0]-par[4]), 2) + par[3] * pow((x[0]-par[4]), 3);
}

void ShapeFitter::FitSgnl(TH1F* histo) {
  std::cout << "ShapeFitter::FitSgnl()\n";
  TFitResultPtr frptr = histo->Fit(sgnl_fit_, "RS0");
  *sgnl_fit_cov_ = frptr->GetCovarianceMatrix();
  sgnl_fit_cov_->Print();
}

//************ double-side crystal ball function ***************************************

double SgnlShape(double* x, double* par) {
  double factor = par[0];
  double shift = par[1];// to be fixed at real peak position
  double mu = par[2];
  double sigma = par[3];
  double a1 = par[4];
  double n1 = TMath::Power(10, par[5]);
  double a2 = par[6];
  double n2 = TMath::Power(10, par[7]);

  double u = (x[0] - shift - mu) / sigma;

  if (u < -a1)
    return factor*TMath::Exp(-a1*a1/2)*TMath::Power(1-a1*(u+a1)/n1, -n1);
  else if (u >= -a1 && u < a2)
    return factor * TMath::Exp(-u * u / 2);
  else if (u >= a2)
    return factor*TMath::Exp(-a2*a2/2)*TMath::Power(1+a2*(u-a2)/n2, -n2);
  else
    return -1.;
}

void ShapeFitter::DefineSgnlFunc(TH1F* histo, float left, float right) {
  const int Npar = 8;
  sgnl_fit_ = new TF1("sgnl_fit", SgnlShape, left, right, Npar);
  sgnl_fit_->SetParameter(0, histo->Interpolate(mu_));
  sgnl_fit_->FixParameter(1, mu_);
  sgnl_fit_->SetParameter(2, 0);
  sgnl_fit_->SetParameter(3, sigma_);
  sgnl_fit_->SetParameter(4, 1.);
  sgnl_fit_->SetParameter(5, 1.);
  sgnl_fit_->SetParameter(6, 1.);
  sgnl_fit_->SetParameter(7, 1.);
  sgnl_fit_->SetParLimits(0, 0, 10 * histo->Interpolate(mu_));
  sgnl_fit_->SetParLimits(4, 0, 100);
  sgnl_fit_->SetParLimits(5, -5, 5);
  sgnl_fit_->SetParLimits(6, 0, 100);
  sgnl_fit_->SetParLimits(7, -5, 5);

  sgnl_fit_cov_ = new TMatrixDSym(sgnl_fit_->GetNpar());
}
//************ double-side crystal ball function ***************************************

// double ShapeFitter::MyGetGradientPar(TF1* f, int i, double x, double eps) const {
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
