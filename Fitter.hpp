#ifndef Fitter_H
#define Fitter_H

#include "ShapeContainer.hpp"

#include "TF1.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TH1.h"
#include "TMatrixDSym.h"
#include "TString.h"

class MyFunctor// TODO rename
{
 public:
  void SetMu(float value) { mu_ = value; }

  MyFunctor(Qn::ShapeContainer* shape) {
    shape_ = shape;
  }

  virtual ~MyFunctor() = default;

  double signal_fit(double* x, double* par) {
    return par[0];
  }

  double bckgr_fit(double* x, double* par) {
    return par[0] + par[1] * (x[0] - mu_);// + par[2]*(x[0]-mu)*(x[0]-mu);
                                          //     return par[0];
  }

  double operator()(double* x, double* par) {
    return (shape_->GetSignal(x[0]) * signal_fit(x, par) + shape_->GetBackground(x[0]) * bckgr_fit(x, &par[1])) / (shape_->GetSignal(x[0]) + shape_->GetBackground(x[0]));
  }

 private:
  float mu_{1.115683};

  Qn::ShapeContainer* shape_{nullptr};
};

class Fitter {
 public:
  Fitter() = default;
  virtual ~Fitter() = default;

  void SetMu(float value) { mu_ = value; }
  void SetShape(Qn::ShapeContainer* shape) { shape_ = shape; };
  void SetGraphToFit(TGraphErrors* graph) { graph_v_ = graph; };
  TGraphErrors* GetGraphToFit() const { return graph_v_; };
  TF1* GetVFit() const { return v_fit_.first; };
  TGraphErrors* GetGraphFit() const { return graph_fit_; };
  TF1* GetBckgrFit() const { return v_fit_bckgr_.first; };
  TGraphErrors* GetBckgrGraph() const { return graph_fit_bckgr_; };
  double GetVSignal() { return fit_params_.at(0); };
  double GetVSignalError() { return fit_params_errors_.at(0); };
  double GetFitChi2() { return fit_chi2_; };
  int GetFitNdf() { return fit_ndf_; };
  double GetFitChi2Ndf() { return fit_chi2_ / fit_ndf_; };
  const std::vector<double>& GetFitParameters() { return fit_params_; };
  const std::vector<double>& GetFitErrors() { return fit_params_errors_; };
//   float EvalError(double* x, std::pair<TF1*, TMatrixDSym*> f_and_cov) const;  // TODO unite with the same function in ShapeFitter
//   TGraphErrors* FuncWithErrors(std::pair<TF1*, TMatrixDSym*> f_and_cov) const;// TODO unite with the same function in ShapeFitter

  void Fit();

 private:
  Qn::ShapeContainer* shape_{nullptr};
  TGraphErrors* graph_v_{nullptr};// to be fitted

  std::pair<TF1*, TMatrixDSym*> v_fit_{nullptr, nullptr};
  TGraphErrors* graph_fit_{nullptr};// result of fit with errors
  std::pair<TF1*, TMatrixDSym*> v_fit_bckgr_{nullptr, nullptr};// only bckgr's contribution to flow
  TGraphErrors* graph_fit_bckgr_{nullptr};// only bckgr's contribution to flow with errors

  std::vector<double> fit_params_;
  std::vector<double> fit_params_errors_;
  double fit_chi2_{-999.};
  int fit_ndf_{-999};
  float mu_{1.115683};
};

#endif//Fitter_H
