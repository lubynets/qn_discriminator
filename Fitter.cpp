#include "Fitter.hpp"

#include "FitHelper.hpp"

#include "TCanvas.h"
#include "TFile.h"
#include "TFitResult.h"

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

  // -------------- begin fitting of BootStrap samples ---------------------
  TF1* fu{nullptr};
  MyFunctor functu(shape_);
  functu.SetMu(mu_);
  for(auto& bs : fbs_) {
    fu = new TF1("fu", functu, graphleft - 0.0001, graphright + 0.0001, Npar);
    bs.bs_graph_v_ -> Fit("fu", "0");
    bs.bs_fit_chi2_ = fu->GetChisquare();
    for (int i = 0; i < Npar; i++) {
      bs.bs_fit_params_.push_back(fu->GetParameter(i));
    }
    delete fu;
  }
  // --------------- end fitting of BootStrap samples ----------------------
}

void Fitter::AddBsGraphToFit(TGraphErrors* graph) {
  fbs_.push_back(FitterBootStrap());
  fbs_.back().bs_graph_v_ = graph;
}

void Fitter::SetBsGraphsToFit(const std::vector<TGraphErrors*> graphs) {
  for(auto& g : graphs) {
    AddBsGraphToFit(g);
  }
}

void Fitter::Print() const {
  std::cout << "\nFitter::Print()\n";
  std::cout << "Main fit:\n";
  std::cout << "Fit parameters:\n";
  for(auto& par : fit_params_){
    std::cout << par << "\t";
  }
  std::cout << "\n";
  std::cout << "Fit parameter errors:\n";
  for(auto& per : fit_params_errors_){
    std::cout << per << "\t";
  }
  std::cout << "\n";
  std::cout << "Chi2 / NDF = " << fit_chi2_ << " / " << fit_ndf_ << "\n\n";

  if(fbs_.size() == 0) return;
  std::cout << "BootStrap samples fit:\n";
  for(int i=0; i<fbs_.at(0).bs_fit_params_.size(); i++) {
    std::cout << i << "-th parameter:\t";
    for(auto& bs : fbs_) {
      std::cout << bs.bs_fit_params_.at(i) << "\t";
    }
    std::cout << "\n";
  }
  std::cout << "Chi2:\t";
  for(auto& bs : fbs_) {
    std::cout << bs.bs_fit_chi2_ << "\t";
  }
  std::cout << "\n";
}

// float Fitter::EvalError(double* x, std::pair<TF1*, TMatrixDSym*> f_and_cov) const {// add check if npar of func is equal to dim cov
//   const int Npar = f_and_cov.first->GetNpar();
//   TMatrixD dfdp(Npar, 1);
//   for (int i = 0; i < Npar; i++) {
//     dfdp[i][0] = f_and_cov.first->GradientPar(i, x);
//   }
//
//   TMatrixD dfdp_T = dfdp;
//   dfdp_T.T();
//
//   float result = std::sqrt((dfdp_T * (*f_and_cov.second) * dfdp)[0][0]);
//
//   if(!std::isfinite(result)) {
//     result = 0.;
//   }
//
//   return result;
// }

// TGraphErrors* Fitter::FuncWithErrors(std::pair<TF1*, TMatrixDSym*> f_and_cov) const {
//   TGraphErrors* graph = new TGraphErrors();
//   const int Nsteps = 1000;
//   const float left = f_and_cov.first->GetXmin();
//   const float right = f_and_cov.first->GetXmax();
//   const float step = (right - left) / Nsteps;
//
//   double x = left;
//   int i = 0;
//   while (x <= right) {
//     const float y = f_and_cov.first->Eval(x);
//     const float ey = EvalError(&x, f_and_cov);
//     graph->SetPoint(i, x, y);
//     graph->SetPointError(i, 0, ey);
//     x += step;
//     i++;
//   }
//
//   return graph;
// }
