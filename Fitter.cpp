#include "Fitter.hpp"

#include "TFile.h"
#include "TCanvas.h"
#include "TFitResult.h"
#include "TMatrixD.h"

void Fitter::Fit()
{
  const int Npar = 3;
  
  const float graphleft = graph_v_ -> GetPointX(0);                     //TODO read this info from Qn::Axis
  const float graphright = graph_v_ -> GetPointX(graph_v_->GetN()-1);

  MyFunctor funct(shape_);
  funct.SetExpectedMu(mu_);
  TF1* f = new TF1("f", funct, graphleft-0.0001, graphright+0.0001, Npar);
  TMatrixDSym* cov = new TMatrixDSym(f->GetNpar());
     
  TFitResultPtr frptr = graph_v_ -> Fit("f", "S0");
  *cov = frptr -> GetCovarianceMatrix();
  
  v_fit_ = {f, cov};
  graph_fit_ = FuncWithErrors(v_fit_);
  
  fit_chi2_ = f->GetChisquare();
  fit_ndf_ = f->GetNDF();

  std::string formula = "[0] + [1]*(x-" + std::to_string(mu_) + ")";          // TODO replace std::to_string(mu_) with another parameter [2]
  TF1* f2 = new TF1("f2", formula.c_str(), graphleft-0.01, graphright+0.01);  // contribution from bckgr to flow
  
  for(int i=0; i<Npar; i++)
  {
    if(i!=0)
      f2->SetParameter(i-1, f->GetParameter(i));
    
    fit_params_.push_back(f->GetParameter(i));
    fit_params_errors_.push_back(f->GetParError(i));
  }
  
  f2 -> SetLineColor(kBlue);
  graph_v_ -> GetListOfFunctions() -> Add(f2);
}

float Fitter::EvalError(double* x, std::pair<TF1*, TMatrixDSym*> f_and_cov) const            // add check if npar of func is equal to dim cov
{
  const int Npar = f_and_cov.first->GetNpar();
  TMatrixD dfdp(Npar, 1);
  for(int i=0; i<Npar; i++)
    dfdp[i][0] = f_and_cov.first->GradientPar(i, x);
  
  TMatrixD dfdp_T = dfdp;
  dfdp_T.T();
  
  return std::sqrt((dfdp_T*(*f_and_cov.second)*dfdp)[0][0]);
}

TGraphErrors* Fitter::FuncWithErrors(std::pair<TF1*, TMatrixDSym*> f_and_cov) const
{
  TGraphErrors* graph = new TGraphErrors();
  const int Nsteps = 1000;
  const float left = f_and_cov.first->GetXmin();
  const float right = f_and_cov.first->GetXmax();
  const float step = (right-left)/Nsteps;
  
  double x = left;
  int i = 0;
  while(x<=right)
  {
    const float y = f_and_cov.first->Eval(x);
    const float ey = EvalError(&x, f_and_cov);
    graph->SetPoint(i, x, y);
    graph->SetPointError(i, 0, ey);
    x += step;
    i++;
  }  
  
  return graph;
}
