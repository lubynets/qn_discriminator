#include "ShapeFitter.hpp"

#include "TMatrixD.h"
#include "TFitResult.h"

#include <array>

ShapeFitter::ShapeFitter(TH1F* histo)
{
  histo_all_ = histo;
}

void ShapeFitter::Fit()
{
  const float left_external = ShapeFitter::mu - 15 * ShapeFitter::sigma;      //TODO make it settable
  const float left_internal = ShapeFitter::mu - 8 * ShapeFitter::sigma;
  const float right_internal = ShapeFitter::mu + 8 * ShapeFitter::sigma;
  const float right_external = ShapeFitter::mu + 15 * ShapeFitter::sigma;

  TH1F* histo_without_peak = ExcludeInterval(histo_all_, left_internal, right_internal);
  
  DefineBckgrFunc(left_external, right_external);
  FitBckgr(histo_without_peak);
  graph_bckgr_ = FuncWithErrors({bckgr_fit_, bckgr_fit_cov_});
  chi2_bckgr_fit_ = bckgr_fit_->GetChisquare() / bckgr_fit_->GetNDF();
  
  histo_sgnl_ = SubtractBckgr(histo_all_, {bckgr_fit_, bckgr_fit_cov_}, left_external, right_external);
  DefineSgnlFunc(histo_sgnl_, left_external, right_external);                 // TODO check if the same ranges as bckgr
  FitSgnl(histo_sgnl_);
  graph_sgnl_ = FuncWithErrors({sgnl_fit_, sgnl_fit_cov_});
  chi2_sgnl_fit_ = sgnl_fit_->GetChisquare() / sgnl_fit_->GetNDF();
  
  DefineAllFunc(left_external, right_external);
  FitAll();
  graph_all_ = FuncWithErrors({all_fit_, all_fit_cov_});
  chi2_all_fit_ = all_fit_->GetChisquare() / all_fit_->GetNDF();
  
  RedefineBckgrAndSgnl(left_external, right_external);
  regraph_bckgr_ = FuncWithErrors({bckgr_refit_, nullptr});
  regraph_sgnl_ = FuncWithErrors({sgnl_refit_, nullptr});  
  rehisto_bckgr_ = Func2Histo(bckgr_refit_);
  rehisto_sgnl_ = Func2Histo(sgnl_refit_);
}

TH1F* ShapeFitter::ExcludeInterval(TH1F* histo, float left, float right) const
{
  TH1F* histo_out = (TH1F*)histo->Clone();
  for(int iBin=histo_out->FindBin(left); iBin<=histo_out->FindBin(right); iBin++)
    histo_out -> SetBinContent(iBin, 0);
  
  return histo_out;
}

TH1F* ShapeFitter::SubtractBckgr(TH1F* histo, std::pair<TF1*, TMatrixDSym*> f_and_cov, float left, float right) const
{
  TH1F* histo_out = (TH1F*)histo->Clone();
  histo_out -> Sumw2();
  histo_out -> Add(f_and_cov.first, -1);
  
  for(int iBin=1; iBin<=histo_out->GetNbinsX(); iBin++)
  {
    if(histo_out->GetBinCenter(iBin)<left || histo_out->GetBinCenter(iBin)>right)
    {
      histo_out->SetBinContent(iBin, 0);
      histo_out->SetBinError(iBin, 0);
    }
    else
    {
      const float eh = histo->GetBinError(iBin);
      double binCenter = histo_out->GetBinCenter(iBin);
      const float ef = EvalError(&binCenter, f_and_cov);              // maybe create EvalErrorSq also to avoid sqrt and then **2
      const float ee = std::sqrt(eh*eh + ef*ef);
      histo_out -> SetBinError(iBin, ee);
    }
  }
      
  return histo_out;
}

float ShapeFitter::EvalError(double* x, std::pair<TF1*, TMatrixDSym*> f_and_cov) const            // add check if npar of func is equal to dim cov
{
  const int Npar = f_and_cov.first->GetNpar();
  TMatrixD dfdp(Npar, 1);
  for(int i=0; i<Npar; i++)
    dfdp[i][0] = f_and_cov.first->GradientPar(i, x);
  
  TMatrixD dfdp_T = dfdp;
  dfdp_T.T();
  
  return std::sqrt((dfdp_T*(*f_and_cov.second)*dfdp)[0][0]);
}

TGraphErrors* ShapeFitter::FuncWithErrors(std::pair<TF1*, TMatrixDSym*> f_and_cov) const
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
    float ey;
    if(f_and_cov.second != nullptr)
      ey = EvalError(&x, f_and_cov);
    else
      ey = 0.;
    graph->SetPoint(i, x, y);
    graph->SetPointError(i, 0, ey);
    x += step;
    i++;
  }  
  
  return graph;
}

TH1F* ShapeFitter::Func2Histo(TF1* func) const
{
  const int nbins = 1000;
  double left, right;
  func -> GetRange(left, right);
  TH1F* h = new TH1F();
  h -> SetBins(nbins, left, right);
  const float binwidth = (right-left)/nbins;
  for(int i=1; i<nbins; i++)
  {
    const float x = left + binwidth*(i-0.5);
    h->SetBinContent(i, func->Eval(x));
  }

  return h;
}

void ShapeFitter::FitAll()
{
  TFitResultPtr frptr = histo_all_ -> Fit(all_fit_, "RS0");
  *all_fit_cov_ = frptr -> GetCovarianceMatrix();
}

void ShapeFitter::DefineAllFunc(float left, float right)
{
  const int Npar_bckgr = 4;
  const int Npar_sgnl = 10;
  const int Npar = Npar_bckgr+Npar_sgnl;                                                        // TODO generalize # of parameters
  
  all_fit_ = new TF1("all_fit", AllShape, left, right, Npar);
  
  for(int i=0; i<Npar_bckgr; i++)
  {
    all_fit_ -> SetParameter(i, bckgr_fit_->GetParameter(i));
    double minlimit, maxlimit;
    bckgr_fit_ -> GetParLimits(i, minlimit, maxlimit);
    all_fit_ -> SetParLimits(i, minlimit, maxlimit);
  }
  
  for(int i=0; i<Npar_sgnl; i++)
  {
    all_fit_ -> SetParameter(Npar_bckgr+i, sgnl_fit_->GetParameter(i));
    double minlimit, maxlimit;
    sgnl_fit_ -> GetParLimits(i, minlimit, maxlimit);
    all_fit_ -> SetParLimits(Npar_bckgr+i, minlimit, maxlimit);
  }
  
  all_fit_cov_ = new TMatrixDSym(all_fit_->GetNpar());
}

void ShapeFitter::RedefineBckgrAndSgnl(float left, float right)
{
//   bckgr_refit_ = (TF1*) bckgr_fit_->Clone();                   // FIXME Clone is not working properly (function is cornery)
//   sgnl_refit_ = (TF1*) sgnl_fit_->Clone();                     // Accordig to root.cern info this function should NOT be cloned
  
  const int Npar_bckgr = bckgr_fit_ -> GetNpar();
  const int Npar_sgnl = sgnl_fit_ -> GetNpar();  

  bckgr_refit_ = new TF1("bckgr_refit", BckgrShape, left, right, Npar_bckgr);
  sgnl_refit_ = new TF1("sgnl_refit", SgnlShape, left, right, Npar_sgnl); 

  for(int i=0; i<Npar_bckgr; i++)
  {
    bckgr_refit_ -> SetParameter(i, all_fit_->GetParameter(i));
    double minlimit, maxlimit;
    all_fit_ -> GetParLimits(i, minlimit, maxlimit);
    bckgr_refit_ -> SetParLimits(i, minlimit, maxlimit);
  }
  
  for(int i=0; i<Npar_sgnl; i++)
  {
    sgnl_refit_ -> SetParameter(i, all_fit_->GetParameter(Npar_bckgr+i));  
    double minlimit, maxlimit;
    all_fit_ -> GetParLimits(Npar_bckgr+i, minlimit, maxlimit);
    sgnl_refit_ -> SetParLimits(i, minlimit, maxlimit);
  }
}

double AllShape(double* x, double* par)
{
  return BckgrShape(x, par) + SgnlShape(x, &par[4]);          // TODO generalize # of parameters smhw. Maybe with functors? They should have getters...
}

void ShapeFitter::FitBckgr(TH1F* histo)
{
  TFitResultPtr frptr = histo -> Fit(bckgr_fit_, "RS0");
  *bckgr_fit_cov_ = frptr -> GetCovarianceMatrix();
}

void ShapeFitter::DefineBckgrFunc(float left, float right)
{
//   bckgr_fit_ =  new TF1("bckgr_fit", "pol3", left, right);
  const int Npar = 4;
  bckgr_fit_ =  new TF1("bckgr_fit", BckgrShape, left, right, Npar);
  bckgr_fit_cov_ = new TMatrixDSym(bckgr_fit_->GetNpar());
}

double BckgrShape(double* x, double* par)
{
  return par[0] + par[1]*x[0] + par[2]*pow(x[0], 2) + par[3]*pow(x[0], 3);
}

void ShapeFitter::FitSgnl(TH1F* histo)
{
  TFitResultPtr frptr = histo -> Fit(sgnl_fit_, "RS0");
  *sgnl_fit_cov_ = frptr -> GetCovarianceMatrix();
}

//************ double-side crystal ball function ***************************************

double SgnlShape(double* x, double* par)
{
  double factor = par[0];
  double shift  = par[1]; // to be fixed at real peak position
  double mu     = par[2];
  double sigma  = par[3];
  double a1     = par[4];
  double n1     = par[5];
  double a2     = par[6];
  double n2     = par[7];
  
  double u = (x[0]-shift-mu)/sigma;
  double A1 = TMath::Power(n1/a1, n1)*TMath::Exp(-a1*a1/2);
  double A2 = TMath::Power(n2/a2, n2)*TMath::Exp(-a2*a2/2);
  double B1 = n1/a1 - a1;
  double B2 = n2/a2 - a2;
  
  if(u<-a1)
    return factor*A1*TMath::Power((B1-u), -n1);
  else if(u>=-a1 && u<a2)
    return factor*TMath::Exp(-u*u/2);
  else if(u>=a2)
    return factor*A2*TMath::Power((B2+u), -n2);
  else
    return -1.;  
}

void ShapeFitter::DefineSgnlFunc(TH1F* histo, float left, float right)
{
  const int Npar = 8;
  sgnl_fit_ = new TF1("sgnl_fit", SgnlShape, left, right, Npar);
  sgnl_fit_ -> SetParameter(0, histo->Interpolate(ShapeFitter::mu));
  sgnl_fit_ -> FixParameter(1, ShapeFitter::mu);
  sgnl_fit_ -> SetParameter(2, 0);
  sgnl_fit_ -> SetParameter(3, ShapeFitter::sigma);
  sgnl_fit_ -> SetParameter(4, 1.);
  sgnl_fit_ -> SetParameter(5, 1.);
  sgnl_fit_ -> SetParameter(6, 1.);
  sgnl_fit_ -> SetParameter(7, 1.);
  sgnl_fit_ -> SetParLimits(0, 0, 10*histo->Interpolate(ShapeFitter::mu));
  sgnl_fit_ -> SetParLimits(4, 0, 10);
  sgnl_fit_ -> SetParLimits(5, 1, 100);
  sgnl_fit_ -> SetParLimits(6, 1, 100);
  sgnl_fit_ -> SetParLimits(7, 1, 100);
    
  sgnl_fit_cov_ = new TMatrixDSym(sgnl_fit_->GetNpar());
}

//**************************************************************************************

// //********** two expos and gauss ******************************************************
// void ShapeFitter::DefineSgnlFunc(TH1F* histo, float left, float right)
// {
//   const int Npar = 8;
//   
// //   MyFunctorShape myfuncsh;
// //   sgnl_fit_ = new TF1("sgnl_fit", myfuncsh, left, right, Npar);
//   sgnl_fit_ = new TF1("sgnl_fit", SgnlShape, left, right, Npar);
//   sgnl_fit_ -> SetParameter(0, histo->Interpolate(ShapeFitter::mu) / TMath::Gaus(0, 0, ShapeFitter::sigma));
//   sgnl_fit_ -> FixParameter(1, ShapeFitter::mu);
//   sgnl_fit_ -> SetParameter(2, 0);
//   sgnl_fit_ -> SetParameter(3, ShapeFitter::sigma);
//   const float x_shift = 1.3e-3;     // point where Exp() is pre-defined
//   sgnl_fit_ -> SetParameter(4, histo->Interpolate(ShapeFitter::mu - x_shift) / TMath::Exp(-1000*x_shift));
//   sgnl_fit_ -> SetParameter(5, histo->Interpolate(ShapeFitter::mu + x_shift) / TMath::Exp(-1000*x_shift));
//   sgnl_fit_ -> SetParameter(6, 1000.);
//   sgnl_fit_ -> SetParameter(7, -1000.);
//     
//   sgnl_fit_cov_ = new TMatrixDSym(sgnl_fit_->GetNpar());
// }
//
// // double MyFunctorShape::operator()(double* x, double* par)
// double SgnlShape(double* x, double* par)
// {
//   constexpr float x0_internal_ = 0.9e-3;
//   constexpr float x0_external_ = 1.3e-3;
//   
//   const double factor_peak  = par[0];
//   const double shift        = par[1]; // to be fixed at real peak position
//   const double mu           = par[2]; // expected to be 0
//   const double sigma        = par[3];
//   const double factor_left  = par[4];
//   const double factor_right = par[5];
//   const double k_left       = par[6];
//   const double k_right      = par[7];
//   
//   const double xx = x[0] - shift;
//   
//   auto alpha = [xx]           // fraction of Peak funkcion in transition region
//   {
//     double ksi = std::abs(xx);
//     if(ksi<x0_internal_)
//       return 1.;
//     else if(ksi>x0_external_)
//       return 0.;
//     else
//       return (x0_external_ - ksi) / (x0_external_ - x0_internal_);
//   };
//   
//   if(xx < -x0_external_)
//     return factor_left*TMath::Exp(k_left * xx);
//   else if(xx > x0_external_)
//     return factor_right*TMath::Exp(k_right * xx);
//   else if(std::abs(xx) < x0_internal_)
//     return factor_peak*TMath::Gaus(xx, mu, sigma);
//   else if(xx>x0_internal_ && xx<x0_external_)
//     return alpha()*factor_peak*TMath::Gaus(xx, mu, sigma) + (1 - alpha())*factor_right*TMath::Exp(k_right * xx);
//   else if(xx<-x0_internal_ && xx>-x0_external_)
//     return alpha()*factor_peak*TMath::Gaus(xx, mu, sigma) + (1 - alpha())*factor_left*TMath::Exp(k_left * xx);
//   else
//     return 0.;
// }
// //*************************************************************************************

// //********** two expos and lorentz ****************************************************
// void ShapeFitter::DefineSgnlFunc(TH1F* histo, float left, float right)
// {
//   const int Npar = 10;
//   
//   sgnl_fit_ = new TF1("sgnl_fit", SgnlShape, left, right, Npar);
//   sgnl_fit_ -> SetParameter(0, histo->Interpolate(ShapeFitter::mu) / TMath::CauchyDist(0, 0, ShapeFitter::sigma));
//   sgnl_fit_ -> FixParameter(1, ShapeFitter::mu);
//   sgnl_fit_ -> SetParameter(2, 0);
//   sgnl_fit_ -> SetParameter(3, ShapeFitter::sigma);
//   const float x_shift = 1.3e-3;     // point where Exp() is pre-defined
//   sgnl_fit_ -> SetParameter(4, histo->Interpolate(ShapeFitter::mu - x_shift) / TMath::Exp(-1000*x_shift));
//   sgnl_fit_ -> SetParameter(5, histo->Interpolate(ShapeFitter::mu + x_shift) / TMath::Exp(-1000*x_shift));
//   sgnl_fit_ -> SetParameter(6, 1000.);
//   sgnl_fit_ -> SetParameter(7, -1000.);
//   sgnl_fit_ -> FixParameter(8, 0.);
//   sgnl_fit_ -> FixParameter(9, 0.);
//   
//   sgnl_fit_cov_ = new TMatrixDSym(sgnl_fit_->GetNpar());
// }
// 
// double SgnlShape(double* x, double* par)
// {
//   constexpr float x0_internal_ = 0.9e-3;
//   constexpr float x0_external_ = 1.3e-3;
//   
//   const double factor_peak  = par[0];
//   const double shift        = par[1]; // to be fixed at real peak position
//   const double mu           = par[2]; // expected to be 0
//   const double sigma        = par[3];
//   const double factor_left  = par[4];
//   const double factor_right = par[5];
//   const double k_left       = par[6];
//   const double k_right      = par[7];
//   const double k_left_sq    = par[8];
//   const double k_right_sq   = par[9];
//   
//   const double xx = x[0] - shift;
//   
//   auto alpha = [xx]           // fraction of Peak funkcion in transition region
//   {
//     double ksi = std::abs(xx);
//     if(ksi<x0_internal_)
//       return 1.;
//     else if(ksi>x0_external_)
//       return 0.;
//     else
//       return (x0_external_ - ksi) / (x0_external_ - x0_internal_);
//   };
//   
//   if(xx < -x0_external_)
//     return factor_left*TMath::Exp(k_left*xx + k_left_sq*xx*xx);
//   else if(xx > x0_external_)
//     return factor_right*TMath::Exp(k_right*xx + k_right_sq*xx*xx);
//   else if(std::abs(xx) < x0_internal_)
//     return factor_peak*TMath::CauchyDist(xx, mu, sigma);
//   else if(xx>x0_internal_ && xx<x0_external_)
//     return alpha()*factor_peak*TMath::CauchyDist(xx, mu, sigma) + (1 - alpha())*factor_right*TMath::Exp(k_right*xx + k_right_sq*xx*xx);
//   else if(xx<-x0_internal_ && xx>-x0_external_)
//     return alpha()*factor_peak*TMath::CauchyDist(xx, mu, sigma) + (1 - alpha())*factor_left*TMath::Exp(k_left*xx + k_left_sq*xx*xx);
//   else
//     return 0.;
// }
// //*************************************************************************************

// //********** two expos and voigt *********************************************************
// void ShapeFitter::DefineSgnlFunc(TH1F* histo, float left, float right)
// {
//   const int Npar = 9;
//   
//   sgnl_fit_ = new TF1("sgnl_fit", SgnlShape, left, right, Npar);
//   sgnl_fit_ -> SetParameter(0, histo->Interpolate(ShapeFitter::mu) / TMath::Voigt(0, ShapeFitter::sigma/2, ShapeFitter::sigma/2));
//   sgnl_fit_ -> FixParameter(1, ShapeFitter::mu);
//   sgnl_fit_ -> SetParameter(2, 0);
//   sgnl_fit_ -> SetParameter(3, ShapeFitter::sigma/2);
//   sgnl_fit_ -> SetParameter(4, ShapeFitter::sigma/2);
//   const float x_shift = 1.3e-3;     // point where Exp() is pre-defined
//   sgnl_fit_ -> SetParameter(5, histo->Interpolate(ShapeFitter::mu - x_shift) / TMath::Exp(-1000*x_shift));
//   sgnl_fit_ -> SetParameter(6, histo->Interpolate(ShapeFitter::mu + x_shift) / TMath::Exp(-1000*x_shift));
//   sgnl_fit_ -> SetParameter(7, 1000.);
//   sgnl_fit_ -> SetParameter(8, -1000.);
//   
//   sgnl_fit_cov_ = new TMatrixDSym(sgnl_fit_->GetNpar());
// }
// 
// double SgnlShape(double* x, double* par)
// {
//   constexpr float x0_internal_ = 0.9e-3;
//   constexpr float x0_external_ = 1.3e-3;
//   
//   const double factor_peak  = par[0];
//   const double shift        = par[1]; // to be fixed at real peak position
//   const double mu           = par[2]; // expected to be 0
//   const double sigma        = par[3];
//   const double gamma        = par[4];
//   const double factor_left  = par[5];
//   const double factor_right = par[6];
//   const double k_left       = par[7];
//   const double k_right      = par[8];
//   
//   const double xx = x[0] - shift;
//   
//   auto alpha = [xx]           // fraction of Peak funkcion in transition region
//   {
//     double ksi = std::abs(xx);
//     if(ksi<x0_internal_)
//       return 1.;
//     else if(ksi>x0_external_)
//       return 0.;
//     else
//       return (x0_external_ - ksi) / (x0_external_ - x0_internal_);
//   };
//   
//   if(xx < -x0_external_)
//     return factor_left*TMath::Exp(k_left * xx);
//   else if(xx > x0_external_)
//     return factor_right*TMath::Exp(k_right * xx);
//   else if(std::abs(xx) < x0_internal_)
//     return factor_peak*TMath::Voigt(xx-mu, sigma, gamma);
//   else if(xx>x0_internal_ && xx<x0_external_)
//     return alpha()*factor_peak*TMath::Voigt(xx-mu, sigma, gamma) + (1 - alpha())*factor_right*TMath::Exp(k_right * xx);
//   else if(xx<-x0_internal_ && xx>-x0_external_)
//     return alpha()*factor_peak*TMath::Voigt(xx-mu, sigma, gamma) + (1 - alpha())*factor_left*TMath::Exp(k_left * xx);
//   else
//     return 0.;
// }
// //*************************************************************************************

// //********** double gaussian **********************************************************
// void ShapeFitter::DefineSgnlFunc(TH1F* histo, float left, float right)
// {
//   const int Npar = 5;
//   
//   sgnl_fit_ = new TF1("sgnl_fit", SgnlShape, left, right, Npar);
//   sgnl_fit_ -> SetParameter(0, histo->Interpolate(ShapeFitter::mu) / TMath::Gaus(0, ShapeFitter::sigma/2, ShapeFitter::sigma) / 2.);
//   sgnl_fit_ -> FixParameter(1, ShapeFitter::mu);
//   sgnl_fit_ -> SetParameter(2, 0);
//   sgnl_fit_ -> SetParameter(3, ShapeFitter::sigma);
//   sgnl_fit_ -> SetParameter(4, ShapeFitter::sigma/2);
//   
//   sgnl_fit_cov_ = new TMatrixDSym(sgnl_fit_->GetNpar());
// }
// 
// double SgnlShape(double* x, double* par)
// {
//   const double factor = par[0];
//   const double shift  = par[1]; // to be fixed at real peak position
//   const double mu     = par[2]; // expected to be 0
//   const double sigma  = par[3];
//   const double k      = par[4];  
//   
//   return factor*(TMath::Gaus(x[0]-shift, mu-k, sigma) + TMath::Gaus(x[0]-shift, mu+k, sigma));
// }
// //*************************************************************************************

// //********** double cauchi ************************************************************
// void ShapeFitter::DefineSgnlFunc(TH1F* histo, float left, float right)
// {
//   const int Npar = 5;
//   
//   sgnl_fit_ = new TF1("sgnl_fit", SgnlShape, left, right, Npar);
//   sgnl_fit_ -> SetParameter(0, histo->Interpolate(ShapeFitter::mu) / TMath::CauchyDist(0, 0, ShapeFitter::sigma) / 2.);
//   sgnl_fit_ -> FixParameter(1, ShapeFitter::mu);
//   sgnl_fit_ -> SetParameter(2, 0);
//   sgnl_fit_ -> SetParameter(3, ShapeFitter::sigma);
//   sgnl_fit_ -> FixParameter(4, 0);
//   
//   sgnl_fit_cov_ = new TMatrixDSym(sgnl_fit_->GetNpar());
// }
// 
// double SgnlShape(double* x, double* par)
// {
//   const double factor = par[0];
//   const double shift  = par[1]; // to be fixed at real peak position
//   const double mu     = par[2]; // expected to be 0
//   const double sigma  = par[3];
//   const double k      = par[4];  
//   
//   return factor*(TMath::CauchyDist(x[0]-shift, mu-k, sigma) + TMath::CauchyDist(x[0]-shift, mu+k, sigma));
// }
// //*************************************************************************************

// //********** double voigt *************************************************************
// void ShapeFitter::DefineSgnlFunc(TH1F* histo, float left, float right)
// {
//   const int Npar = 6;
//   
//   sgnl_fit_ = new TF1("sgnl_fit", SgnlShape, left, right, Npar);
//   sgnl_fit_ -> SetParameter(0, histo->Interpolate(ShapeFitter::mu) / TMath::Voigt(ShapeFitter::sigma/4, ShapeFitter::sigma/2, ShapeFitter::sigma/2) / 2.);
//   sgnl_fit_ -> FixParameter(1, ShapeFitter::mu);
//   sgnl_fit_ -> SetParameter(2, 0);
//   sgnl_fit_ -> SetParameter(3, ShapeFitter::sigma/2);
//   sgnl_fit_ -> SetParameter(4, ShapeFitter::sigma/2);
//   sgnl_fit_ -> SetParameter(5, ShapeFitter::sigma/4);
//   
//   sgnl_fit_cov_ = new TMatrixDSym(sgnl_fit_->GetNpar());
// }
// 
// double SgnlShape(double* x, double* par)
// {
//   const double factor = par[0];
//   const double shift  = par[1]; // to be fixed at real peak position
//   const double mu     = par[2]; // expected to be 0
//   const double sigma  = par[3];
//   const double gamma  = par[4];
//   const double k      = par[5];  
//   
//   return factor*(TMath::Voigt(x[0]-shift-(mu-k), sigma, gamma) + TMath::Voigt(x[0]-shift-(mu+k), sigma, gamma));
// }
// //*************************************************************************************

// //********** voigt-revoigt *************************************************************
// void ShapeFitter::DefineSgnlFunc(TH1F* histo, float left, float right)
// {
//   const int Npar = 6;
//   
//   sgnl_fit_ = new TF1("sgnl_fit", SgnlShape, left, right, Npar);
//   sgnl_fit_ -> SetParameter(0, histo->Interpolate(ShapeFitter::mu) / TMath::Voigt(0, ShapeFitter::sigma/2, ShapeFitter::sigma/2) / 2);
//   sgnl_fit_ -> FixParameter(1, 1);
//   sgnl_fit_ -> FixParameter(2, ShapeFitter::mu);
//   sgnl_fit_ -> SetParameter(3, 0);
//   sgnl_fit_ -> SetParameter(4, ShapeFitter::sigma/2);
//   sgnl_fit_ -> SetParameter(5, ShapeFitter::sigma/2);
//   
//   sgnl_fit_cov_ = new TMatrixDSym(sgnl_fit_->GetNpar());
// }
// 
// double SgnlShape(double* x, double* par)
// {
//   const double factor_ext = par[0];
//   const double factor_int = par[1];
//   const double shift      = par[2]; // to be fixed at real peak position
//   const double mu         = par[3]; // expected to be 0
//   const double sigma      = par[4];
//   const double gamma      = par[5];
//   
//   return factor_ext*(TMath::Voigt(x[0]-shift-mu, sigma, gamma) + factor_int*TMath::Voigt(x[0]-shift-mu, gamma, sigma));
// }
// //*************************************************************************************
