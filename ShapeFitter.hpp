#ifndef ShapeFitter_H
#define ShapeFitter_H

#include "TF1.h"
#include "TGraphErrors.h"
#include "TH1F.h"
#include "TMatrixDSym.h"

double BckgrShape(double* x, double* par);
double SgnlShape(double* x, double* par);
double AllShape(double* x, double* par);

class ShapeFitter {
 public:
  ShapeFitter(TH1F* histo);
  virtual ~ShapeFitter() = default;

  void SetExpectedMu(float value) { mu_ = value; }
  void SetExpectedSigma(float value) { sigma_ = value; }

  TH1F* GetHistoSgnl() const { return histo_sgnl_; };
  TF1* GetFuncBckgr() const { return bckgr_fit_; };
  TF1* GetReFuncBckgr() const { return bckgr_refit_; };
  //   TH1F* GetReHistoBckgr() const { return rehisto_bckgr_; };
  TGraphErrors* GetGraphBckgr() const { return graph_bckgr_; };
  TGraphErrors* GetReGraphBckgr() const { return regraph_bckgr_; };
  TF1* GetFuncSgnl() const { return sgnl_fit_; };
  TF1* GetReFuncSgnl() const { return sgnl_refit_; };
  //   TH1F* GetReHistoSgnl() const { return rehisto_sgnl_; };
  TMatrixDSym* GetCovSgnl() const { return sgnl_fit_cov_; };// TODO need after debugging?
  TGraphErrors* GetGraphSgnl() const { return graph_sgnl_; };
  TGraphErrors* GetReGraphSgnl() const { return regraph_sgnl_; };
  TF1* GetFuncAll() const { return all_fit_; };
  TF1* GetReFuncAll() const { return all_refit_; };
  TGraphErrors* GetGraphAll() const { return graph_all_; };
  TGraphErrors* GetReGraphAll() const { return regraph_all_; };
  void Fit();
  float GetChi2BckgrFit() const { return chi2_bckgr_fit_; };
  float GetChi2SgnlFit() const { return chi2_sgnl_fit_; };
  float GetChi2AllFit() const { return chi2_all_fit_; };

 private:
  void FitSgnl(TH1F* histo);                                                  // TODO move to private after debugging
  void DefineSgnlFunc(TH1F* histo, float left, float right);                  // TODO move to private after debugging
  TGraphErrors* FuncWithErrors(std::pair<TF1*, TMatrixDSym*> f_and_cov) const;// TODO move to private after debugging

  // private:

  void DefineBckgrFunc(float left, float right);
  void DefineAllFunc(float left, float right);
  void RedefineBckgrAndSgnl(float left, float right);
  TH1F* ExcludeInterval(TH1F* histo, float left, float right) const;
  void FitBckgr(TH1F* histo);
  void FitAll();
  TH1F* SubtractBckgr(TH1F* histo, std::pair<TF1*, TMatrixDSym*> f_and_cov, float left, float right) const;
  float EvalError(double* x, std::pair<TF1*, TMatrixDSym*> f_and_cov) const;
  //   TH1F* Func2Histo(TF1* func) const;

  TH1F* histo_all_{nullptr};
  TH1F* histo_sgnl_{nullptr};
  TF1* bckgr_fit_{nullptr};
  TF1* bckgr_refit_{nullptr};
  //   TH1F* rehisto_bckgr_{nullptr};        //TODO temporary solution. For some reasons functions are saved bad...
  TMatrixDSym* bckgr_fit_cov_{nullptr};
  TGraphErrors* graph_bckgr_{nullptr};
  TGraphErrors* regraph_bckgr_{nullptr};
  TF1* sgnl_fit_{nullptr};
  TF1* sgnl_refit_{nullptr};
  //   TH1F* rehisto_sgnl_{nullptr};        //TODO temporary solution. For some reasons functions are saved bad...
  TMatrixDSym* sgnl_fit_cov_{nullptr};
  TGraphErrors* graph_sgnl_{nullptr};
  TGraphErrors* regraph_sgnl_{nullptr};
  TF1* all_fit_{nullptr};
  TF1* all_refit_{nullptr};
  TMatrixDSym* all_refit_cov_{nullptr};
  TGraphErrors* graph_all_{nullptr};
  TGraphErrors* regraph_all_{nullptr};
  float chi2_bckgr_fit_{-799.};
  float chi2_sgnl_fit_{-799.};
  float chi2_all_fit_{-799.};

  float mu_{1.115683};
  float sigma_{0.00145786};
};
#endif// ShapeFitter_H
