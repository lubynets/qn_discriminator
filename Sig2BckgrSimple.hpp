#ifndef Sig2BckgrSimple_H
#define Sig2BckgrSimple_H

#include "TH1F.h"

class Sig2BckgrSimple {
public:
  Sig2BckgrSimple(TH1F* histo);
  virtual ~Sig2BckgrSimple() = default;

  void SetMu(float value) { mu_ = value; }
  void SetLeftShift(float value) { left_shift_ = value; }
  void SetRightShift(float value) { right_shift_ = value; }
  void SetLeftWidth(float value) { left_width_ = value; }
  void SetRightWidth(float value) { right_width_ = value; }
  void SetMidWidth(float value) { mid_width_ = value; }

  float GetS2B() const { return sig_yield_/bckgr_yield_; }
  float GetPurity() const { return sig_yield_/(sig_yield_+bckgr_yield_); }

  void Calculate();

private:
  TH1F* histo_{nullptr};
  float mu_;
  float left_shift_;
  float right_shift_;
  float left_width_;
  float right_width_;
  float mid_width_;

  float sig_yield_;
  float bckgr_yield_;
};

#endif // Sig2BckgrSimple_H
