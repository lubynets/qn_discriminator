#include "Sig2BckgrSimple.hpp"

#include "Helper.hpp"

Sig2BckgrSimple::Sig2BckgrSimple(TH1F* histo) {
  histo_ = histo;
}

void Sig2BckgrSimple::Calculate() {
  const float left_yield = HistoIntegral (histo_, mu_ - left_shift_ - left_width_/2,   mu_ - left_shift_ + left_width_/2);
  const float right_yield = HistoIntegral(histo_, mu_ + right_shift_ - right_width_/2, mu_ + right_shift_ + right_width_/2);
  const float mid_yield = HistoIntegral  (histo_, mu_ - mid_width_/2,                  mu_ + mid_width_/2);

  bckgr_yield_ = (left_yield  / left_width_  * right_shift_ +
                             right_yield / right_width_ * left_shift_ ) /
                            (left_shift_ + right_shift_) * mid_width_;

  sig_yield_ = mid_yield - bckgr_yield_;
}
