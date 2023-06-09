#include "FitterSimple.hpp"

void FitterSimple::SetDataContainer(Qn::DataContainerStatCollect* datacontainer) {
  data_container_ = new Qn::DataContainerStatCalculate(*datacontainer);
}

void FitterSimple::SetEdgesOfLMR(std::array<float, 6> ar) {
  left_ranges_ = std::make_pair(ar.at(0), ar.at(1));
  mid_ranges_ = std::make_pair(ar.at(2), ar.at(3));
  right_ranges_ = std::make_pair(ar.at(4), ar.at(5));
}

void FitterSimple::AutoSetEdgesOfLMR() {
  auto binedges = data_container_->GetAxis(select_axis_).GetBinEdges();
  if(binedges.size() != 6) {
    throw std::runtime_error("Error: FitterSimple::AutoSetEdgesOfLMR() is not possible");
  } else {
    this->SetEdgesOfLMR({binedges.at(0), binedges.at(1), binedges.at(2), binedges.at(3), binedges.at(4), binedges.at(5)});
  }
}

void FitterSimple::DefineLRM() {
  dcLeft_ = data_container_->Rebin({select_axis_, 1, left_ranges_.first, left_ranges_.second});
  dcLeft_ = dcLeft_.Select({select_axis_, 1, left_ranges_.first, left_ranges_.second});

  dcRight_ = data_container_->Rebin({select_axis_, 1, right_ranges_.first, right_ranges_.second});
  dcRight_ = dcRight_.Select({select_axis_, 1, right_ranges_.first, right_ranges_.second});

  dcMid_ = data_container_->Rebin({select_axis_, 1, mid_ranges_.first, mid_ranges_.second});
  dcMid_ = dcMid_.Select({select_axis_, 1, mid_ranges_.first, mid_ranges_.second});

  out_container_ = new Qn::DataContainerStatDiscriminator();
  out_container_->AddAxes(purity_container_->GetAxes());
}

void FitterSimple::Calculate() {
  this->DefineLRM();

  const float shiftL = (mid_ranges_.first + mid_ranges_.second)/2 - (left_ranges_.first + left_ranges_.second)/2;
  const float shiftR = (right_ranges_.first + right_ranges_.second)/2 - (mid_ranges_.first + mid_ranges_.second)/2;

  const int Nsamples = dcMid_.At(0).GetSampleMeans().size();

  for(int i=0; i<dcRight_.size(); i++) {
    //-------- Mean ----------------------------
    const float vL = dcLeft_.At(i).Mean();
    const float vR = dcRight_.At(i).Mean();
    const float vMid = dcMid_.At(i).Mean();
    const float eL = dcLeft_.At(i).StdDevOfMeanFromPropagation();
    const float eR = dcRight_.At(i).StdDevOfMeanFromPropagation();
    const float eMid = dcMid_.At(i).StdDevOfMeanFromPropagation();
    const float purity = purity_container_->At(i).MeanFromPropagation();

    const float vS = CalculateValue(vL, vR, vMid, purity, shiftL, shiftR);
    const float eS = CalculateError(eL, eR, eMid, purity, shiftL, shiftR);
    const float wS = dcMid_.At(i).SumWeights() * purity;

    out_container_->At(i).SetVEW(vS, eS, wS);
    //------------------------------------------

    //------- Bootstrap ------------------------
    std::vector<double> bsL = dcLeft_.At(i).GetSampleMeans();
    std::vector<double> bsR = dcRight_.At(i).GetSampleMeans();
    std::vector<double> bsMid = dcMid_.At(i).GetSampleMeans();
    std::vector<double> bsWeight = dcMid_.At(i).GetSampleWeights();
    for(int jSample=0; jSample<Nsamples; jSample++) {
      const float bsS = CalculateValue(bsL.at(jSample), bsR.at(jSample), bsMid.at(jSample), purity, shiftL, shiftR);
      const float bsW = bsWeight.at(jSample)*purity;
      out_container_->At(i).AddSampleMean(bsS);
      out_container_->At(i).AddSampleWeight(bsW);
    }
    //------------------------------------------
  }
}

float FitterSimple::CalculateValue(float vL, float vR, float vMid, float purity, float shiftL, float shiftR) const {
  const float value_bckgr = (vL*shiftR + vR*shiftL) / (shiftL + shiftR);
  const float value_sgnl = (vMid - (1-purity)*value_bckgr) / purity;

  return value_sgnl;
}

float FitterSimple::CalculateError(float errL, float errR, float errMid, float purity, float shiftL, float shiftR) const {
  const float err_bckgr2 = ( shiftR*shiftR*errL*errL + shiftL*shiftL*errR*errR ) /(shiftL+shiftR)/(shiftL+shiftR) ;
  const float err_sgnl2 = (errMid*errMid + (1-purity)*(1-purity)*err_sgnl2) / purity/purity;

  return std::sqrt(err_sgnl2);
}
