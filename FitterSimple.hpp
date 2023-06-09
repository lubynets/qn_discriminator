#ifndef FitterSimple_H
#define FitterSimple_H

#include "DataContainer.hpp"

class FitterSimple {
public:

  FitterSimple() = default;
  virtual ~FitterSimple() = default;

  void SetDataContainer(Qn::DataContainerStatCalculate* datacontainer) { data_container_ = datacontainer; }
  void SetDataContainer(Qn::DataContainerStatCollect* datacontainer);
  void SetSelectAxis(std::string select_axis) { select_axis_ = select_axis; }
  void SetPurityContainer(Qn::DataContainerStatDiscriminator* puritycontainer) { purity_container_ = puritycontainer; }
  void SetEdgesOfLMR(std::array<float, 6> ar);
  void AutoSetEdgesOfLMR();

  void Calculate();

  Qn::DataContainerStatDiscriminator* GetResult() const { return out_container_; }

private:
  void DefineLRM();
  float CalculateValue(float vL, float vR, float vMid, float purity, float shiftL=1.f, float shiftR=1.f) const;
  float CalculateError(float errL, float errR, float errMid, float purity, float shiftL=1.f, float shiftR=1.f) const;

  Qn::DataContainerStatCalculate* data_container_{nullptr};
  Qn::DataContainerStatCalculate dcLeft_;
  Qn::DataContainerStatCalculate dcRight_;
  Qn::DataContainerStatCalculate dcMid_;

  Qn::DataContainerStatDiscriminator* purity_container_{nullptr};

  Qn::DataContainerStatDiscriminator* out_container_{nullptr};

  std::string select_axis_;
  std::pair<float, float> left_ranges_;
  std::pair<float, float> right_ranges_;
  std::pair<float, float> mid_ranges_;
};


#endif // FitterSimple_H
