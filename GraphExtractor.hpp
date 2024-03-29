#ifndef GraphExtractor_H
#define GraphExtractor_H

#include "DataContainer.hpp"

#include "TGraphErrors.h"

#include <string>
#include <vector>

class GraphExtractor {
 public:
  GraphExtractor() = default;
  virtual ~GraphExtractor() = default;

  void SetDataContainer(Qn::DataContainer<Qn::StatDiscriminator, Qn::Axis<double>>* datacontainer) { data_container_ = datacontainer; };
  void SetDataContainer(Qn::DataContainer<Qn::StatCollect, Qn::Axis<double>>* datacontainer);
  void SetDataContainer(Qn::DataContainer<Qn::StatCalculate, Qn::Axis<double>>* datacontainer);
//   void SetNamesAxesToExclude(std::vector<std::string> names_axes_to_exclude) { names_axes_to_exclude_ = names_axes_to_exclude; };
  void SetSelectAxis(std::string select_axis);

  TGraphErrors* GetGraph(std::vector<size_t> bin);
  std::vector<TGraphErrors*> GetSamplesGraphs(std::vector<size_t> bin);
  std::vector<double> GetSamplesWeights(std::vector<size_t> bin);
  TGraphErrors* GetGraph() const;
  std::vector<TGraphErrors*> GetSamplesGraphs() const;
  std::vector<double> GetSamplesWeights() const;
  void ReduceDataContainerToBin(std::vector<size_t> bin);

 private:
  Qn::DataContainer<Qn::StatDiscriminator, Qn::Axis<double>> data_container_reduced_;
  Qn::DataContainer<Qn::StatDiscriminator, Qn::Axis<double>>* data_container_{nullptr};
  std::string select_axis_;
};

#endif//GraphExtractor_H
