#ifndef GraphExtractor_H
#define GraphExtractor_H

#include "DataContainer.hpp"

#include "TGraphErrors.h"

#include <vector>
#include <string>

class GraphExtractor
{
public:
  
  GraphExtractor() = default;
  virtual ~GraphExtractor() = default;
  
  void SetDataContainer(Qn::DataContainer<Qn::StatCalculate,Qn::Axis<double>>* datacontainer) { data_container_ = (Qn::DataContainer<Qn::StatCalculate,Qn::Axis<double>>*)datacontainer; };
  void SetNamesAxesToExclude(std::vector<std::string> names_axes_to_exclude) { names_axes_to_exclude_ = names_axes_to_exclude; };
  void SetSelectAxis(std::string select_axis);
  
  TGraphErrors* GetGraph(std::vector<size_t> bin) const;
  std::vector<int> GetAxesSizes();
  std::vector<std::vector<double>> GetAxesBinEdges();
  
private:
   
  Qn::DataContainer<Qn::StatCalculate,Qn::Axis<double>>* data_container_{nullptr};
  std::vector<std::string> names_axes_to_exclude_;
  std::string select_axis_;
//   int number_of_select_axis_;
};

#endif  //GraphExtractor_H