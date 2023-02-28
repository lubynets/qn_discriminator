#include "GraphExtractor.hpp"

void GraphExtractor::SetSelectAxis(std::string select_axis) {
  select_axis_ = select_axis;
}

TGraphErrors* GraphExtractor::GetGraph(std::vector<size_t> bin) const {
  if (data_container_->GetAxes().size() - 1 != bin.size())
    throw std::runtime_error("Bin vector's dimensionality must equal to DataContainer's dimensionality - 1");

  Qn::DataContainer<Qn::StatCalculate, Qn::Axis<double>> data_container_reduced;
  data_container_reduced = *data_container_;// failure with memory exceed (std::bad_alloc)

  int j = 0;
  for (auto& ax : data_container_->GetAxes()) {
    std::string ax_name = ax.Name();
    if (ax_name == select_axis_) continue;
    double left_edge = ax.GetLowerBinEdge(bin.at(j));
    double right_edge = ax.GetLowerBinEdge(bin.at(j) + 1);
    Qn::Axis<double> axis(ax_name, 1, left_edge, right_edge);
    data_container_reduced = data_container_reduced.Select(axis);
    j++;
  }

  TGraphErrors* graph = Qn::DataContainerHelper::ToTGraph(data_container_reduced);

  return graph;
}

std::vector<int> GraphExtractor::GetAxesSizes() {
  std::vector<int> sizes;
  for (auto& axisname : names_axes_to_exclude_)
    sizes.push_back(data_container_->GetAxis(axisname.c_str()).GetNBins());

  return sizes;
}

std::vector<std::vector<double>> GraphExtractor::GetAxesBinEdges() {
  std::vector<std::vector<double>> v_binedges;
  for (auto& axisname : names_axes_to_exclude_) {
    std::vector<double> binedges;
    for (int i = 0; i <= data_container_->GetAxis(axisname.c_str()).GetNBins(); i++)
      binedges.push_back(data_container_->GetAxis(axisname.c_str()).GetLowerBinEdge(i));
    v_binedges.push_back(binedges);
  }

  return v_binedges;
}

