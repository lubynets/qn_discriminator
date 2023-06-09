#include "GraphExtractor.hpp"

#include "QnToolsHelper.hpp"

void GraphExtractor::SetDataContainer(Qn::DataContainer<Qn::StatCollect, Qn::Axis<double>>* datacontainer) {
  data_container_ = new Qn::DataContainer<Qn::StatDiscriminator, Qn::Axis<double>>(*datacontainer);
}

void GraphExtractor::SetDataContainer(Qn::DataContainer<Qn::StatCalculate, Qn::Axis<double>>* datacontainer) {
  data_container_ = new Qn::DataContainer<Qn::StatDiscriminator, Qn::Axis<double>>(*datacontainer);
}

void GraphExtractor::SetSelectAxis(std::string select_axis) {
  select_axis_ = select_axis;
}

TGraphErrors* GraphExtractor::GetGraph(std::vector<size_t> bin) {
  ReduceDataContainerToBin(bin);
  return GetGraph();
}

std::vector<TGraphErrors*> GraphExtractor::GetSamplesGraphs(std::vector<size_t> bin) {
  ReduceDataContainerToBin(bin);
  return GetSamplesGraphs();
}

std::vector<double> GraphExtractor::GetSamplesWeights(std::vector<size_t> bin) {
  ReduceDataContainerToBin(bin);
  return GetSamplesWeights();
}

TGraphErrors* GraphExtractor::GetGraph() const {
  if(data_container_reduced_.size() == 0 || data_container_reduced_.GetDimension() != 1) {
    throw std::runtime_error("GraphExtractor::GetGraph() - data_container_reduced_ must be reduced using ReduceDataContainerToBin() function");
  }
  TGraphErrors* graph = Qn::DataContainerHelper::ToTGraph(data_container_reduced_);
  return graph;
}

std::vector<TGraphErrors*> GraphExtractor::GetSamplesGraphs() const {
  if(data_container_reduced_.size() == 0 || data_container_reduced_.GetDimension() != 1) {
    throw std::runtime_error("GraphExtractor::GetSamplesGraphs() - data_container_reduced_ must be reduced using ReduceDataContainerToBin() function");
  }
  std::vector<TGraphErrors*> graphs = SamplesToTGraph(data_container_reduced_);
  return graphs;
}

std::vector<double> GraphExtractor::GetSamplesWeights() const {
  if(data_container_reduced_.size() == 0 || data_container_reduced_.GetDimension() != 1) {
    throw std::runtime_error("GraphExtractor::GetSamplesWeights() - data_container_reduced_ must be reduced using ReduceDataContainerToBin() function");
  }
  std::vector<double> weights = SamplesWeights(data_container_reduced_);
  return weights;
}

void GraphExtractor::ReduceDataContainerToBin(std::vector<size_t> bin) {
  if (data_container_->GetAxes().size() - 1 != bin.size())
    throw std::runtime_error("Bin vector's dimensionality must equal to DataContainer's dimensionality - 1");

  data_container_reduced_ = *data_container_;// failure with memory exceed (std::bad_alloc)

  int j = 0;
  for (auto& ax : data_container_->GetAxes()) {
    std::string ax_name = ax.Name();
    if (ax_name == select_axis_) continue;
    double left_edge = ax.GetLowerBinEdge(bin.at(j));
    double right_edge = ax.GetLowerBinEdge(bin.at(j) + 1);
    Qn::Axis<double> axis(ax_name, 1, left_edge, right_edge);
    data_container_reduced_ = data_container_reduced_.Select(axis);
    j++;
  }
}
