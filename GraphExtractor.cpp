#include "GraphExtractor.hpp"

void GraphExtractor::SetSelectAxis(std::string select_axis) {
  select_axis_ = select_axis;
  
//   number_of_select_axis_ = 0;
//   while(data_container_->GetAxes().at(number_of_select_axis_).Name() != select_axis_)
//     number_of_select_axis_++;
}

// TGraph* GraphExtractor::GetGraph(std::vector<int> bin)
// {
//   if(bin.size() != names_axes_to_exclude_.size())
//     throw std::runtime_error("Bin vector's dimensionality must equal to number of axes to exclude");
//     
//   std::vector<Qn::DataContainer<Qn::StatCalculate,Qn::Axis<double>>> data_container_reduced;
//   data_container_reduced.push_back(*data_container_);
//   
//   for(int j=0; j<bin.size(); j++)
//   {
//     double left_edge = data_container_ -> GetAxis(names_axes_to_exclude_.at(j).c_str()).GetLowerBinEdge(bin.at(j));
//     double right_edge = data_container_ -> GetAxis(names_axes_to_exclude_.at(j).c_str()).GetLowerBinEdge(bin.at(j)+1);
//     Qn::Axis<double> axis(names_axes_to_exclude_.at(j), 1, left_edge, right_edge);
//     data_container_reduced.push_back(data_container_reduced.back().Select(axis));
//   }
//   
//   TGraph* graph = Qn::DataContainerHelper::ToTGraph(data_container_reduced.back());
//   
//   return graph;
// }

TGraphErrors* GraphExtractor::GetGraph(std::vector<size_t> bin) const
{
  if(data_container_->GetAxes().size()-1 != bin.size())
    throw std::runtime_error("Bin vector's dimensionality must equal to DataContainer's dimensionality - 1");
    
  Qn::DataContainer<Qn::StatCalculate, Qn::Axis<double>> data_container_reduced;
  data_container_reduced = *data_container_;                                         // failure with memory exceed (std::bad_alloc)
  
  int j = 0;
  for(auto& ax : data_container_->GetAxes())
  {
    std::string ax_name = ax.Name();
    if(ax_name == select_axis_) continue;
    double left_edge = ax.GetLowerBinEdge(bin.at(j));
    double right_edge = ax.GetLowerBinEdge(bin.at(j)+1);
    Qn::Axis<double> axis(ax_name, 1, left_edge, right_edge);
    data_container_reduced = data_container_reduced.Select(axis);
    j++;
  }
  
  TGraphErrors* graph = Qn::DataContainerHelper::ToTGraph(data_container_reduced);
  
  return graph;
}

// TGraphErrors* GraphExtractor::GetGraph(std::vector<size_t> bin) const
// {
//   const auto& axis = data_container_->GetAxis(select_axis_.c_str());
//   TGraphErrors* graph = new TGraphErrors();
// 
//   for(int iselect=0; iselect<axis.size(); iselect++) {
//     std::vector<size_t> bin4D = bin;
//     bin4D.insert(bin4D.begin()+number_of_select_axis_, iselect);
//     Qn::StatCalculate scalc = data_container_->At(bin4D);
//     const double x = axis.GetBinCenter(iselect);
//     const double y = scalc.Mean();
//     const double ey = scalc.StandardErrorOfMean();
//     graph->SetPoint(iselect, x, y);
//     graph->SetPointError(iselect, 0, ey);
//     graph->GetXaxis()->SetTitle(select_axis_.c_str());
//   }
//   
//   return graph;
// }

std::vector<int> GraphExtractor::GetAxesSizes()
{
  std::vector<int> sizes;
  for(auto& axisname : names_axes_to_exclude_)
    sizes.push_back(data_container_->GetAxis(axisname.c_str()).GetNBins());
  
  return sizes;
}

std::vector<std::vector<double>> GraphExtractor::GetAxesBinEdges()
{
  std::vector<std::vector<double>> v_binedges;
  for(auto& axisname : names_axes_to_exclude_)
  {
    std::vector<double> binedges;
    for(int i=0; i<=data_container_->GetAxis(axisname.c_str()).GetNBins(); i++)
      binedges.push_back(data_container_->GetAxis(axisname.c_str()).GetLowerBinEdge(i));
    v_binedges.push_back(binedges);
  }
  
  return v_binedges;
}