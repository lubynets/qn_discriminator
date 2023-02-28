#ifndef QnToolsHelper_H
#define QnToolsHelper_H

#include "DataContainer.hpp"

using namespace Qn;

inline std::vector<TGraph*> SamplesToTGraph(const DataContainer<StatCalculate, AxisD> &data) {
  if (data.GetAxes().size() > 1) {
    std::cout << "Data container has more than one dimension. " << std::endl;
    std::cout << "Cannot draw as Graph. Use Projection() to make it one dimensional." << std::endl;
    throw;
  }
  const int Nsamples = data.At(0).GetSampleMeans().size();
  std::vector<TGraph*> graphs;
  graphs.resize(Nsamples);
  for(auto& gr : graphs) {
    gr = new TGraph();
  }
  unsigned int ibin = 0;
  for (const auto &bin : data) {
    if (bin.SumWeights() <= 0.) {
      ibin++;
      continue;
    }
    auto xhi = data.GetAxes().front().GetUpperBinEdge(ibin);
    auto xlo = data.GetAxes().front().GetLowerBinEdge(ibin);
    auto xhalfwidth = (xhi - xlo)/2.;
    auto x = xlo + xhalfwidth;

    std::vector<double> ys = bin.GetSampleMeans();
    for(int i_sample = 0; i_sample<Nsamples; i_sample++) {
      double y = ys.at(i_sample);
      graphs.at(i_sample)->SetPoint(graphs.at(i_sample)->GetN(), x, y);
      graphs.at(i_sample)->SetMarkerStyle(kFullCircle);
    }
    ibin++;
  }
  return graphs;
}

inline std::vector<double> SamplesWeights(const DataContainer<StatCalculate, AxisD> &data) {
  if (data.GetAxes().size() > 1) {
    std::cout << "Data container has more than one dimension. " << std::endl;
    std::cout << "Cannot draw as Graph. Use Projection() to make it one dimensional." << std::endl;
    throw;
  }
  const int Nsamples = data.At(0).GetSampleMeans().size();
  std::vector<double> weights;
  weights.resize(Nsamples);
  for(auto& we : weights) {
    we = 0.;
  }

  unsigned int ibin = 0;
  for (const auto &bin : data) {
    if (bin.SumWeights() <= 0.) {
      ibin++;
      continue;
    }
    auto binweights = bin.GetSampleWeights();
    for(int i_sample = 0; i_sample<Nsamples; i_sample++) {
      weights.at(i_sample) += binweights.at(i_sample);
    }
    ibin++;
  }

  return weights;
}

#endif // QnToolsHelper_H
