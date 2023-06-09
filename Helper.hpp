#ifndef HELPER_H
#define HELPER_H

#include <iostream>
#include <sstream>
#include <string>

#include <TFile.h>

template<typename T>
inline std::string to_string_with_precision(const T a_value, const int n = 6) {
  std::ostringstream out;
  out.precision(n);
  out << std::fixed << a_value;
  return out.str();
}

inline std::string StringBinNumber(int number) {
  if (number < 10)
    return "0" + std::to_string(number);
  else
    return std::to_string(number);
}

inline void CD(TFile* file, std::string dir) {
  if(file == nullptr) throw std::runtime_error("Helper::CD() - file is nullptr");

  if(file->GetDirectory(dir.c_str()) == nullptr) file->mkdir(dir.c_str());
  file->cd(dir.c_str());
}

inline float HistoIntegral(TH1F* histo, float low, float up)
{
  if(up<low) {
    std::cout << "Warning: Helper::HistoIntegral: up<low\n";
    return -HistoIntegral(histo, up, low);
  }

  if(low<histo->GetXaxis()->GetBinLowEdge(0)) {
    std::cout << "Warning: Helper::HistoIntegral: low<histo->GetXaxis()\n";
    low = histo->GetXaxis()->GetBinLowEdge(0);
  }

  if(up>histo->GetXaxis()->GetBinLowEdge(histo->GetXaxis()->GetNbins()+1)) {
    std::cout << "Warning: Helper::HistoIntegral: up>histo->GetXaxis()->GetBinLowEdge(histo->GetXaxis()->GetNbins()+1)\n";
    up = histo->GetXaxis()->GetBinLowEdge(histo->GetXaxis()->GetNbins()+1);
  }

  const int lowbin = histo->GetXaxis()->FindBin(low);
  const int upbin = histo->GetXaxis()->FindBin(up);

  float integral = histo->Integral(lowbin, upbin);
  integral -= histo->GetBinContent(lowbin)*(low-histo->GetXaxis()->GetBinLowEdge(lowbin))/histo->GetXaxis()->GetBinWidth(lowbin);
  integral -= histo->GetBinContent(upbin)*(histo->GetXaxis()->GetBinLowEdge(upbin+1)-up)/histo->GetXaxis()->GetBinWidth(upbin);

  return integral;
}

#endif//HELPER_H
