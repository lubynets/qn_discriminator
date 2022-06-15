#include "Helper.hpp"
#include "ShapeFitter.hpp"

#include "DataContainer.hpp"

#include "TString.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TH3F.h"
#include "TPaveText.h"
// 
// #include <string>

std::string StringBinNumber(int number);
float GetChi2TH1FTF1(TH1F* histo, TF1* func);
float GetChi2TH1FTH1F(TH1F* histo1, TH1F* histo2);
void SetAxesNames(TH3F* histo,
                  TString xaxisname="centrality, %",
                  TString yaxisname="rapidity",
                  TString zaxisname="p_{T}, GeV");

int main(int argc, char** argv)
{
//   TString allfilename="/home/user/cbmdir/working/qna/shapes/massDC.apr20.dcmqgsm.nopid.lightcuts1.set4.all.root";
//   TString allfilename="/home/user/cbmdir/working/qna/shapes/massDC.apr20.dcmqgsm.12agev.recpid.lightcuts1.set4.310.all.root";
  TString allfilename = "/home/user/cbmdir/working/qna/shapes/massDC.apr20.dcmqgsm.3.3agev.recpid.lightcuts1.set4.3122.all.root";

  TFile* allfile = TFile::Open(allfilename, "read");

//   Qn::DataContainer<TH1F, Qn::Axis<double>> dcall = *(allfile -> Get<Qn::DataContainer<TH1F,Qn::Axis<double>>>("dcmass"));

  Qn::DataContainer<TH1F, Qn::Axis<double>> dcprimary = *(allfile -> Get<Qn::DataContainer<TH1F,Qn::Axis<double>>>("dcmass"));
  auto dcall = dcprimary.Rebin({"centrality", {0,20,40,70}});

  const float mu = 1.115683;
  const float sigma = 0.00145786;
  
//   const float mu = 0.497;
//   const float sigma = 0.0037;
  
  

  //
  // all              shape of sgnl&bckgr together (the only distribution we will operate with in real life)
//   // sgnl_mc          shape of sgnl (MC-true)
//   // bckgr_mc         shape of bckgr (MC-bckgr)
//   // sgnl_mc_fit      fit of sgnl_mc
//   // bckgr_rec        exclude peak region from all
//   // bckgr_rec_fit    fit of bckgr_rec
//   // sgnl_rec         all minus bckgr_fit
//   // sgnl_rec_fit     fit of sgnl_rec
//   // all_fit          fit of all with sum of s&b
//   // bckgr_refit      bckgr part of all_fit
//   // sgnl_refit       sgnl part of all_fit
//   // We consider 2 types of chi2's: chi2 of fit and chi2 of difference between 2 histograms / functions / histogram and function
//   //
  //***************************************************

  Qn::DataContainerStatDiscriminator dc_chi2_all_and_all_fit;
  dc_chi2_all_and_all_fit.AddAxes(dcall.GetAxes());
  
  Qn::DataContainerStatDiscriminator dc_chi2_all_fit;
  dc_chi2_all_fit.AddAxes(dcall.GetAxes());
  
  Qn::DataContainerShapeContainer sct;
  sct.AddAxes(dcall.GetAxes());
  
//   ShapeContainerTensor sct;
//   sct.SetFrame({C_nbins, y_nbins, pT_nbins});
//   
  TFile* fileOut = TFile::Open("shapetensor_fit.root", "recreate");
//   fileOut -> mkdir("bckgr_mc_and_bckgr_rec_fit"); // is chi2 between (hchi2_bckgr_mc_and_bckgr_rec_fit)
//   fileOut -> mkdir("sgnl_mc_and_sgnl_rec");       // is chi2 between (hchi2_sgnl_mc_and_sgnl_rec)
//   fileOut -> mkdir("sgnl_mc_and_sgnl_mc_fit");    // is chi2 of fit  (hchi2_sgnl_mc_fit)
//   fileOut -> mkdir("sgnl_rec_and_sgnl_rec_fit");  // is chi2 of fit  (hchi2_sgnl_rec_fit)
//   fileOut -> mkdir("sgnl_mc_and_sgnl_rec_fit");   // is chi2 between (hchi2_sgnl_mc_and_sgnl_rec_fit)
  fileOut -> mkdir("all_and_all_fit");            // is chi2 of fit  (hchi2_all_fit)
//   fileOut -> mkdir("bckgr_mc_and_bckgr_refit");   // is chi2 between (hchi2_bckgr_mc_and_bckgr_refit)
//   fileOut -> mkdir("sgnl_mc_and_sgnl_refit");     // is chi2 between (hchi2_sgnl_mc_and_sgnl_refit)
//   
  for(int i=0; i<dcall.size(); i++) {
    ShapeFitter sftr(&dcall[i]);
    sftr.SetExpectedMu(mu);
    sftr.SetExpectedSigma(sigma);
    sftr.Fit();
    sct[i].SetShape(sftr.GetReFuncSgnl(), sftr.GetReFuncBckgr());
    sct[i].SetBinWidth(dcall[i].GetBinWidth(0));
    sct[i].SetChi2BckgrFit(sftr.GetChi2BckgrFit());
    
    std::vector indices = dcall.GetIndex(i);
    std::string binname = "C" + StringBinNumber(indices.at(0)+1) + "_y" + StringBinNumber(indices.at(1)+1) + "_pT" + StringBinNumber(indices.at(2)+1);
    const float C_lo = dcall.GetAxis("centrality").GetLowerBinEdge(indices.at(0));
    const float C_hi = dcall.GetAxis("centrality").GetUpperBinEdge(indices.at(0));
    const float pT_lo = dcall.GetAxis("pT").GetLowerBinEdge(indices.at(1));
    const float pT_hi = dcall.GetAxis("pT").GetUpperBinEdge(indices.at(1));
    const float y_lo = dcall.GetAxis("y").GetLowerBinEdge(indices.at(2));
    const float y_hi = dcall.GetAxis("y").GetUpperBinEdge(indices.at(2));
    
    fileOut -> cd("all_and_all_fit");
    TCanvas c6("", "", 1500, 900);
    c6.cd();
    dcall[i].SetTitle(binname.c_str());
    dcall[i].Draw();
    sftr.GetGraphAll()->SetFillStyle(3001);
    sftr.GetGraphAll()->SetFillColor(kRed-4);
    sftr.GetGraphAll()->SetLineColor(kRed);
    sftr.GetGraphAll()->SetLineWidth(2);
    sftr.GetGraphAll() -> Draw("l e3 same");
    sftr.GetReGraphBckgr()->SetFillStyle(3001);
    sftr.GetReGraphBckgr()->SetFillColor(kGreen-4);
    sftr.GetReGraphBckgr()->SetLineColor(kGreen+2);
    sftr.GetReGraphBckgr()->SetLineWidth(2);
    sftr.GetReGraphBckgr() -> Draw("l e3 same");
    
    TPaveText binedges(0.15, 0.70, 0.30, 0.85, "brNDC");
    binedges.AddText(("C: " + to_string_with_precision(C_lo, 2) + " - " + to_string_with_precision(C_hi, 2) + " %").c_str());
    binedges.AddText(("p_{T}: " + to_string_with_precision(pT_lo, 2) + " - " + to_string_with_precision(pT_hi, 2) + " GeV/c").c_str());
    binedges.AddText(("y_{LAB}: " + to_string_with_precision(y_lo, 2) + " - " + to_string_with_precision(y_hi, 2)).c_str());
    binedges.SetFillColor(0);
    binedges.SetTextSize(0.03);
    binedges.SetTextFont(22);
    binedges.Draw("same");
    
    c6.Write(binname.c_str());

    dc_chi2_all_and_all_fit[i].SetVEW(GetChi2TH1FTF1(&dcall[i], sftr.GetFuncAll()));
    dc_chi2_all_fit[i].SetVEW(sftr.GetChi2AllFit());
  }
//   
  fileOut -> cd();
  sct.Write("dcshape");
  dc_chi2_all_and_all_fit.Write("chi2_all_and_all_fit");
  dc_chi2_all_fit.Write("chi2_all_fit");
  fileOut -> Close();
  allfile -> Close();
  
  return 0;
}

float GetChi2TH1FTF1(TH1F* histo, TF1* func)
{
  int firstbin = histo -> FindBin(func->GetXmin());
  if(histo->GetBinCenter(firstbin) < func->GetXmin())
    firstbin++;
  
  int lastbin = histo -> FindBin(func->GetXmax());
  if(histo->GetBinCenter(lastbin) > func->GetXmax())
    lastbin--;
  
  int ndf = 0;
  float chi2 = 0.f;
  for(int iBin=firstbin; iBin<=lastbin; iBin++)
  {
    if(histo->GetBinError(iBin) == 0.) continue;
    const float delta = (func->Eval(histo->GetBinCenter(iBin)) - histo->GetBinContent(iBin)) / histo->GetBinError(iBin);
    chi2 += delta*delta;
    ndf++;
  }
  
  ndf -= func->GetNumberFreeParameters();
  
//   std::cout << "# of free parameters = " << func->GetNumberFreeParameters() << "\n";
  
  std::cout << "chi2/ndf h-f = " << chi2 << " / " << ndf << "\n";
  
  return chi2/ndf;
}

float GetChi2TH1FTH1F(TH1F* histo1, TH1F* histo2)
{
  int ndf = 0;
  float chi2 = 0.f;
  for(int iBin=1; iBin<=histo1->GetNbinsX(); iBin++)
  {
    if(histo1->GetBinError(iBin) == 0. || histo2->GetBinError(iBin) == 0.) continue;
    const float v1 = histo1->GetBinContent(iBin);
    const float v2 = histo2->GetBinContent(iBin);
    const float e1 = histo1->GetBinError(iBin);
    const float e2 = histo2->GetBinError(iBin);
    chi2 += (v1-v2)*(v1-v2)/(e1*e1 + e2*e2);
    ndf++;
  }
  std::cout << "chi2/ndf h-h = " << chi2 << " / " << ndf << "\n";
  
  return chi2/ndf;
}

std::string StringBinNumber(int number)
{
  if(number<10)
    return "0" + std::to_string(number);
  else
    return std::to_string(number);
}

void SetAxesNames(TH3F* histo, TString xaxisname, TString yaxisname, TString zaxisname)
{
  histo -> GetXaxis() -> SetTitle(xaxisname);
  histo -> GetYaxis() -> SetTitle(yaxisname);
  histo -> GetZaxis() -> SetTitle(zaxisname);
}