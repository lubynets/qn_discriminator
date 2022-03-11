#include "ShapeFitter.hpp"

#include "DataContainer.hpp"

#include "TString.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TH3F.h"
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
  TString allfilename="/home/user/cbmdir/working/qna/shapes/massDC.apr20.dcmqgsm.nopid.lightcuts1.set4.all.root";
//   TString sgnlfilename="/home/user/cbmdir/working/qna/OLD/shapes/out.mass3D.apr20.dcmqgsm.nopid.lightcuts1.set4.sgnl_12.root";
//   TString bckgrfilename="/home/user/cbmdir/working/qna/OLD/shapes/out.mass3D.apr20.dcmqgsm.nopid.lightcuts1.set4.bckgr.root";
//   
  TFile* allfile = TFile::Open(allfilename, "read");
//   TFile* sgnlfile = TFile::Open(sgnlfilename, "read");
//   TFile* bckgrfile = TFile::Open(bckgrfilename, "read");
//   
  Qn::DataContainer<Qn::H1F, Qn::Axis<double>> dcall = *(allfile -> Get<Qn::DataContainer<Qn::H1F,Qn::Axis<double>>>("dcmass"));
//   TH1F* histoall = nullptr;
//   TH1F* histosgnl = nullptr;
//   TH1F* histobckgr = nullptr;
//   
//   const float mu = 1.115683;
//   const float sigma = 0.00145786;
//   
//   const int C_nbins = 3;                                                    // TODO rm hardcoded nbins
//   const int y_nbins = 4;
//   const int pT_nbins = 4;
//   
//   double C_edges_array[] = {0, 20, 40, 100};                                // TODO rm hardcoded binranges
//   double y_edges_array[] = {1.02179, 1.42179, 1.82179, 2.22179, 2.62179};
//   double pT_edges_array[] = {0.2, 0.5, 0.8, 1.1, 1.4};
//   
//   double* C_edges = &C_edges_array[0];
//   double* y_edges = &y_edges_array[0];
//   double* pT_edges = &pT_edges_array[0];
//   
  //***** definition of terms *************************
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
//   
//   TH3F hchi2_bckgr_rec_fit("hchi2_bckgr_rec_fit", "", C_nbins, C_edges, y_nbins, y_edges, pT_nbins, pT_edges);
//   SetAxesNames(&hchi2_bckgr_rec_fit);
//   
//   TH3F hchi2_bckgr_mc_and_bckgr_rec_fit("hchi2_bckgr_mc_and_bckgr_rec_fit", "", C_nbins, C_edges, y_nbins, y_edges, pT_nbins, pT_edges);
//   SetAxesNames(&hchi2_bckgr_mc_and_bckgr_rec_fit);
//   
//   TH3F hchi2_sgnl_mc_and_sgnl_rec("hchi2_sgnl_mc_and_sgnl_rec", "", C_nbins, C_edges, y_nbins, y_edges, pT_nbins, pT_edges);
//   SetAxesNames(&hchi2_sgnl_mc_and_sgnl_rec);
//     
//   TH3F hchi2_sgnl_mc_fit("hchi2_sgnl_mc_fit", "", C_nbins, C_edges, y_nbins, y_edges, pT_nbins, pT_edges);
//   SetAxesNames(&hchi2_sgnl_mc_fit);
//   
//   TH3F hchi2_sgnl_rec_fit("hchi2_sgnl_rec_fit", "", C_nbins, C_edges, y_nbins, y_edges, pT_nbins, pT_edges);
//   SetAxesNames(&hchi2_sgnl_rec_fit);
//   
//   TH3F hchi2_sgnl_mc_and_sgnl_rec_fit("hchi2_sgnl_mc_and_sgnl_rec_fit", "", C_nbins, C_edges, y_nbins, y_edges, pT_nbins, pT_edges);
//   SetAxesNames(&hchi2_sgnl_mc_and_sgnl_rec_fit);
//   
//   TH3F hchi2_all_and_all_fit("hchi2_all_and_all_fit", "", C_nbins, C_edges, y_nbins, y_edges, pT_nbins, pT_edges);
//   SetAxesNames(&hchi2_all_and_all_fit);
//   
//   TH3F hchi2_all_fit("hchi2_all_fit", "", C_nbins, C_edges, y_nbins, y_edges, pT_nbins, pT_edges);
//   SetAxesNames(&hchi2_all_fit);
  Qn::DataContainerStatDiscriminator dc_chi2_all_and_all_fit;
  dc_chi2_all_and_all_fit.AddAxes(dcall.GetAxes());
  
  Qn::DataContainerStatDiscriminator dc_chi2_all_fit;
  dc_chi2_all_fit.AddAxes(dcall.GetAxes());
  
//   
//   TH3F hchi2_bckgr_mc_and_bckgr_refit("hchi2_bckgr_mc_and_bckgr_refit", "", C_nbins, C_edges, y_nbins, y_edges, pT_nbins, pT_edges);
//   SetAxesNames(&hchi2_bckgr_mc_and_bckgr_refit);
//   
//   TH3F hchi2_sgnl_mc_and_sgnl_refit("hchi2_sgnl_mc_and_sgnl_refit", "", C_nbins, C_edges, y_nbins, y_edges, pT_nbins, pT_edges);
//   SetAxesNames(&hchi2_sgnl_mc_and_sgnl_refit);
//   
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
    ShapeFitter sftr((TH1F*) &dcall[i]);
    sftr.Fit();
    sct[i].SetShape(sftr.GetReFuncSgnl(), sftr.GetReFuncBckgr());
    sct[i].SetChi2BckgrFit(sftr.GetChi2BckgrFit());
    
    std::vector indices = dcall.GetIndex(i);
    std::string binname = "C" + StringBinNumber(indices.at(0)+1) + "_y" + StringBinNumber(indices.at(1)+1) + "_pT" + StringBinNumber(indices.at(2)+1);
    fileOut -> cd("all_and_all_fit");
    TCanvas c6("", "", 1500, 900);
    c6.cd();
    dcall[i].Draw();
    sftr.GetGraphAll()->SetFillStyle(3001);
    sftr.GetGraphAll()->SetFillColor(kRed-4);
    sftr.GetGraphAll()->SetLineColor(kRed);
    sftr.GetGraphAll()->SetLineWidth(2);
    sftr.GetGraphAll() -> Draw("l e3 same");
    c6.Write(binname.c_str());
  
  
//   for(int iC=0; iC<C_nbins; iC++)
//     for(int iy=0; iy<y_nbins; iy++)
//       for(int ipT=0; ipT<pT_nbins; ipT++)
//       {
//         std::string binname = "C" + StringBinNumber(iC+1) + "_y" + StringBinNumber(iy+1) + "_pT" + StringBinNumber(ipT+1);
//         histoall = (TH1F*) allfile -> Get(binname.c_str());
//         histosgnl = (TH1F*) sgnlfile -> Get(binname.c_str());
//         histobckgr = (TH1F*) bckgrfile -> Get(binname.c_str());
//         ShapeFitter sftr(histoall);
//         sftr.Fit();
// //         sct.SetShape(histosgnl, histobckgr, {iC, iy, ipT});
// //         sct.SetShape(sftr.GetReHistoSgnl(), sftr.GetReHistoBckgr(), {iC, iy, ipT});                     // Problems with visualization of fitting - function is 0
//         sct.SetShape(sftr.GetReFuncSgnl(), sftr.GetReFuncBckgr(), {iC, iy, ipT});                     // Problems with visualization of fitting
//         sct.SetChi2BckgrFit(sftr.GetChi2BckgrFit(), {iC, iy, ipT});                                     // Why do I need it?
//                 
//         fileOut -> cd("bckgr_mc_and_bckgr_rec_fit");
//         TCanvas c1("", "", 1500, 900);
//         c1.cd();
//         histobckgr -> Draw();
//         sftr.GetGraphBckgr()->SetLineColor(kRed);
//         sftr.GetGraphBckgr()->SetLineWidth(2);
//         sftr.GetGraphBckgr()->SetFillStyle(3001);
//         sftr.GetGraphBckgr()->SetFillColor(kRed-4);
//         sftr.GetGraphBckgr() -> Draw("l e3 same");
//         c1.Write(binname.c_str());
//         
//         fileOut -> cd("sgnl_mc_and_sgnl_rec");
//         TCanvas c2("", "", 1500, 900);
//         c2.cd();
//         histosgnl -> Draw();
//         sftr.GetHistoSgnl() -> SetLineColor(kRed);
//         sftr.GetHistoSgnl() -> Draw("same");
//         c2.Write(binname.c_str());
//         
//         fileOut -> cd("sgnl_mc_and_sgnl_mc_fit");
//         TCanvas c3("", "", 1500, 900);
//         c3.cd();
//         histosgnl -> Draw();
//         ShapeFitter sftr_sgnl_mc(histosgnl);  // can be initialized with any histogram, not important
//         sftr_sgnl_mc.DefineSgnlFunc(histosgnl, mu - 15*sigma, mu + 15*sigma);
//         sftr_sgnl_mc.FitSgnl(histosgnl);
//         TGraphErrors* gr_sgnl_fit = sftr_sgnl_mc.FuncWithErrors({sftr_sgnl_mc.GetFuncSgnl(), sftr_sgnl_mc.GetCovSgnl()});
//         gr_sgnl_fit->SetLineColor(kRed);
//         gr_sgnl_fit->SetLineWidth(2);
//         gr_sgnl_fit->SetFillStyle(3001);
//         gr_sgnl_fit->SetFillColor(kRed-4);
//         gr_sgnl_fit -> Draw("l e3 same");
//         c3.Write(binname.c_str());
//         
//         fileOut -> cd("sgnl_rec_and_sgnl_rec_fit");
//         TCanvas c4("", "", 1500, 900);
//         c4.cd();
//         sftr.GetHistoSgnl() -> SetLineColor(kBlue);
//         sftr.GetHistoSgnl() -> Draw();
//         sftr.GetGraphSgnl()->SetFillStyle(3001);
//         sftr.GetGraphSgnl()->SetFillColor(kRed-4);
//         sftr.GetGraphSgnl()->SetLineColor(kRed);
//         sftr.GetGraphSgnl()->SetLineWidth(2);
//         sftr.GetGraphSgnl() -> Draw("l e3 same");
//         c4.Write(binname.c_str());
//         
//         fileOut -> cd("sgnl_mc_and_sgnl_rec_fit");
//         TCanvas c5("", "", 1500, 900);
//         c5.cd();
//         histosgnl -> Draw();
//         sftr.GetGraphSgnl()->SetFillStyle(3001);
//         sftr.GetGraphSgnl()->SetFillColor(kRed-4);
//         sftr.GetGraphSgnl()->SetLineColor(kRed);
//         sftr.GetGraphSgnl()->SetLineWidth(2);
//         sftr.GetGraphSgnl() -> Draw("l e3 same");
//         c5.Write(binname.c_str());
//         
//         fileOut -> cd("all_and_all_fit");
//         TCanvas c6("", "", 1500, 900);
//         c6.cd();
//         histoall -> Draw();
//         sftr.GetGraphAll()->SetFillStyle(3001);
//         sftr.GetGraphAll()->SetFillColor(kRed-4);
//         sftr.GetGraphAll()->SetLineColor(kRed);
//         sftr.GetGraphAll()->SetLineWidth(2);
//         sftr.GetGraphAll() -> Draw("l e3 same");
//         c6.Write(binname.c_str());
//         
//         fileOut -> cd("bckgr_mc_and_bckgr_refit");
//         TCanvas c7("", "", 1500, 900);
//         c7.cd();
//         histobckgr -> Draw();
// //         sftr.GetReHistoBckgr() -> Draw("same");
//         sftr.GetReGraphBckgr()->SetFillStyle(3001);
//         sftr.GetReGraphBckgr()->SetFillColor(kRed-4);
//         sftr.GetReGraphBckgr()->SetLineColor(kRed);
//         sftr.GetReGraphBckgr()->SetLineWidth(2);
//         sftr.GetReGraphBckgr() -> Draw("l e3 same");
//         c7.Write(binname.c_str());
//         
//         fileOut -> cd("sgnl_mc_and_sgnl_refit");
//         TCanvas c8("", "", 1500, 900);
//         c8.cd();
//         histosgnl -> Draw();
// //         sftr.GetReHistoSgnl() -> Draw("same");
//         sftr.GetReGraphSgnl()->SetFillStyle(3001);
//         sftr.GetReGraphSgnl()->SetFillColor(kRed-4);
//         sftr.GetReGraphSgnl()->SetLineColor(kRed);
//         sftr.GetReGraphSgnl()->SetLineWidth(2);
//         sftr.GetReGraphSgnl() -> Draw("l e3 same");
//         c8.Write(binname.c_str());
//                 
//         hchi2_bckgr_rec_fit.SetBinContent(iC+1, iy+1, ipT+1, sftr.GetChi2BckgrFit());
//         hchi2_bckgr_mc_and_bckgr_rec_fit.SetBinContent(iC+1, iy+1, ipT+1, GetChi2TH1FTF1(histobckgr, sftr.GetFuncBckgr()));
//         hchi2_sgnl_mc_and_sgnl_rec.SetBinContent(iC+1, iy+1, ipT+1, GetChi2TH1FTH1F(histosgnl, sftr.GetHistoSgnl()));
//         hchi2_sgnl_mc_fit.SetBinContent(iC+1, iy+1, ipT+1, sftr_sgnl_mc.GetFuncSgnl()->GetChisquare() / sftr_sgnl_mc.GetFuncSgnl()->GetNDF());
//         hchi2_sgnl_rec_fit.SetBinContent(iC+1, iy+1, ipT+1, sftr.GetChi2SgnlFit());
//         hchi2_sgnl_mc_and_sgnl_rec_fit.SetBinContent(iC+1, iy+1, ipT+1, GetChi2TH1FTF1(histosgnl, sftr.GetFuncSgnl()));
//     hchi2_all_and_all_fit.SetBinContent(indices.at(0)+1, indices.at(1)+1, indices.at(2)+1, GetChi2TH1FTF1(histoall, sftr.GetFuncAll()));
//     hchi2_all_fit.SetBinContent(indices.at(0)+1, indices.at(1)+1, indices.at(2)+1, sftr.GetChi2AllFit());
    dc_chi2_all_and_all_fit[i].SetVEW(GetChi2TH1FTF1((TH1F*) &dcall[i], sftr.GetFuncAll()));
    dc_chi2_all_fit[i].SetVEW(sftr.GetChi2AllFit());
//         hchi2_bckgr_mc_and_bckgr_refit.SetBinContent(iC+1, iy+1, ipT+1, GetChi2TH1FTF1(histobckgr, sftr.GetReFuncBckgr()));
//         hchi2_sgnl_mc_and_sgnl_refit.SetBinContent(iC+1, iy+1, ipT+1, GetChi2TH1FTF1(histosgnl, sftr.GetReFuncSgnl()));
  }
//   
  fileOut -> cd();
  sct.Write("dcshape");
//   hchi2_bckgr_rec_fit.Write();
//   hchi2_bckgr_mc_and_bckgr_rec_fit.Write();
//   hchi2_sgnl_mc_and_sgnl_rec.Write();
//   hchi2_sgnl_mc_fit.Write();
//   hchi2_sgnl_rec_fit.Write();
//   hchi2_sgnl_mc_and_sgnl_rec_fit.Write();
//   hchi2_all_and_all_fit.Write();
//   hchi2_all_fit.Write();
  dc_chi2_all_and_all_fit.Write("chi2_all_and_all_fit");
  dc_chi2_all_fit.Write("chi2_all_fit");
//   hchi2_bckgr_mc_and_bckgr_refit.Write();
//   hchi2_sgnl_mc_and_sgnl_refit.Write();
  fileOut -> Close();
//   
  allfile -> Close();
//   sgnlfile -> Close();
//   bckgrfile -> Close();
  
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