#include "Helper.hpp"
#include "GraphExtractor.hpp"
#include "Fitter.hpp"

#include "ShapeContainer.hpp"

#include "TFile.h"
#include "TDirectory.h"
#include "TH3.h"
#include "TCanvas.h"
#include "TPaveText.h"

#include <iostream>


std::string StringBinNumber(int number);
void SetAxesNames(TH3F* histo,
                  TString xaxisname="centrality, %",
                  TString yaxisname="rapidity",
                  TString zaxisname="p_{T}, GeV");

int main(int argc, char** argv)
{
  gErrorIgnoreLevel = 6001;
  
//   TString shapefilename="/home/user/cbmdir/working/qna/shapes/shapetensor_fit.apr20.dcmqgsm.12agev.recpid.lightcuts1.set4.310.all.root";
//   TString shapefilename="/home/user/cbmdir/working/qna/shapes/shapetensor_fit.apr20.dcmqgsm.3.3agev.recpid.lightcuts1.set4.3122.all.root";
  TString shapefilename="/home/user/cbmdir/working/qna/shapes/shapetensor_fit.apr20.dcmqgsm.12agev.recpid.lightcuts1.set4.3122.all.root";
  TFile* shapefile = TFile::Open(shapefilename, "read");
  Qn::DataContainer<Qn::ShapeContainer, Qn::Axis<double>>* shcntr = (Qn::DataContainer<Qn::ShapeContainer, Qn::Axis<double>>*)shapefile -> Get("dcshape");
  
  TString v1filename="/home/user/cbmdir/working/qna/correlations/aXmass/v1andR1.dcmqgsm.3.3agev.apr20.recpid.lightcuts1.3122.set4.root";
//   TString v1filename="/home/user/cbmdir/working/qna/correlations/aXmass/v1andR1.dcmqgsm.apr20.recpid.lightcuts1.3122.set4.all.root";
  TFile* v1file = TFile::Open(v1filename, "read");
  
  std::string invmassaxis = "ReconstructedParticles_mass";
  
  std::vector<std::string> components = {"x1x1", "y1y1"};
  
  TFile* fileOut = TFile::Open("out.fitter.root", "recreate");
  TDirectory* dirFit = fileOut->mkdir("fit");
  TDirectory* dirPar = fileOut->mkdir("parameters");
  
  Qn::DataContainerStatDiscriminator dc_entries_sgnl;
  dc_entries_sgnl.AddAxes(shcntr->GetAxes());
  
  Qn::DataContainerStatDiscriminator dc_entries_bckgr;
  dc_entries_bckgr.AddAxes(shcntr->GetAxes());
  
  bool common_dc_for_all_components{false};
      
  for(auto& co : components) {
        
//     Qn::DataContainer<Qn::StatCalculate,Qn::Axis<double>> lambda_psi = 
//     Qn::DataContainerStatCalculate(*(Qn::DataContainer<Qn::StatCollect,Qn::Axis<double>>*)v1file -> Get(("rec/RESCALED/u_rec_RESCALED.Q_psi_PLAIN." + co).c_str()));
    
    
    Qn::DataContainer<Qn::StatCalculate,Qn::Axis<double>> lambda_psi_pre = 
    Qn::DataContainerStatCalculate(*(Qn::DataContainer<Qn::StatCalculate,Qn::Axis<double>>*)v1file -> Get(("v1/uQ_R1_sub4_sts_pipos/v1.u_rec_RESCALED.psd2_RECENTERED.res_sub4_sts_pipos." + co).c_str()));
    
    auto lambda_psi = lambda_psi_pre.Rebin({"AnaEventHeader_centrality_tracks", {0,20,40,70}});
    
    const double invmass_lo = lambda_psi.GetAxis(invmassaxis.c_str()).GetFirstBinEdge();
    const double invmass_hi = lambda_psi.GetAxis(invmassaxis.c_str()).GetLastBinEdge();
    
    Qn::DataContainer<Qn::StatCalculate,Qn::Axis<double>> lambda_psi_rebinned = lambda_psi.Rebin({invmassaxis, 1, invmass_lo, invmass_hi});
    Qn::DataContainer<Qn::StatCalculate,Qn::Axis<double>> lambda_psi_reduced = lambda_psi_rebinned.Select({invmassaxis, 1, invmass_lo, invmass_hi});
    
    Qn::DataContainerStatDiscriminator dc_signal;
    dc_signal.AddAxes(shcntr->GetAxes());
    
    Qn::DataContainerStatDiscriminator dc_bckgr_0;
    dc_bckgr_0.AddAxes(shcntr->GetAxes());
    
    Qn::DataContainerStatDiscriminator dc_bckgr_1;
    dc_bckgr_1.AddAxes(shcntr->GetAxes());
    
    Qn::DataContainerStatDiscriminator dc_fit_chi2ndf;
    dc_fit_chi2ndf.AddAxes(shcntr->GetAxes());
        
    GraphExtractor gex;
    gex.SetDataContainer(&lambda_psi);
    gex.SetSelectAxis(invmassaxis.c_str());
        
    dirFit->cd();
        
    for(int i=0; i<shcntr->size(); i++) {
      TGraphErrors* gr = gex.GetGraph(shcntr->GetIndex(i));
      std::vector indices = shcntr->GetIndex(i);
      std::string binname = "C" + StringBinNumber(indices.at(0)+1) + "_y" + StringBinNumber(indices.at(1)+1) + "_pT" + StringBinNumber(indices.at(2)+1) + "." + co;
      const float C_lo = shcntr->GetAxis("centrality").GetLowerBinEdge(indices.at(0));
      const float C_hi = shcntr->GetAxis("centrality").GetUpperBinEdge(indices.at(0));
      const float pT_lo = shcntr->GetAxis("pT").GetLowerBinEdge(indices.at(1));
      const float pT_hi = shcntr->GetAxis("pT").GetUpperBinEdge(indices.at(1));
      const float y_lo = shcntr->GetAxis("y").GetLowerBinEdge(indices.at(2));
      const float y_hi = shcntr->GetAxis("y").GetUpperBinEdge(indices.at(2));
      
      gr -> SetName(binname.c_str());
      gr -> GetXaxis() -> SetTitle("m_{inv}, GeV");
      
      const float sumweights = lambda_psi_reduced[i].SumWeights();
      const float sgnl_integral = shcntr->At(i).GetSignalIntegral(invmass_lo, invmass_hi);
      const float bckgr_integral = shcntr->At(i).GetBackgroundIntegral(invmass_lo, invmass_hi);
      
      const float sgnl_weight = sumweights * sgnl_integral / (sgnl_integral + bckgr_integral);
      const float bckgr_weight = sumweights * bckgr_integral / (sgnl_integral + bckgr_integral);
      
      Fitter fitter;
      fitter.SetShape(&shcntr->At(i));
      fitter.SetGraphToFit(gr);
      fitter.Fit();
      
      std::string vsignal = "v1 = " + std::to_string(fitter.GetVSignal()) + " #pm " + std::to_string(fitter.GetVSignalError());
              
      gr -> SetTitle(vsignal.c_str());
      TCanvas c1("", "", 1500, 900);
      c1.cd();
      gr -> Draw("AP");
      fitter.GetGraphFit()->SetFillStyle(3001);
      fitter.GetGraphFit()->SetFillColor(kRed-4);
      fitter.GetGraphFit()->SetLineColor(kRed);
      fitter.GetGraphFit()->SetLineWidth(2);
      fitter.GetGraphFit() -> Draw("l e3 same");
      
      TPaveText binedges(0.15, 0.70, 0.30, 0.85, "brNDC");
      binedges.AddText(("C: " + to_string_with_precision(C_lo, 2) + " - " + to_string_with_precision(C_hi, 2) + " %").c_str());
      binedges.AddText(("p_{T}: " + to_string_with_precision(pT_lo, 2) + " - " + to_string_with_precision(pT_hi, 2) + " GeV/c").c_str());
      binedges.AddText(("y_{LAB}: " + to_string_with_precision(y_lo, 2) + " - " + to_string_with_precision(y_hi, 2)).c_str());
      binedges.SetFillColor(0);
      binedges.SetTextSize(0.03);
      binedges.SetTextFont(22);
      binedges.Draw("same");      
      
      c1.Write(gr->GetName());
      
      dc_signal[i].SetValue(fitter.GetFitParameters().at(0));
      dc_signal[i].SetError(fitter.GetFitErrors().at(0));
      dc_signal[i].SetWeight(sgnl_weight);
      
      dc_bckgr_0[i].SetValue(fitter.GetFitParameters().at(1));
      dc_bckgr_0[i].SetError(fitter.GetFitErrors().at(1));
      dc_bckgr_0[i].SetWeight(bckgr_weight);
      
      dc_bckgr_1[i].SetValue(fitter.GetFitParameters().at(2));
      dc_bckgr_1[i].SetError(fitter.GetFitErrors().at(2));
      dc_bckgr_1[i].SetWeight(bckgr_weight);
      
      dc_fit_chi2ndf[i].SetVEW(fitter.GetFitChi2Ndf());
      
      if(common_dc_for_all_components) continue;
      
      dc_entries_sgnl[i].SetValue(shcntr->At(i).GetSignalIntegral(invmass_lo, invmass_hi));
      dc_entries_sgnl[i].SetError(std::sqrt(dc_entries_sgnl[i].Mean()));
      dc_entries_sgnl[i].SetWeight(dc_entries_sgnl[i].Mean());
      
      dc_entries_bckgr[i].SetValue(shcntr->At(i).GetBackgroundIntegral(invmass_lo, invmass_hi));
      dc_entries_bckgr[i].SetError(std::sqrt(dc_entries_bckgr[i].Mean()));
      dc_entries_bckgr[i].SetWeight(dc_entries_bckgr[i].Mean());
      
      common_dc_for_all_components = true;
      
    }
    
    dirPar->cd();
    
    dc_signal.Write(("signal." + co).c_str());
    dc_bckgr_0.Write(("bckgr_0." + co).c_str());
    dc_bckgr_1.Write(("bckgr_1." + co).c_str());
    dc_fit_chi2ndf.Write(("fit_chi2ndf." + co).c_str());
   
  }
  
  dc_entries_sgnl.Write("entries_sgnl");
  dc_entries_bckgr.Write("entries_bckgr");
    
  fileOut -> Close();
  
  return 0;
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