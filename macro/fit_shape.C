float GetChi2TH1FTF1(TH1F* histo, TF1* func);

#include "Helper.hpp"

void fit_shape() {

  gROOT->Macro( "/home/oleksii/cbmdir/flow_drawing_tools/example/style_1.cc" );

//   std::string fileName = "/home/oleksii/cbmdir/working/qna/shapes/massDC.apr20.dcmqgsm.12agev.recpid.lightcuts1.set4.3122.all.root";
//   const float mu = 1.115683;
//   const float sigma = 0.00145786;
//   std::string particle = "#Lambda";

  std::string fileName = "/home/oleksii/cbmdir/working/qna/shapes/massDC.apr20.dcmqgsm.12agev.recpid.lightcuts1.set4.310.all.root";
  const float mu = 0.497611;
  const float sigma = 0.0037;
  std::string particle = "K^{0}_{S}";

  TFile* fileIn = TFile::Open(fileName.c_str(), "read");

  Qn::DataContainer<TH1F, Qn::Axis<double>> dcprimary = *(fileIn->Get<Qn::DataContainer<TH1F, Qn::Axis<double>>>("dcmass"));
  auto dcIn = dcprimary.Rebin({"centrality", {0, 10, 20, 40, 70}});

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

  Qn::DataContainerShapeContainer dcReFit;
  dcReFit.AddAxes(dcIn.GetAxes());

  Qn::DataContainerShapeContainer dcPreFit;
  dcPreFit.AddAxes(dcIn.GetAxes());

  std::vector<Qn::DataContainerStatDiscriminator> dcParamsSgnl;
  dcParamsSgnl.resize(8);
  for (auto& dc : dcParamsSgnl) {
    dc.AddAxes(dcIn.GetAxes());
  }

  std::vector<Qn::DataContainerStatDiscriminator> dcParamsBckgr;
  dcParamsBckgr.resize(5);
  for (auto& dc : dcParamsBckgr) {
    dc.AddAxes(dcIn.GetAxes());
  }

  Qn::DataContainerStatDiscriminator dc_chi2_prefit;
  dc_chi2_prefit.AddAxes(dcIn.GetAxes());

  Qn::DataContainerStatDiscriminator dc_chi2_fit;
  dc_chi2_fit.AddAxes(dcIn.GetAxes());

  TFile* fileOut = TFile::Open("fileOut.root", "recreate");
  fileOut->mkdir("Fits");
  fileOut->mkdir("Params/sgnl");
  fileOut->mkdir("Params/bckgr");
  fileOut->mkdir("Chi2s");

  for (int i = 0; i < dcIn.size(); i++) {

    std::vector<size_t> indices = dcIn.GetIndex(i);
    std::string binname = "C" + StringBinNumber(indices.at(0) + 1) + "_pT" + StringBinNumber(indices.at(1) + 1) + "_y" + StringBinNumber(indices.at(2) + 1);
    const float C_lo = dcIn.GetAxis("centrality").GetLowerBinEdge(indices.at(0));
    const float C_hi = dcIn.GetAxis("centrality").GetUpperBinEdge(indices.at(0));
    const float pT_lo = dcIn.GetAxis("pT").GetLowerBinEdge(indices.at(1));
    const float pT_hi = dcIn.GetAxis("pT").GetUpperBinEdge(indices.at(1));
    const float y_lo = dcIn.GetAxis("y").GetLowerBinEdge(indices.at(2));
    const float y_hi = dcIn.GetAxis("y").GetUpperBinEdge(indices.at(2));

    std::cout << "DataContainer cell # " << i << "\n";
    std::cout << binname << "\n";
    std::cout << ("C: " + to_string_with_precision(C_lo, 2) + " - " + to_string_with_precision(C_hi, 2) + " %").c_str() << "\n";
    std::cout << ("p_{T}: " + to_string_with_precision(pT_lo, 2) + " - " + to_string_with_precision(pT_hi, 2) + " GeV/c").c_str() << "\n";
    std::cout << ("y_{LAB}: " + to_string_with_precision(y_lo, 2) + " - " + to_string_with_precision(y_hi, 2)).c_str() << "\n";

    ShapeFitter shFtr(&dcIn[i]);
    shFtr.SetExpectedMu(mu);
    shFtr.SetExpectedSigma(sigma);
    shFtr.Fit();

    dcReFit[i].SetInputHisto(&dcIn[i]);
    dcReFit[i].SetShape(shFtr.GetReFuncSgnl(), shFtr.GetReFuncBckgr());

    dcPreFit[i].SetInputHisto(&dcIn[i]);
    dcPreFit[i].SetShape(shFtr.GetFuncSgnl(), shFtr.GetFuncBckgr());

    fileOut->cd("Fits");
    TCanvas cc("", "", 1500, 900);
    cc.cd();

    dcIn[i].SetTitle(binname.c_str());
    dcIn[i].GetXaxis()->SetRangeUser(mu - 16 * sigma, mu + 16 * sigma);
    dcIn[i].SetLineWidth(2);
    dcIn[i].SetStats(kFALSE);
    dcIn[i].Draw("");

    shFtr.GetReGraphAll()->SetFillStyle(3001);
    shFtr.GetReGraphAll()->SetFillColor(kRed - 4);
    shFtr.GetReGraphAll()->SetLineColor(kRed);
    shFtr.GetReGraphAll()->Draw("l e3 same");
    shFtr.GetReGraphBckgr()->SetLineColor(kGreen + 2);
    shFtr.GetReGraphBckgr()->Draw("l same");

    shFtr.GetGraphAll()->SetLineColor(kRed);
    shFtr.GetGraphAll()->SetLineStyle(2);
    shFtr.GetGraphAll()->Draw("l x same");
    shFtr.GetGraphBckgr()->SetLineColor(kGreen + 2);
    shFtr.GetGraphBckgr()->SetLineStyle(2);
    shFtr.GetGraphBckgr()->Draw("l x same");

    TLegend legend(0.12, 0.52, 0.27, 0.73);
    legend.SetBorderSize(0);
    legend.AddEntry(shFtr.GetReGraphAll(), "All Fit", "L");
    legend.AddEntry(shFtr.GetReGraphBckgr(), "BG Fit", "L");
    legend.AddEntry(shFtr.GetGraphAll(), "All Pre-Fit", "L");
    legend.AddEntry(shFtr.GetGraphBckgr(), "BG Pre-Fit", "L");
    legend.SetTextSize(0.03);
    legend.SetTextFont(22);
    legend.Draw("same");

    TPaveText binedges(0.12, 0.75, 0.27, 0.92, "brNDC");
    binedges.AddText(particle.c_str());
    binedges.AddText(("C: " + to_string_with_precision(C_lo, 2) + " - " + to_string_with_precision(C_hi, 2) + " %").c_str());
    binedges.AddText(("p_{T}: " + to_string_with_precision(pT_lo, 2) + " - " + to_string_with_precision(pT_hi, 2) + " GeV/c").c_str());
    binedges.AddText(("y_{LAB}: " + to_string_with_precision(y_lo, 2) + " - " + to_string_with_precision(y_hi, 2)).c_str());
    binedges.SetFillColor(0);
    binedges.SetTextSize(0.03);
    binedges.SetTextFont(22);
    binedges.Draw("same");

    TF1* func_sgnl = shFtr.GetReFuncSgnl();
    TF1* func_bckgr = shFtr.GetReFuncBckgr();
    const int npar_sgnl = func_sgnl->GetNpar();
    const int npar_bckgr = func_bckgr->GetNpar();
    std::vector<float> par_sgnl, parerr_sgnl;
    par_sgnl.resize(npar_sgnl);
    parerr_sgnl.resize(npar_sgnl);
    std::vector<float> par_bckgr, parerr_bckgr;
    par_bckgr.resize(npar_bckgr);
    parerr_bckgr.resize(npar_bckgr);
    for (int j = 0; j < npar_sgnl; j++) {
      par_sgnl.at(j) = func_sgnl->GetParameter(j);
      parerr_sgnl.at(j) = func_sgnl->GetParError(j);
      dcParamsSgnl.at(j)[i].SetVEW(par_sgnl.at(j), parerr_sgnl.at(j));
    }
    for (int j = 0; j < npar_bckgr; j++) {
      par_bckgr.at(j) = func_bckgr->GetParameter(j);
      parerr_bckgr.at(j) = func_bckgr->GetParError(j);
      dcParamsBckgr.at(j)[i].SetVEW(par_bckgr.at(j), parerr_bckgr.at(j));
    }

    TPaveText ptpar(0.77, 0.51, 0.95, 0.92, "brNDC");
    ptpar.AddText("DSCB parameters");
    ptpar.AddText(("Height = " + to_string_with_precision(par_sgnl.at(0), 2) + " #pm " + to_string_with_precision(parerr_sgnl.at(0), 2)).c_str());
//     ptpar.AddText(("#mu_{ref} = " + to_string_with_precision(mu, 4)).c_str());
    ptpar.AddText(("#mu - #mu_{ref} = (" + to_string_with_precision(par_sgnl.at(2) * 1e4, 3) + " #pm " + to_string_with_precision(parerr_sgnl.at(2) * 1e4, 3) + ") #times 10^{-4}").c_str());
    ptpar.AddText(("#sigma = (" + to_string_with_precision(par_sgnl.at(3) * 1e3, 3) + " #pm " + to_string_with_precision(parerr_sgnl.at(3) * 1e3, 3) + ") #times 10^{-3}").c_str());
    ptpar.AddText(("a_{1} = " + to_string_with_precision(par_sgnl.at(4), 2) + " #pm " + to_string_with_precision(parerr_sgnl.at(4), 2)).c_str());
    ptpar.AddText(("lg n_{1} = " + to_string_with_precision(par_sgnl.at(5), 2) + " #pm " + to_string_with_precision(parerr_sgnl.at(5), 2)).c_str());
    ptpar.AddText(("a_{2} = " + to_string_with_precision(par_sgnl.at(6), 2) + " #pm " + to_string_with_precision(parerr_sgnl.at(6), 2)).c_str());
    ptpar.AddText(("lg n_{2} = " + to_string_with_precision(par_sgnl.at(7), 2) + " #pm " + to_string_with_precision(parerr_sgnl.at(7), 2)).c_str());
    ptpar.AddText("pol3 parameters");
    ptpar.AddText(("p_{0} = " + to_string_with_precision(par_bckgr.at(0), 2) + " #pm " + to_string_with_precision(parerr_bckgr.at(0), 2)).c_str());
    ptpar.AddText(("p_{1} = " + to_string_with_precision(par_bckgr.at(1), 2) + " #pm " + to_string_with_precision(parerr_bckgr.at(1), 2)).c_str());
    ptpar.AddText(("p_{2} = " + to_string_with_precision(par_bckgr.at(2), 2) + " #pm " + to_string_with_precision(parerr_bckgr.at(2), 2)).c_str());
    ptpar.AddText(("p_{3} = " + to_string_with_precision(par_bckgr.at(3), 2) + " #pm " + to_string_with_precision(parerr_bckgr.at(3), 2)).c_str());
    ptpar.SetFillColor(0);
    ptpar.SetTextSize(0.025);
    ptpar.SetTextFont(22);
    ptpar.Draw("same");

    const float chi2_fit = shFtr.GetChi2AllFit();
    const float chi2_prefit = GetChi2TH1FTF1(&dcIn[i], shFtr.GetFuncAll());
    TPaveText ptchi2(0.30, 0.79, 0.45, 0.89, "brNDC");
    ptchi2.AddText(("#chi^{2}/_{ndf} Fit = " + to_string_with_precision(chi2_fit, 2)).c_str());
    ptchi2.AddText(("#chi^{2}/_{ndf} Pre-Fit = " + to_string_with_precision(chi2_prefit, 2)).c_str());
    ptchi2.SetFillColor(0);
    ptchi2.SetTextSize(0.03);
    ptchi2.SetTextFont(22);
    ptchi2.Draw("same");

    const float yield_sgnl = dcReFit[i].GetSignalIntegral(mu - 3 * sigma, mu + 3 * sigma);
    const float yield_bckgr = dcReFit[i].GetBackgroundIntegral(mu - 3 * sigma, mu + 3 * sigma);
    TPaveText ptyield(0.12, 0.40, 0.27, 0.50, "brNDC");
    ptyield.AddText(("N_{sgnl}^{*} = " + to_string_with_precision(yield_sgnl, 0)).c_str());
    ptyield.AddText(("N_{bckgr} = " + to_string_with_precision(yield_bckgr, 0)).c_str());
    ptyield.AddText(("S/B = " + to_string_with_precision(yield_sgnl / yield_bckgr, 2)).c_str());
    ptyield.SetFillColor(0);
    ptyield.SetTextSize(0.03);
    ptyield.SetTextFont(22);
    ptyield.Draw("same");

    cc.Write(binname.c_str());

    dc_chi2_prefit[i].SetVEW(chi2_prefit);
    dc_chi2_fit[i].SetVEW(chi2_fit);

    std::cout << "\n\n";
  }
  //
  fileOut->cd();
  dcReFit.Write("ReFit");
  dcPreFit.Write("PreFit");
  std::vector<std::string> SgnlParamNames = {"Height", "mu", "mu_shift", "sigma", "a1", "lgn1", "a2", "lgn2"};
  std::vector<std::string> BckgrParamNames = {"p0", "p1", "p2", "p3"};

  fileOut->cd("Params/sgnl");
  for (int j = 0; j < SgnlParamNames.size(); j++) {
    dcParamsSgnl.at(j).Write(SgnlParamNames.at(j).c_str());
  }

  fileOut->cd("Params/bckgr");
  for (int j = 0; j < BckgrParamNames.size(); j++) {
    dcParamsBckgr.at(j).Write(BckgrParamNames.at(j).c_str());
  }

  fileOut->cd("Chi2s");
  dc_chi2_prefit.Write("chi2_prefit");
  dc_chi2_fit.Write("chi2_fit");

  fileOut->Close();
  fileIn->Close();

  return 0;
}

float GetChi2TH1FTF1(TH1F* histo, TF1* func) {
  int firstbin = histo->FindBin(func->GetXmin());
  if (histo->GetBinCenter(firstbin) < func->GetXmin())
    firstbin++;

  int lastbin = histo->FindBin(func->GetXmax());
  if (histo->GetBinCenter(lastbin) > func->GetXmax())
    lastbin--;

  int ndf = 0;
  float chi2 = 0.f;
  for (int iBin = firstbin; iBin <= lastbin; iBin++) {
    if (histo->GetBinError(iBin) == 0.) continue;
    const float delta = (func->Eval(histo->GetBinCenter(iBin)) - histo->GetBinContent(iBin)) / histo->GetBinError(iBin);
    chi2 += delta * delta;
    ndf++;
  }

  ndf -= func->GetNumberFreeParameters();

//   std::cout << "chi2/ndf h-f = " << chi2 << " / " << ndf << "\n";

  return chi2 / ndf;
}
