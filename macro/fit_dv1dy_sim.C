void fit_dv1dy_sim() {
  std::string evegen = "dcmqgsm";
//   std::string evegen = "urqmd";

  bool is_rebin_centrality = true;
//   bool is_rebin_centrality = false;

  std::string fileName = "/home/oleksii/cbmdir/working/qna/simtracksflow/" + evegen + "/v1andR1.stf." + evegen + ".root";

  TFile* fileIn = TFile::Open(fileName.c_str());

  std::vector<std::string> particles{"lambda", "kshort", "xi", "pipos", "pineg"};
  std::vector<std::string> components;
  std::vector<std::string> subevents{"psd1", "psd2", "psd3",
                                     "etacut_1_charged", "etacut_2_charged", "etacut_3_charged",
                                     "etacut_1_all", "etacut_2_all", "etacut_3_all"};
  std::string step;
  bool average_comp;

  std::string axistofit = "SimParticles_rapidity";
  const float midrapidity = 1.6217901;

  std::string fileOutName;
  if(!is_rebin_centrality) {
    fileOutName = "dv1dy.stf." + evegen + ".root";
  } else {
    fileOutName = "dv1dy.rebinned.stf." + evegen + ".root";
  }
  TFile* fileOut = TFile::Open(fileOutName.c_str(), "recreate");
  fileOut->cd();
  for(auto& pa : particles) {
    fileOut->mkdir(pa.c_str());
    fileOut->cd(pa.c_str());

    std::vector<double> pt_bin_edges;

    if(evegen == "dcmqgsm" && !is_rebin_centrality) {
      if(pa == "lambda") pt_bin_edges = {0, 0.4, 0.8, 1.2, 1.6};
      if(pa == "kshort") pt_bin_edges = {0, 0.4, 0.8, 1.6};
      if(pa == "pipos" || pa == "pineg") pt_bin_edges = {0, 0.4, 0.6, 1.0, 1.4, 2.0};
      if(pa == "xi") pt_bin_edges = {0, 0.4, 0.8, 1.2, 1.6};
    }
    if(evegen == "dcmqgsm" && is_rebin_centrality) {
      if(pa == "lambda") pt_bin_edges = {0, 0.4, 1.0, 1.6};
      if(pa == "kshort") pt_bin_edges = {0, 0.4, 0.8, 1.6};
      if(pa == "pipos" || pa == "pineg") pt_bin_edges = {0, 0.4, 0.8, 2.0};
      if(pa == "xi") pt_bin_edges = {0, 0.4, 1.0, 1.6};
    }

    if(evegen == "urqmd") {
      if(pa == "lambda") pt_bin_edges = {0, 0.8, 1.2, 1.6};
      if(pa == "kshort") pt_bin_edges = {0, 0.8, 1.2, 1.6};
      if(pa == "pipos" || pa == "pineg") pt_bin_edges = {0, 1.0, 1.4, 2.0};
      if(pa == "xi") pt_bin_edges = {0, 0.8, 1.2, 1.6};
    }

    Qn::DataContainerStatCalculate v1sim_in = *((Qn::DataContainerStatCalculate*) fileIn->Get(("v1/" + pa + "/uPsi/v1.uPsi.x1x1").c_str())) +
                                              *((Qn::DataContainerStatCalculate*) fileIn->Get(("v1/" + pa + "/uPsi/v1.uPsi.y1y1").c_str()));
    const int Nsamples = v1sim_in.At(0).GetSampleMeans().size();
    v1sim_in = v1sim_in/2.;
    Qn::DataContainerStatCalculate v1sim = v1sim_in.Rebin({"SimParticles_pT", pt_bin_edges});
    if(is_rebin_centrality) v1sim = v1sim.Rebin({"SimEventHeader_centrality_impactpar", {0, 15, 40, 70}});

    const double fitaxis_lo = v1sim.GetAxis(axistofit.c_str()).GetFirstBinEdge();
    const double fitaxis_hi = v1sim.GetAxis(axistofit.c_str()).GetLastBinEdge();

    Qn::DataContainerStatCalculate v1sim_rebinned = v1sim.Rebin({axistofit, 1, fitaxis_lo, fitaxis_hi});
    Qn::DataContainerStatCalculate v1sim_reduced = v1sim_rebinned.Select({axistofit, 1, fitaxis_lo, fitaxis_hi});

    Qn::DataContainerStatDiscriminator v1sim_slope;
    Qn::DataContainerStatDiscriminator v1sim_intercept;

    v1sim_slope.AddAxes(v1sim_reduced.GetAxes());
    v1sim_intercept.AddAxes(v1sim_reduced.GetAxes());

    GraphExtractor gex_sim;
    gex_sim.SetDataContainer(&v1sim);
    gex_sim.SetSelectAxis(axistofit.c_str());

    for (int i = 0; i < v1sim_reduced.size(); i++) {
      gex_sim.ReduceDataContainerToBin(v1sim_reduced.GetIndex(i));
      TGraphErrors* gr_sim = gex_sim.GetGraph();

      TF1* fsim = new TF1("fsim", "[0]+[1]*(x-[2])", fitaxis_lo, fitaxis_hi);
      fsim->FixParameter(2, midrapidity);

      gr_sim->Fit(fsim, "0");

      v1sim_intercept[i].SetVEW(fsim->GetParameter(0), fsim->GetParError(0));
      v1sim_slope[i].SetVEW(fsim->GetParameter(1), fsim->GetParError(1));
      delete fsim;

      std::vector<TGraph*> gr_sims = gex_sim.GetSamplesGraphs();
      std::vector<double> samples_weights = gex_sim.GetSamplesWeights();
      for(int isample = 0; isample<Nsamples; isample++) {
        fsim = new TF1("fsim", "[0]+[1]*(x-[2])", fitaxis_lo, fitaxis_hi);
        fsim->FixParameter(2, midrapidity);

        gr_sims.at(isample)->Fit(fsim, "0");

        v1sim_intercept[i].AddSampleMean(fsim->GetParameter(0));
        v1sim_slope[i].AddSampleMean(fsim->GetParameter(1));
        v1sim_intercept[i].AddSampleWeight(samples_weights.at(isample));
        v1sim_slope[i].AddSampleWeight(samples_weights.at(isample));
        delete fsim;
        delete gr_sims.at(isample);
      }
    }
    v1sim_intercept.Write("v1sim_intercept.psi.ave");
    v1sim_slope.Write("v1sim_slope.psi.ave");

    for (auto& se : subevents) {

      if(se[0] == 'p') {
        step = "_RECENTERED";
        average_comp = false;
        components = {"x1x1", "y1y1"};
      }
      if(se[0] == 'e') {
        step = "_PLAIN";
        average_comp = true;
        components = {"ave"};
      }

      for(auto& co : components) {

        Qn::DataContainerStatCalculate v1_rec_in;
        if(average_comp) {
          v1_rec_in = *((Qn::DataContainerStatCalculate*) fileIn->Get(("v1/" + pa + "/uQ_R1/v1.uQ_R1." + se + step +".x1x1").c_str())) +
                      *((Qn::DataContainerStatCalculate*) fileIn->Get(("v1/" + pa + "/uQ_R1/v1.uQ_R1." + se + step +".y1y1").c_str()));
          v1_rec_in = v1_rec_in/2.;
        } else {
          v1_rec_in = *((Qn::DataContainerStatCalculate*) fileIn->Get(("v1/" + pa + "/uQ_R1/v1.uQ_R1." + se + step +"." + co).c_str()));
        }

        Qn::DataContainerStatCalculate v1rec = v1_rec_in.Rebin({"SimParticles_pT", pt_bin_edges});
        if(is_rebin_centrality) v1rec = v1rec.Rebin({"SimEventHeader_centrality_impactpar", {0, 15, 40, 70}});

        const double fitaxis_lo = v1rec.GetAxis(axistofit.c_str()).GetFirstBinEdge();
        const double fitaxis_hi = v1rec.GetAxis(axistofit.c_str()).GetLastBinEdge();

        Qn::DataContainerStatCalculate v1rec_rebinned = v1rec.Rebin({axistofit, 1, fitaxis_lo, fitaxis_hi});
        Qn::DataContainerStatCalculate v1rec_reduced = v1rec_rebinned.Select({axistofit, 1, fitaxis_lo, fitaxis_hi});

        Qn::DataContainerStatDiscriminator v1rec_slope;
        Qn::DataContainerStatDiscriminator v1rec_intercept;

        v1rec_slope.AddAxes(v1rec_reduced.GetAxes());
        v1rec_intercept.AddAxes(v1rec_reduced.GetAxes());

        GraphExtractor gex_rec;
        gex_rec.SetDataContainer(&v1rec);
        gex_rec.SetSelectAxis(axistofit.c_str());

        for (int i = 0; i < v1sim_reduced.size(); i++) {
          gex_rec.ReduceDataContainerToBin(v1rec_reduced.GetIndex(i));
          TGraphErrors* gr_rec = gex_rec.GetGraph();

          TF1* frec = new TF1("frec", "[0]+[1]*(x-[2])", fitaxis_lo, fitaxis_hi);
          frec->FixParameter(2, midrapidity);

          gr_rec->Fit(frec, "0");

          v1rec_intercept[i].SetVEW(frec->GetParameter(0), frec->GetParError(0));
          v1rec_slope[i].SetVEW(frec->GetParameter(1), frec->GetParError(1));
          delete frec;

          std::vector<TGraph*> gr_recs = gex_rec.GetSamplesGraphs();
          std::vector<double> samples_weights = gex_rec.GetSamplesWeights();
          for(int isample = 0; isample<Nsamples; isample++) {
            frec = new TF1("frec", "[0]+[1]*(x-[2])", fitaxis_lo, fitaxis_hi);
            frec->FixParameter(2, midrapidity);

            gr_recs.at(isample)->Fit(frec, "0");

            v1rec_intercept[i].AddSampleMean(frec->GetParameter(0));
            v1rec_slope[i].AddSampleMean(frec->GetParameter(1));
            v1rec_intercept[i].AddSampleWeight(samples_weights.at(isample));
            v1rec_slope[i].AddSampleWeight(samples_weights.at(isample));
            delete frec;
            delete gr_recs.at(isample);
          }
        }

        v1rec_intercept.Write(("v1rec_intercept." + se + "." + co).c_str());
        v1rec_slope.Write(("v1rec_slope." + se + "." + co).c_str());
      }
    }
  }
  fileOut->Close();
}
