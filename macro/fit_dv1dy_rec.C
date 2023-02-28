void fit_dv1dy_rec() {
  std::string evegen = "dcmqgsm";
//   std::string evegen = "urqmd";

  std::string weightstatus = "wei";
//   std::string weightstatus = "now";

//   bool is_rebin_centrality = true;
  bool is_rebin_centrality = false;

  std::string fileName = "/home/oleksii/cbmdir/working/qna/rectrackspsi/" + evegen + "/cl.rtp." + evegen + "." + weightstatus + ".root";

  TFile* fileIn = TFile::Open(fileName.c_str());
  if(!fileIn) throw std::runtime_error("fileIn absent!");

  std::vector<std::string> particles{"lambda", "kshort"};
  std::vector<std::string> steps{"PLAIN", "RECENTERED", "TWIST", "RESCALED"};

  bool average_comp{true}; std::vector<std::string> components{"ave"};

//   bool average_comp{false}; std::vector<std::string> components{"x1x1", "y1y1"};

  std::string axistofit_sim = "SimParticles_rapidity";
  std::string axistofit_rec = "ReconstructedParticles_rapidity";
  const float midrapidity = 1.6217901;

  std::string fileOutName;
  if(!is_rebin_centrality) {
    fileOutName = "dv1dy.rtp." + evegen + "." + weightstatus + ".root";
  } else {
    fileOutName = "dv1dy.rebinned.rtp." + evegen + "." + weightstatus + ".root";
  }
  if(average_comp) {
    fileOutName.erase(fileOutName.length()-4, 4);
    fileOutName += "ave.root";
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
    }
    if(evegen == "dcmqgsm" && is_rebin_centrality) {
      if(pa == "lambda") pt_bin_edges = {0, 0.4, 1.0, 1.6};
      if(pa == "kshort") pt_bin_edges = {0, 0.4, 0.8, 1.6};
    }

    if(evegen == "urqmd") {
      if(pa == "lambda") pt_bin_edges = {0, 0.8, 1.2, 1.6};
      if(pa == "kshort") pt_bin_edges = {0, 0.8, 1.2, 1.6};
    }

    Qn::DataContainerStatCalculate v1sim_in = Qn::DataContainerStatCalculate(*((Qn::DataContainerStatCollect*) fileIn->Get(("sim/u_" + pa + "_sim_PLAIN.Q_psi_PLAIN.x1x1").c_str()))) +
                                              Qn::DataContainerStatCalculate(*((Qn::DataContainerStatCollect*) fileIn->Get(("sim/u_" + pa + "_sim_PLAIN.Q_psi_PLAIN.y1y1").c_str())));
    v1sim_in = v1sim_in/2.;

    Qn::DataContainerStatCalculate v1sim = v1sim_in.Rebin({"SimParticles_pT", pt_bin_edges});
    if(is_rebin_centrality) v1sim = v1sim.Rebin({"RecEventHeader_centrality_tracks", {0, 15, 40, 70}});

    const double fitaxis_lo = v1sim.GetAxis(axistofit_sim.c_str()).GetFirstBinEdge();
    const double fitaxis_hi = v1sim.GetAxis(axistofit_sim.c_str()).GetLastBinEdge();

    Qn::DataContainerStatCalculate v1sim_rebinned = v1sim.Rebin({axistofit_sim, 1, fitaxis_lo, fitaxis_hi});
    Qn::DataContainerStatCalculate v1sim_reduced = v1sim_rebinned.Select({axistofit_sim, 1, fitaxis_lo, fitaxis_hi});

    Qn::DataContainerStatDiscriminator v1sim_slope;
    Qn::DataContainerStatDiscriminator v1sim_intercept;

    v1sim_slope.AddAxes(v1sim_reduced.GetAxes());
    v1sim_intercept.AddAxes(v1sim_reduced.GetAxes());

    GraphExtractor gex_sim;
    gex_sim.SetDataContainer(&v1sim);
    gex_sim.SetSelectAxis(axistofit_sim.c_str());

    for (int i = 0; i < v1sim_reduced.size(); i++) {
      TGraphErrors* gr_sim = gex_sim.GetGraph(v1sim_reduced.GetIndex(i));

      TF1* fsim = new TF1("fsim", "[0]+[1]*(x-[2])", fitaxis_lo, fitaxis_hi);
      fsim->FixParameter(2, midrapidity);

      gr_sim->Fit(fsim, "0");

      v1sim_intercept[i].SetVEW(fsim->GetParameter(0), fsim->GetParError(0));
      v1sim_slope[i].SetVEW(fsim->GetParameter(1), fsim->GetParError(1));
    }
    v1sim_intercept.Write("v1sim_intercept.psi.ave");
    v1sim_slope.Write("v1sim_slope.psi.ave");

    for (auto& step : steps) {
      for(auto& co : components) {

        Qn::DataContainerStatCalculate v1_rec_in;
        if(average_comp) {
          v1_rec_in = Qn::DataContainerStatCalculate(*((Qn::DataContainerStatCollect*) fileIn->Get(("rec/u_" + pa + "_rec_" + step + ".Q_psi_PLAIN.x1x1").c_str()))) +
                      Qn::DataContainerStatCalculate(*((Qn::DataContainerStatCollect*) fileIn->Get(("rec/u_" + pa + "_rec_" + step + ".Q_psi_PLAIN.y1y1").c_str())));
          v1_rec_in = v1_rec_in/2.;
        } else {
          v1_rec_in = Qn::DataContainerStatCalculate(*((Qn::DataContainerStatCollect*) fileIn->Get(("rec/u_" + pa + "_rec_" + step + ".Q_psi_PLAIN." + co).c_str())));
        }

        Qn::DataContainerStatCalculate v1rec = v1_rec_in.Rebin({"ReconstructedParticles_pT", pt_bin_edges});
        if(is_rebin_centrality) v1rec = v1rec.Rebin({"RecEventHeader_centrality_tracks", {0, 15, 40, 70}});

        const double fitaxis_lo = v1rec.GetAxis(axistofit_rec.c_str()).GetFirstBinEdge();
        const double fitaxis_hi = v1rec.GetAxis(axistofit_rec.c_str()).GetLastBinEdge();

        Qn::DataContainerStatCalculate v1rec_rebinned = v1rec.Rebin({axistofit_rec, 1, fitaxis_lo, fitaxis_hi});
        Qn::DataContainerStatCalculate v1rec_reduced = v1rec_rebinned.Select({axistofit_rec, 1, fitaxis_lo, fitaxis_hi});

        Qn::DataContainerStatDiscriminator v1rec_slope;
        Qn::DataContainerStatDiscriminator v1rec_intercept;

        v1rec_slope.AddAxes(v1rec_reduced.GetAxes());
        v1rec_intercept.AddAxes(v1rec_reduced.GetAxes());

        GraphExtractor gex_rec;
        gex_rec.SetDataContainer(&v1rec);
        gex_rec.SetSelectAxis(axistofit_rec.c_str());

        for (int i = 0; i < v1sim_reduced.size(); i++) {
          TGraphErrors* gr_rec = gex_rec.GetGraph(v1rec_reduced.GetIndex(i));

          TF1* frec = new TF1("frec", "[0]+[1]*(x-[2])", fitaxis_lo, fitaxis_hi);
          frec->FixParameter(2, midrapidity);

//           gr_rec->Fit(frec, "0", "", midrapidity-0.01, fitaxis_hi);
          gr_rec->Fit(frec, "0", "");

          v1rec_intercept[i].SetVEW(frec->GetParameter(0), frec->GetParError(0));
          v1rec_slope[i].SetVEW(frec->GetParameter(1), frec->GetParError(1));
        }

        v1rec_intercept.Write(("v1rec_intercept." + step + "." + co).c_str());
        v1rec_slope.Write(("v1rec_slope." + step + "." + co).c_str());
      }
    }
  }
  fileOut->Close();
}
