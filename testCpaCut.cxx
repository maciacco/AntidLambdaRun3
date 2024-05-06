double cpaCut[]{0.99995, 0.9999, 0.9995, 0.999, 0.995, 0.99};
void testCpaCut() {
  TFile out("testCpaCutOut.root", "recreate");
  TFile cpa_("testYields_cpa_.root");
  TFile cpa0_("testYields_cpa0_.root");
  TFile cpa1_("testYields_cpa1_.root");
  TFile cpa2_("testYields_cpa2_.root");
  TFile cpa3_("testYields_cpa3_.root");
  TFile cpa4_("testYields_cpa4_.root");
  TH1F* h_cpa_[6];
  h_cpa_[0] = (TH1F*)cpa_.Get("k1");
  h_cpa_[1] = (TH1F*)cpa0_.Get("k1");
  h_cpa_[2] = (TH1F*)cpa1_.Get("k1");
  h_cpa_[3] = (TH1F*)cpa2_.Get("k1");
  h_cpa_[4] = (TH1F*)cpa3_.Get("k1");
  h_cpa_[5] = (TH1F*)cpa4_.Get("k1");
  TGraphErrors* hCpaCut[8];
  TF1 *ff[8];
  for (int i{0}; i < 8; ++i) {
    hCpaCut[i] = new TGraphErrors();
    hCpaCut[i]->SetName(Form("hCpa_%d", i));
    hCpaCut[i]->SetTitle(";cosPA;Yields");
    for (int iP{0}; iP < 6; ++iP) {
      hCpaCut[i]->AddPoint(cpaCut[iP], h_cpa_[iP]->GetBinContent(i + 1));
      hCpaCut[i]->SetPointError(hCpaCut[i]->GetN() - 1, 1.e-6, h_cpa_[iP]->GetBinError(i + 1));
    }
    ff[i] = new TF1(Form("ff_%d", i), "[0]-[0]/(1. + TMath::Exp([1]*x - [2]))");
    ff[i]->SetParLimits(0, hCpaCut[i]->GetPointY(0) * 0.9, hCpaCut[i]->GetPointY(0) * 1.1);
    ff[i]->SetParLimits(2, 0., 1.e7);
    ff[i]->SetParLimits(1, -1.e7, 0.);
    hCpaCut[i]->Fit(ff[i]);
    out.cd();
    hCpaCut[i]->Write();
  }
  out.Close();
}