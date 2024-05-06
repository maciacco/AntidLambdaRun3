#include "Constants.h"

void Systematics(){
  TFile* files[10];
  TH1D* rho[10];
  TH1D* k2k1[10];
  TH1D* sys[6];
  TFile out("outSys_newPtBins.root", "recreate");
  for (int i{0}; i < 10; ++i) {
    auto cutSet = cutSets[i];
    files[i] = TFile::Open(Form("%s_%s", cutSet, fileIO::outFileName.data()));
    rho[i] = (TH1D*)files[i]->Get("hAntidNetL");
    k2k1[i] = (TH1D*)files[i]->Get("hNetL_k2k1");
    rho[i]->SetName(Form("hAntidNetL%s", cutSet));
    k2k1[i]->SetName(Form("hNetL_k2k1%s", cutSet));
  }
  for (int i{0}; i < 5; ++i) {
    sys[i] = new TH1D(*rho[2 * i + 1]);
    if (i > 0 && i < 4) {
      for (int iB{1}; iB < 9; ++iB) {
        sys[i]->SetBinContent(iB, std::abs(rho[2 * i + 1]->GetBinContent(iB) - rho[2 * (i + 1)]->GetBinContent(iB)) * 0.5);
        sys[i]->SetBinError(iB, 0.);
      }
    }
    else {
      for (int iB{1}; iB < 9; ++iB) {
        sys[i]->SetBinContent(iB, std::abs(rho[2 * i + 1]->GetBinContent(iB) - rho[0]->GetBinContent(iB)));
        sys[i]->SetBinError(iB, 0.);
      }
    }
    out.cd();
    sys[i]->Write();
  }
  sys[5] = new TH1D("sysTot", ";Centrality (%);#rho", 10, 0, 100);
  for (int iB{1}; iB < 9; ++iB) {
    double sum = 0.;
    for (int i{0}; i < 5; ++i) {
      sum += std::pow(sys[i]->GetBinContent(iB), 2.);
    }
    sys[5]->SetBinContent(iB, std::sqrt(sum));
  }
  out.cd();
  sys[5]->Write();
}