#include <TFile.h>
#include <THnSparse.h>
#include <TH1D.h>

#include "Constants.h"

void subsample(THnSparse *nPart, TH1D *hOut){
  TH1D* nPartVsCent = nullptr;
  auto nSubsamples = nPart->GetAxis(0)->GetNbins();
  auto nCent = nPart->GetAxis(1)->GetNbins();
  for (int iC{1}; iC < nCent + 1; ++iC) {
    for (int iS{0}; iS < nSubsamples; ++iS) {
      auto binc = nPart->GetBinContent((int[]){iS + 1, iC, 4});
      auto binc_prev = hOut->GetBinContent(iC);
      auto bine_prev = hOut->GetBinError(iC);
      std::cout << binc << std::endl;
      hOut->SetBinContent(iC, binc_prev + binc);
      hOut->SetBinError(iC, bine_prev + binc * binc);
    }
    hOut->SetBinContent(iC, hOut->GetBinContent(iC) / nSubsamples);
    auto stdDev = hOut->GetBinError(iC) / nSubsamples - hOut->GetBinContent(iC) * hOut->GetBinContent(iC);
    stdDev = std::sqrt(stdDev / (nSubsamples - 1));
    hOut->SetBinError(iC, stdDev);
  }
}

void setAxisRanges(THnSparse* h, std::initializer_list<int> const xmin, std::initializer_list<int> const xmax){
  int size = xmin.size();
  if (xmax.size() != size) {
    std::cout << "fatal: size inconsistency" << std::endl;
    return;
  }
  auto imin = xmin.begin();
  auto imax = xmax.begin();
  int iA{0};
  for (; imin != xmin.end(); ++imin, ++imax, ++iA) {
    // std::cout << *imin << " - " << *imax << std::endl;
    h->GetAxis(iA)->SetRange(*imin, *imax);
  }
}

void RawCorrelation(){
  auto f = TFile::Open(fileIO::inFileName.data());
  auto of = TFile::Open(fileIO::outFileName.data(), "recreate");
  auto nEv = dynamic_cast<THnSparse*>(f->Get("antid-lambda-ebye/nEv"));
  auto nAntid = dynamic_cast<THnSparse*>(f->Get("antid-lambda-ebye/nAntid"));
  auto nAntiL = dynamic_cast<THnSparse*>(f->Get("antid-lambda-ebye/nAntiL"));
  auto nL = dynamic_cast<THnSparse*>(f->Get("antid-lambda-ebye/nL"));
  auto nSqAntiL = dynamic_cast<THnSparse*>(f->Get("antid-lambda-ebye/nSqAntiL"));
  auto nSqL = dynamic_cast<THnSparse*>(f->Get("antid-lambda-ebye/nSqL"));
  auto nLAntiL = dynamic_cast<THnSparse*>(f->Get("antid-lambda-ebye/nLantiL"));

  auto nSubsamples = nAntid->GetAxis(0)->GetNbins();
  TH1D *hAntiL = new TH1D("hAntiL", ";Centrality (%);#kappa_{2}/#kappa_{1}(#bar{#Lambda})", 100, 0., 100.);
  TH1D *hAntiL_k1_check = new TH1D("hAntiL_k1_check", ";Centrality (%);#kappa_{1}(#bar{#Lambda})", 100, 0., 100.);
  TH1D *hAntiL_k2_check = new TH1D("hAntiL_k2_check", ";Centrality (%);#kappa_{2}(#bar{#Lambda})", 100, 0., 100.);
  // TH1D *hAntiL1 = new TH1D("hAntiL1", ";Centrality (%);#kappa_{1}(#bar{#Lambda})", 100, 0., 100.);
  TH1D *hL = new TH1D("hL", ";Centrality (%);#kappa_{2}/#kappa_{1}(#Lambda)", 100, 0., 100.);
  TH1D *hNetL = new TH1D("hNetL", ";Centrality (%);#kappa_{2}/#kappa_{1}(#Delta#Lambda)", 100, 0., 100.);

  int bins[3]{5, 100, 4};
  double xmin[3]{0., 0., 0.};
  double xmax[3]{5., 100., 0.8};
  THnSparseD hk2k1("hk2k1", ";Subsample;Centrality (%);#Delta#eta;", 3, bins, xmin, xmax);
  THnSparseD hk1_check("hk1_check", ";Subsample;Centrality (%);#Delta#eta;", 3, bins, xmin, xmax);
  THnSparseD hk2_check("hk2_check", ";Subsample;Centrality (%);#Delta#eta;", 3, bins, xmin, xmax);
  //TH3D hk2k13("hk2k13", ";Subsample;Centrality (%);#Delta#eta;", bins[0], xmin[0], xmax[0], bins[1], xmin[1], xmax[1], bins[2], xmin[2], xmax[2]);

  TH1D *proj_nev = nullptr, *proj_antiL = nullptr;
  TH2D *proj_sqAntiL = nullptr;
  for (int iS{1}; iS < bins[0] + 1; ++iS) {
    setAxisRanges(nEv, {iS}, {iS});
    proj_nev = dynamic_cast<TH1D*>(nEv->Projection(1));
    std::cout << " --- S U B S A M P L E    N . " << iS << " --- \n";
    for (int iC{1}; iC < bins[1] + 1; ++iC) {
      auto nev = proj_nev->GetBinContent(iC);
      std::cout << "Processing centrality " << iC << "...\n";
      for (int iE{1}; iE < bins[2] + 1; ++iE) {

        setAxisRanges(nAntiL, {iS, iC, iE}, {iS, iC, iE});
        setAxisRanges(nL, {iS, iC, iE}, {iS, iC, iE});
        setAxisRanges(nSqAntiL, {iS, iC, iE}, {iS, iC, iE});
        setAxisRanges(nSqL, {iS, iC, iE}, {iS, iC, iE});
        setAxisRanges(nLAntiL, {iS, iC, iE}, {iS, iC, iE});

        proj_antiL = dynamic_cast<TH1D*>(nAntiL->Projection(3));
        proj_sqAntiL = dynamic_cast<TH2D*>(nSqAntiL->Projection(3, 4));

        double sum = 0.;
        double sumSq = 0.;
        for (int iP{1}; iP < proj_antiL->GetNbinsX() + 1; ++iP) { // loop over pT bins
          double n = proj_antiL->GetBinContent(iP);
          for (int iP2{1}; iP2 < proj_antiL->GetNbinsX() + 1; ++iP2) {
            double nsq = proj_sqAntiL->GetBinContent(iP, iP2);
            sumSq += nsq;
          }
          sum += n;
        }
        if (nev > 1.e-5) {
          auto k1 = sum / nev;
          auto k2 = sumSq / nev - std::pow(k1, 2.);
          // std::cout << k2 << std::endl;
          hk2_check.SetBinContent((int[]){iS, iC, iE}, k2);
          hk1_check.SetBinContent((int[]){iS, iC, iE}, k1);
          if (k1 > 0) {
            hk2k1.SetBinContent((int[]){iS, iC, iE}, k2 / k1);
          }
          else {
            hk2k1.SetBinContent((int[]){iS, iC, iE}, 0.);
          }
        }
        else {
          hk2k1.SetBinContent((int[]){iS, iC, iE}, 0.);
          hk1_check.SetBinContent((int[]){iS, iC, iE}, 0.);
          hk2_check.SetBinContent((int[]){iS, iC, iE}, 0.);
        }
        //hk2k13.SetBinContent(iS, iC, iE, k1);
        proj_antiL->Clear();
        proj_sqAntiL->Clear();
      }
    }
    proj_nev->Clear();
  }

  subsample(&hk2k1, hAntiL);
  subsample(&hk1_check, hAntiL_k1_check);
  subsample(&hk2_check, hAntiL_k2_check);
  // subsample1(&hk2k1, hAntiL1);

  of->cd();
  hAntiL->Write();
  hAntiL_k1_check->Write();
  hAntiL_k2_check->Write();
  //hAntiL1->Write();
  //hk2k1.Write();
  //hk2k13.Write();
  // hL->Write();

  of->Close();
}
