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
  auto nAntip = dynamic_cast<THnSparse*>(f->Get("antid-lambda-ebye/nAntip"));
  auto nAntiL = dynamic_cast<THnSparse*>(f->Get("antid-lambda-ebye/nAntiL"));
  auto nL = dynamic_cast<THnSparse*>(f->Get("antid-lambda-ebye/nL"));
  auto nSqAntid = dynamic_cast<THnSparse*>(f->Get("antid-lambda-ebye/nSqAntid"));
  auto nSqAntip = dynamic_cast<THnSparse*>(f->Get("antid-lambda-ebye/nSqAntip"));
  auto nSqAntiL = dynamic_cast<THnSparse*>(f->Get("antid-lambda-ebye/nSqAntiL"));
  auto nSqL = dynamic_cast<THnSparse*>(f->Get("antid-lambda-ebye/nSqL"));
  auto nLAntiL = dynamic_cast<THnSparse*>(f->Get("antid-lambda-ebye/nLantiL"));
  auto nLAntid = dynamic_cast<THnSparse*>(f->Get("antid-lambda-ebye/nLantid"));
  auto nAntiLAntid = dynamic_cast<THnSparse*>(f->Get("antid-lambda-ebye/nAntiLantid"));
  auto nAntipAntid = dynamic_cast<THnSparse*>(f->Get("antid-lambda-ebye/nAntipAntid"));

  auto nSubsamples = nAntid->GetAxis(0)->GetNbins();

  // output 1d histograms
  TH1D *hAntid_k2k1 = new TH1D("hAntid_k2k1", ";Centrality (%);#kappa_{2}/#kappa_{1}(#bar{d})", bins::bins[1], bins::xmin[1], bins::xmax[1]);
  TH1D *hAntid_k2k1check = new TH1D("hAntid_k2k1check", ";Centrality (%);(#kappa_{2} - #kappa_{1})/#kappa_{1}^{2}(#bar{d})", bins::bins[1], bins::xmin[1], bins::xmax[1]);
  TH1D *hAntid_k1 = new TH1D("hAntid_k1", ";Centrality (%);#kappa_{1}(#bar{d})", bins::bins[1], bins::xmin[1], bins::xmax[1]);
  TH1D *hAntid_k2 = new TH1D("hAntid_k2", ";Centrality (%);#kappa_{2}(#bar{d})", bins::bins[1], bins::xmin[1], bins::xmax[1]);

  TH1D *hAntip_k2k1 = new TH1D("hAntip_k2k1", ";Centrality (%);#kappa_{2}/#kappa_{1}(#bar{p})", bins::bins[1], bins::xmin[1], bins::xmax[1]);
  TH1D *hAntip_k2k1check = new TH1D("hAntip_k2k1check", ";Centrality (%);(#kappa_{2} - #kappa_{1})/#kappa_{1}^{2}(#bar{p})", bins::bins[1], bins::xmin[1], bins::xmax[1]);
  TH1D *hAntip_k1 = new TH1D("hAntip_k1", ";Centrality (%);#kappa_{1}(#bar{p})", bins::bins[1], bins::xmin[1], bins::xmax[1]);
  TH1D *hAntip_k2 = new TH1D("hAntip_k2", ";Centrality (%);#kappa_{2}(#bar{p})", bins::bins[1], bins::xmin[1], bins::xmax[1]);

  TH1D *hL_k2k1 = new TH1D("hL_k2k1", ";Centrality (%);#kappa_{2}/#kappa_{1}(#Lambda)", bins::bins[1], bins::xmin[1], bins::xmax[1]);
  TH1D *hL_k2k1check = new TH1D("hL_k2k1check", ";Centrality (%);(#kappa_{2} - #kappa_{1})/#kappa_{1}^{2}(#Lambda)", bins::bins[1], bins::xmin[1], bins::xmax[1]);
  TH1D *hL_k1 = new TH1D("hL_k1", ";Centrality (%);#kappa_{1}(#Lambda)", bins::bins[1], bins::xmin[1], bins::xmax[1]);
  TH1D *hL_k2 = new TH1D("hL_k2", ";Centrality (%);#kappa_{2}(#Lambda)", bins::bins[1], bins::xmin[1], bins::xmax[1]);

  TH1D *hAntiL_k2k1 = new TH1D("hAntiL_k2k1", ";Centrality (%);#kappa_{2}/#kappa_{1}(#bar{#Lambda})", bins::bins[1], bins::xmin[1], bins::xmax[1]);
  TH1D *hAntiL_k2k1check = new TH1D("hAntiL_k2k1check", ";Centrality (%);(#kappa_{2} - #kappa_{1})/#kappa_{1}^{2}(#bar{#Lambda})", bins::bins[1], bins::xmin[1], bins::xmax[1]);
  TH1D *hAntiL_k1 = new TH1D("hAntiL_k1", ";Centrality (%);#kappa_{1}(#bar{#Lambda})", bins::bins[1], bins::xmin[1], bins::xmax[1]);
  TH1D *hAntiL_k2 = new TH1D("hAntiL_k2", ";Centrality (%);#kappa_{2}(#bar{#Lambda})", bins::bins[1], bins::xmin[1], bins::xmax[1]);

  TH1D *hNetL_k2k1 = new TH1D("hNetL_k2k1", ";Centrality (%);#kappa_{2}/#kappa_{1}(#Delta#Lambda)", bins::bins[1], bins::xmin[1], bins::xmax[1]);
  TH1D *hAntidNetL = new TH1D("hAntidNetL", ";Centrality (%);#rho_{#bar{d}#Delta#Lambda}", bins::bins[1], bins::xmin[1], bins::xmax[1]);
  TH1D *hLantid = new TH1D("hLantid", ";Centrality (%);#rho_{#bar{d}#Lambda}", bins::bins[1], bins::xmin[1], bins::xmax[1]);
  TH1D *hAntiLantid = new TH1D("hAntiLantid", ";Centrality (%);#rho_{#bar{d}#bar{#Lambda}}", bins::bins[1], bins::xmin[1], bins::xmax[1]);
  TH1D *hAntipAntid = new TH1D("hAntipAntid", ";Centrality (%);#rho_{#bar{d}#bar{p}}", bins::bins[1], bins::xmin[1], bins::xmax[1]);

  // temporary thn for calculations
  THnSparseD hnAntid_k2k1("hnAntid_k2k1", ";Subsample;Centrality (%);#Delta#eta;", 3, bins::bins, bins::xmin, bins::xmax);
  THnSparseD hnAntid_k2k1check("hnAntid_k2k1_check", ";Subsample;Centrality (%);#Delta#eta;", 3, bins::bins, bins::xmin, bins::xmax);
  THnSparseD hnAntid_k1("hnAntid_k1", ";Subsample;Centrality (%);#Delta#eta;", 3, bins::bins, bins::xmin, bins::xmax);
  THnSparseD hnAntid_k2("hnAntid_k2", ";Subsample;Centrality (%);#Delta#eta;", 3, bins::bins, bins::xmin, bins::xmax);

  THnSparseD hnAntip_k2k1("hnAntip_k2k1", ";Subsample;Centrality (%);#Delta#eta;", 3, bins::bins, bins::xmin, bins::xmax);
  THnSparseD hnAntip_k2k1check("hnAntip_k2k1_check", ";Subsample;Centrality (%);#Delta#eta;", 3, bins::bins, bins::xmin, bins::xmax);
  THnSparseD hnAntip_k1("hnAntip_k1", ";Subsample;Centrality (%);#Delta#eta;", 3, bins::bins, bins::xmin, bins::xmax);
  THnSparseD hnAntip_k2("hnAntip_k2", ";Subsample;Centrality (%);#Delta#eta;", 3, bins::bins, bins::xmin, bins::xmax);

  THnSparseD hnL_k2k1("hnL_k2k1", ";Subsample;Centrality (%);#Delta#eta;", 3, bins::bins, bins::xmin, bins::xmax);
  THnSparseD hnL_k2k1check("hnL_k2k1_check", ";Subsample;Centrality (%);#Delta#eta;", 3, bins::bins, bins::xmin, bins::xmax);
  THnSparseD hnL_k1("hnL_k1", ";Subsample;Centrality (%);#Delta#eta;", 3, bins::bins, bins::xmin, bins::xmax);
  THnSparseD hnL_k2("hnL_k2", ";Subsample;Centrality (%);#Delta#eta;", 3, bins::bins, bins::xmin, bins::xmax);

  THnSparseD hnAntiL_k2k1("hnAntiL_k2k1", ";Subsample;Centrality (%);#Delta#eta;", 3, bins::bins, bins::xmin, bins::xmax);
  THnSparseD hnAntiL_k2k1check("hnAntiL_k2k1_check", ";Subsample;Centrality (%);#Delta#eta;", 3, bins::bins, bins::xmin, bins::xmax);
  THnSparseD hnAntiL_k1("hnAntiL_k1", ";Subsample;Centrality (%);#Delta#eta;", 3, bins::bins, bins::xmin, bins::xmax);
  THnSparseD hnAntiL_k2("hnAntiL_k2", ";Subsample;Centrality (%);#Delta#eta;", 3, bins::bins, bins::xmin, bins::xmax);

  THnSparseD hnNetL("hNetL", ";Subsample;Centrality (%);#Delta#eta;", 3, bins::bins, bins::xmin, bins::xmax);
  THnSparseD hnLantid("hLantid", ";Subsample;Centrality (%);#Delta#eta;", 3, bins::bins, bins::xmin, bins::xmax);
  THnSparseD hnAntidNetL("hAntidNetL", ";Subsample;Centrality (%);#Delta#eta;", 3, bins::bins, bins::xmin, bins::xmax);
  THnSparseD hnAntiLantid("hAntiLantid", ";Subsample;Centrality (%);#Delta#eta;", 3, bins::bins, bins::xmin, bins::xmax);
  THnSparseD hnAntipAntid("hAntipAntid", ";Subsample;Centrality (%);#Delta#eta;", 3, bins::bins, bins::xmin, bins::xmax);

  // loop over bins
  TH1D *proj_nev = nullptr, *proj_antiL = nullptr, *proj_L = nullptr, *proj_antid = nullptr, *proj_antip = nullptr;
  TH2D *proj_sqAntiL = nullptr, *proj_sqL = nullptr, *proj_sqAntip = nullptr, *proj_sqAntid = nullptr;
  TH2D *proj_LantiL = nullptr, *proj_Lantid = nullptr, *proj_antiLantid = nullptr, *proj_antipAntid = nullptr;
  for (int iS{1}; iS < bins::bins[0] + 1; ++iS) {
    setAxisRanges(nEv, {iS}, {iS});
    proj_nev = dynamic_cast<TH1D*>(nEv->Projection(1));
    std::cout << " --- S U B S A M P L E    N . " << iS << " --- \n";
    for (int iC{1}; iC < bins::bins[1] + 1; ++iC) {
      int iCsmallMin = nEv->GetAxis(1)->FindBin(hAntid_k1->GetXaxis()->GetBinLowEdge(iC) + constants::epsilon);
      int iCsmallMax = nEv->GetAxis(1)->FindBin(hAntid_k1->GetXaxis()->GetBinUpEdge(iC) - constants::epsilon);
      double nev_tot = proj_nev->Integral(iCsmallMin, iCsmallMax); // number of events in large centrality bin
      std::cout << "nev_tot = " << nev_tot << "\n";
      std::cout << "Processing centrality " << iC << "...\n";
      for (int iE{1}; iE < bins::bins[2] + 1; ++iE) {

        // Centrality bin width correction (CBWC)
        double k1_L = 0.;
        double k2_L = 0.;
        double k1_antiL = 0.;
        double k2_antiL = 0.;
        double k1_antid = 0.;
        double k2_antid = 0.;
        double k1_antip = 0.;
        double k2_antip = 0.;
        double k11_LantiL = 0.;
        double k11_Lantid = 0.;
        double k11_antiLantid = 0.;
        double k11_antipAntid = 0.;

        for (int iC_small{iCsmallMin}; iC_small < iCsmallMax + 1; ++iC_small) {
          double nev = proj_nev->GetBinContent(iC); // number of events in small centrality bin

          // select axis ranges
          setAxisRanges(nAntid, {iS, iC_small, iE}, {iS, iC_small, iE});
          setAxisRanges(nAntip, {iS, iC_small, iE}, {iS, iC_small, iE});
          setAxisRanges(nAntiL, {iS, iC_small, iE}, {iS, iC_small, iE});
          setAxisRanges(nL, {iS, iC_small, iE}, {iS, iC_small, iE});
          setAxisRanges(nSqAntid, {iS, iC_small, iE}, {iS, iC_small, iE});
          setAxisRanges(nSqAntip, {iS, iC_small, iE}, {iS, iC_small, iE});
          setAxisRanges(nSqAntiL, {iS, iC_small, iE}, {iS, iC_small, iE});
          setAxisRanges(nSqL, {iS, iC_small, iE}, {iS, iC_small, iE});
          setAxisRanges(nLAntiL, {iS, iC_small, iE}, {iS, iC_small, iE});
          setAxisRanges(nLAntid, {iS, iC_small, iE}, {iS, iC_small, iE});
          setAxisRanges(nAntiLAntid, {iS, iC_small, iE}, {iS, iC_small, iE});
          setAxisRanges(nAntipAntid, {iS, iC_small, iE}, {iS, iC_small, iE});

          double sum_antip = 0.;
          double sumSq_antip = 0.;
          double sum_antid = 0.;
          double sumSq_antid = 0.;
          double sum_antiL = 0.;
          double sumSq_antiL = 0.;
          double sum_L = 0.;
          double sumSq_L = 0.;
          double sum_LantiL = 0.;
          double sum_Lantid = 0.;
          double sum_antiLantid = 0.;
          double sum_antipAntid = 0.;

          // project on pT axes
          proj_antid = dynamic_cast<TH1D*>(nAntid->Projection(3));
          proj_sqAntid = dynamic_cast<TH2D*>(nSqAntid->Projection(3, 4));
          proj_antip = dynamic_cast<TH1D*>(nAntip->Projection(3));
          proj_sqAntip = dynamic_cast<TH2D*>(nSqAntip->Projection(3, 4));
          proj_antiL = dynamic_cast<TH1D*>(nAntiL->Projection(3));
          proj_sqAntiL = dynamic_cast<TH2D*>(nSqAntiL->Projection(3, 4));
          proj_L = dynamic_cast<TH1D*>(nL->Projection(3));
          proj_sqL = dynamic_cast<TH2D*>(nSqL->Projection(3, 4));
          proj_LantiL = dynamic_cast<TH2D*>(nLAntiL->Projection(3, 4));
          proj_Lantid = dynamic_cast<TH2D*>(nLAntid->Projection(3, 4));
          proj_antiLantid = dynamic_cast<TH2D*>(nAntiLAntid->Projection(3, 4));
          proj_antipAntid = dynamic_cast<TH2D*>(nAntipAntid->Projection(3, 4));

          for (int iPL{1}; iPL < proj_antiL->GetNbinsX() + 1; ++iPL) { // loop over pT bins (lambda)
            double n_L = proj_L->GetBinContent(iPL);
            double n_antiL = proj_antiL->GetBinContent(iPL);
            sum_L += n_L;
            sum_antiL += n_antiL;
            for (int iPL2{1}; iPL2 < proj_antiL->GetNbinsX() + 1; ++iPL2) { // loop over pT bins (lambda, 2)
              double nsq_L = proj_sqL->GetBinContent(iPL, iPL2);
              double nsq_antiL = proj_sqAntiL->GetBinContent(iPL, iPL2);
              double n_LantiL = proj_LantiL->GetBinContent(iPL, iPL2);
              sumSq_L += nsq_L;
              sumSq_antiL += nsq_antiL;
              sum_LantiL += n_LantiL;
            }
            for (int iPD{1}; iPD < proj_antid->GetNbinsX() + 1; ++iPD) { // loop over pT bins (deuteron)
              double n_Lantid = proj_Lantid->GetBinContent(iPL, iPD);
              double n_antiLantid = proj_antiLantid->GetBinContent(iPL, iPD);
              // std::cout << "n_antiLantid = " << n_antiLantid << std::endl;
              sum_Lantid += n_Lantid;
              sum_antiLantid += n_antiLantid;
            }
          }
          for (int iPD{1}; iPD < proj_antid->GetNbinsX() + 1; ++iPD) { // loop over pT bins (deuteron)
            double n_antid = proj_antid->GetBinContent(iPD);
            sum_antid += n_antid;
            for (int iPD2{1}; iPD2 < proj_antid->GetNbinsX() + 1; ++iPD2) { // loop over pT bins (deuteron, 2)
              double nsq_antid = proj_sqAntid->GetBinContent(iPD, iPD2);
              sumSq_antid += nsq_antid;
            }
          }
          for (int iPP{1}; iPP < proj_antid->GetNbinsX() + 1; ++iPP) { // loop over pT bins (proton)
            double n_antip = proj_antip->GetBinContent(iPP);
            sum_antip += n_antip;
            for (int iPP2{1}; iPP2 < proj_antip->GetNbinsX() + 1; ++iPP2) { // loop over pT bins (proton, 2)
              double nsq_antip = proj_sqAntip->GetBinContent(iPP, iPP2);
              sumSq_antip += nsq_antip;
            }
            for (int iPD2{1}; iPD2 < proj_antid->GetNbinsX() + 1; ++iPD2) { // loop over pT bins (deuteron, 2)
              double n_antipAntid = proj_antipAntid->GetBinContent(iPP, iPD2);
              sum_antipAntid += n_antipAntid;
            }
          }

          if (nev > 1.e-5) {
            k1_L += sum_L;
            k2_L += sumSq_L - std::pow(sum_L, 2.) / nev;
            k1_antiL += sum_antiL;
            k2_antiL += sumSq_antiL - std::pow(sum_antiL, 2.) / nev;
            k1_antid += sum_antid;
            k2_antid += sumSq_antid - std::pow(sum_antid, 2.) / nev;
            k1_antip += sum_antip;
            k2_antip += sumSq_antip - std::pow(sum_antip, 2.) / nev;
            k11_LantiL += sum_LantiL - sum_L * sum_antiL / nev;
            k11_Lantid += sum_Lantid - sum_L * sum_antid / nev;
            k11_antiLantid += sum_antiLantid - sum_antiL * sum_antid / nev;
            k11_antipAntid += sum_antipAntid - sum_antip * sum_antid / nev;

            // std::cout << "k1_L = " << k1_L << "\n";
          }

          proj_L->Clear();
          proj_sqL->Clear();
          proj_antiL->Clear();
          proj_sqAntiL->Clear();
          proj_antid->Clear();
          proj_sqAntid->Clear();
          proj_antip->Clear();
          proj_sqAntip->Clear();
          proj_LantiL->Clear();
          proj_Lantid->Clear();
          proj_antiLantid->Clear();
          proj_antipAntid->Clear();
        }

        if (nev_tot > 1.e-5) {
          k1_L /= nev_tot;
          k2_L /= nev_tot;
          k1_antiL /= nev_tot;
          k2_antiL /= nev_tot;
          k1_antid /= nev_tot;
          k2_antid /= nev_tot;
          k1_antip /= nev_tot;
          k2_antip /= nev_tot;
          k11_LantiL /= nev_tot;
          k11_Lantid /= nev_tot;
          k11_antiLantid /= nev_tot;
          k11_antipAntid /= nev_tot;

          std::cout << "k1_L = " << k1_L << ", k11_Lantid = " << k11_Lantid << ", k11_antiLantid = " << k11_antiLantid << "\n";

          // std::cout << "k11_antiLantid = " << k11_antiLantid << ", k2_antiL * k2_antid = " << k2_antiL * k2_antid << std::endl;
          // std::cout << "k11_Lantid = " << k11_Lantid << ", k2_L * k2_antid = " << k2_L * k2_antid << std::endl;
          // std::cout << "k11_antipAntid = " << k11_antipAntid << std::endl;

          // lambda
          hnL_k2.SetBinContent((int[]){iS, iC, iE}, k2_L);
          hnL_k1.SetBinContent((int[]){iS, iC, iE}, k1_L);
          if (k1_L > 0) {
            hnL_k2k1.SetBinContent((int[]){iS, iC, iE}, k2_L / k1_L);
            hnL_k2k1check.SetBinContent((int[]){iS, iC, iE}, (k2_L - k1_L) / std::pow(k1_L, 2.));
          }
          else {
            hnL_k2k1.SetBinContent((int[]){iS, iC, iE}, 0.);
            hnL_k2k1check.SetBinContent((int[]){iS, iC, iE}, 0.);
          }

          // anti-lambda
          hnAntiL_k2.SetBinContent((int[]){iS, iC, iE}, k2_antiL);
          hnAntiL_k1.SetBinContent((int[]){iS, iC, iE}, k1_antiL);
          if (k1_antiL > 0) {
            hnAntiL_k2k1.SetBinContent((int[]){iS, iC, iE}, k2_antiL / k1_antiL);
            hnAntiL_k2k1check.SetBinContent((int[]){iS, iC, iE}, (k2_antiL - k1_antiL) / std::pow(k1_antiL, 2.));
          }
          else {
            hnAntiL_k2k1.SetBinContent((int[]){iS, iC, iE}, 0.);
            hnAntiL_k2k1check.SetBinContent((int[]){iS, iC, iE}, 0.);
          }

          // anti-proton
          hnAntip_k2.SetBinContent((int[]){iS, iC, iE}, k2_antip);
          hnAntip_k1.SetBinContent((int[]){iS, iC, iE}, k1_antip);
          if (k1_antip > 0) {
            hnAntip_k2k1.SetBinContent((int[]){iS, iC, iE}, k2_antip / k1_antip);
            hnAntip_k2k1check.SetBinContent((int[]){iS, iC, iE}, (k2_antip - k1_antip) / std::pow(k1_antip, 2.));
          }
          else {
            hnAntip_k2k1.SetBinContent((int[]){iS, iC, iE}, 0.);
            hnAntip_k2k1check.SetBinContent((int[]){iS, iC, iE}, 0.);
          }

          // anti-deuteron
          hnAntid_k2.SetBinContent((int[]){iS, iC, iE}, k2_antid);
          hnAntid_k1.SetBinContent((int[]){iS, iC, iE}, k1_antid);
          if (k1_antid > 0) {
            hnAntid_k2k1.SetBinContent((int[]){iS, iC, iE}, k2_antid / k1_antid);
            hnAntid_k2k1check.SetBinContent((int[]){iS, iC, iE}, (k2_antid - k1_antid) / std::pow(k1_antid, 2.));
          }
          else {
            hnAntid_k2k1.SetBinContent((int[]){iS, iC, iE}, 0.);
            hnAntid_k2k1check.SetBinContent((int[]){iS, iC, iE}, 0.);
          }

          // net-lambda
          auto k2_netL = k2_L + k2_antiL - 2 * (k11_LantiL);
          if (k1_L + k1_antiL > 1.e-5) {
            hnNetL.SetBinContent((int[]){iS, iC, iE}, ( k2_netL ) / (k1_L + k1_antiL));
          }
          else {
            hnNetL.SetBinContent((int[]){iS, iC, iE}, 0.);
          }

          // net-lambda - deuteron
          if (k2_netL * k2_antid > 1.e-5) {
            hnAntidNetL.SetBinContent((int[]){iS, iC, iE}, ( k11_Lantid - k11_antiLantid ) / std::sqrt( k2_antid * k2_netL ));
          }
          else {
            hnAntidNetL.SetBinContent((int[]){iS, iC, iE}, 0.);
          }

          // Lambda-antideuteron
          if (k2_L * k2_antid > 0.) {
            hnLantid.SetBinContent((int[]){iS, iC, iE}, ( k11_Lantid ) / std::sqrt( k2_L * k2_antid ));
            // std::cout << "L-antid" << ( k11_Lantid - k1_L * k1_antid ) / std::sqrt( k2_L * k2_antid ) << std::endl;
          }
          else {
            hnLantid.SetBinContent((int[]){iS, iC, iE}, 0.);
          }

          // antiLambda-antideuteron
          if (k2_antiL * k2_antid > 0.) {
            hnAntiLantid.SetBinContent((int[]){iS, iC, iE}, ( k11_antiLantid ) / std::sqrt( k2_antiL * k2_antid ));
            // std::cout << "antiL-antid" << ( k11_antiLantid - k1_antiL * k1_antid ) / std::sqrt( k2_antiL * k2_antid ) << std::endl;
          }
          else {
            hnAntiLantid.SetBinContent((int[]){iS, iC, iE}, 0.);
          }

          // antiproton-antideuteron
          if (k2_antip * k2_antid > 1.e-9) {
            hnAntipAntid.SetBinContent((int[]){iS, iC, iE}, ( k11_antipAntid ) / std::sqrt( k2_antid * k2_antip ));
            // std::cout << "antip-antid" << ( k11_antipAntid - k1_antip * k1_antid ) / std::sqrt( k2_antid * k2_antip ) << std::endl;
          }
          else {
            hnAntipAntid.SetBinContent((int[]){iS, iC, iE}, 0.);
          }
        }
        else {
          hnL_k2k1.SetBinContent((int[]){iS, iC, iE}, 0.);
          hnL_k2k1check.SetBinContent((int[]){iS, iC, iE}, 0.);
          hnL_k1.SetBinContent((int[]){iS, iC, iE}, 0.);
          hnL_k2.SetBinContent((int[]){iS, iC, iE}, 0.);

          hnAntiL_k2k1.SetBinContent((int[]){iS, iC, iE}, 0.);
          hnAntiL_k2k1check.SetBinContent((int[]){iS, iC, iE}, 0.);
          hnAntiL_k1.SetBinContent((int[]){iS, iC, iE}, 0.);
          hnAntiL_k2.SetBinContent((int[]){iS, iC, iE}, 0.);

          hnAntip_k2k1.SetBinContent((int[]){iS, iC, iE}, 0.);
          hnAntip_k2k1check.SetBinContent((int[]){iS, iC, iE}, 0.);
          hnAntip_k1.SetBinContent((int[]){iS, iC, iE}, 0.);
          hnAntip_k2.SetBinContent((int[]){iS, iC, iE}, 0.);

          hnAntid_k2k1.SetBinContent((int[]){iS, iC, iE}, 0.);
          hnAntid_k2k1check.SetBinContent((int[]){iS, iC, iE}, 0.);
          hnAntid_k1.SetBinContent((int[]){iS, iC, iE}, 0.);
          hnAntid_k2.SetBinContent((int[]){iS, iC, iE}, 0.);

          hnNetL.SetBinContent((int[]){iS, iC, iE}, 0.);
          hnAntidNetL.SetBinContent((int[]){iS, iC, iE}, 0.);
          hnLantid.SetBinContent((int[]){iS, iC, iE}, 0.);
          hnAntiLantid.SetBinContent((int[]){iS, iC, iE}, 0.);
          hnAntipAntid.SetBinContent((int[]){iS, iC, iE}, 0.);
        }
      }
    }
    proj_nev->Clear();
  }



  subsample(&hnL_k2k1, hL_k2k1);
  subsample(&hnL_k2k1check, hL_k2k1check);
  subsample(&hnL_k1, hL_k1);
  subsample(&hnL_k2, hL_k2);
  subsample(&hnAntiL_k2k1, hAntiL_k2k1);
  subsample(&hnAntiL_k2k1check, hAntiL_k2k1check);
  subsample(&hnAntiL_k1, hAntiL_k1);
  subsample(&hnAntiL_k2, hAntiL_k2);
  subsample(&hnAntip_k2k1, hAntip_k2k1);
  subsample(&hnAntip_k2k1check, hAntip_k2k1check);
  subsample(&hnAntip_k1, hAntip_k1);
  subsample(&hnAntip_k2, hAntip_k2);
  subsample(&hnAntid_k2k1, hAntid_k2k1);
  subsample(&hnAntid_k2k1check, hAntid_k2k1check);
  subsample(&hnAntid_k1, hAntid_k1);
  subsample(&hnAntid_k2, hAntid_k2);
  subsample(&hnNetL, hNetL_k2k1);
  subsample(&hnAntidNetL, hAntidNetL);
  subsample(&hnLantid, hLantid);
  subsample(&hnAntiLantid, hAntiLantid);
  subsample(&hnAntipAntid, hAntipAntid);

  of->cd();
  hL_k2k1->Write();
  hL_k2k1check->Write();
  hL_k1->Write();
  hL_k2->Write();
  hAntiL_k2k1->Write();
  hAntiL_k2k1check->Write();
  hAntiL_k1->Write();
  hAntiL_k2->Write();
  hAntip_k2k1->Write();
  hAntip_k2k1check->Write();
  hAntip_k1->Write();
  hAntip_k2->Write();
  hAntid_k2k1->Write();
  hAntid_k2k1check->Write();
  hAntid_k1->Write();
  hAntid_k2->Write();
  hNetL_k2k1->Write();
  hAntidNetL->Write();
  hLantid->Write();
  hAntiLantid->Write();
  hAntipAntid->Write();

  of->Close();
}
