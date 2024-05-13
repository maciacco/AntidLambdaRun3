#include <TFile.h>
#include <THnSparse.h>
#include <TH1D.h>
#include <TColor.h>
#include <TStyle.h>
#include <TH3F.h>
#include <TH2F.h>
#include <Riostream.h>

#include "Constants.h"

const double scale_eff = 0.95;

void subsample(THnSparse *nPart, TH1D *hOut){
  TH1D* nPartVsCent = nullptr;
  auto nSubsamples = nPart->GetAxis(0)->GetNbins();
  auto nCent = nPart->GetAxis(1)->GetNbins();
  for (int iC{1}; iC < nCent + 1; ++iC) {
    for (int iS{0}; iS < nSubsamples; ++iS) {
      int idx2[]{iS + 1, iC, 1};
      auto binc = nPart->GetBinContent(idx2);
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

void RawCorrelation(const int cut = 0, const bool useRescaledMB = false){
  const char* cutSet = cutSets[cut];
  auto f = TFile::Open(fileIO::inFileName.data());
  auto of = TFile::Open(Form("%s_%s", cutSet, fileIO::outFileName.data()), "recreate");
  auto file_eff = TFile::Open(Form("efficiency_%s_%s.root", period, cutSet));
  // auto file_eff_MB = TFile::Open("efficiencyMB2018.root");
  auto nEv = dynamic_cast<THnSparse*>(f->Get(Form("antid-lambda-ebye%s/nEv", cutSet)));
  auto nAntid = dynamic_cast<THnSparse*>(f->Get(Form("antid-lambda-ebye%s/n%sAntid",cutSet, genRec)));
  auto nAntip = dynamic_cast<THnSparse*>(f->Get(Form("antid-lambda-ebye%s/n%sAntip",cutSet, genRec)));
  auto nAntiL = dynamic_cast<THnSparse*>(f->Get(Form("antid-lambda-ebye%s/n%sAntiL",cutSet, genRec)));
  auto nL = dynamic_cast<THnSparse*>(f->Get(Form("antid-lambda-ebye%s/n%sL",cutSet, genRec)));
  auto nSqAntid = dynamic_cast<THnSparse*>(f->Get(Form("antid-lambda-ebye%s/n%sSqAntid",cutSet, genRec)));
  auto nSqAntip = dynamic_cast<THnSparse*>(f->Get(Form("antid-lambda-ebye%s/n%sSqAntip",cutSet, genRec)));
  auto nSqAntiL = dynamic_cast<THnSparse*>(f->Get(Form("antid-lambda-ebye%s/n%sSqAntiL",cutSet, genRec)));
  auto nSqL = dynamic_cast<THnSparse*>(f->Get(Form("antid-lambda-ebye%s/n%sSqL",cutSet, genRec)));
  auto nLAntiL = dynamic_cast<THnSparse*>(f->Get(Form("antid-lambda-ebye%s/n%sLantiL",cutSet, genRec)));
  auto nLAntid = dynamic_cast<THnSparse*>(f->Get(Form("antid-lambda-ebye%s/n%sLantid",cutSet, genRec)));
  auto nAntiLAntid = dynamic_cast<THnSparse*>(f->Get(Form("antid-lambda-ebye%s/n%sAntiLantid",cutSet, genRec)));
  auto nAntipAntid = dynamic_cast<THnSparse*>(f->Get(Form("antid-lambda-ebye%s/n%sAntipAntid",cutSet, genRec)));

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
  TH1D *hNetL_k2 = new TH1D("hNetL_k2", ";Centrality (%);#kappa_{2}", bins::bins[1], bins::xmin[1], bins::xmax[1]);
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
  THnSparseD hnNetL_k2("hNetL_k2", ";Subsample;Centrality (%);#Delta#eta;", 3, bins::bins, bins::xmin, bins::xmax);
  THnSparseD hnLantid("hLantid", ";Subsample;Centrality (%);#Delta#eta;", 3, bins::bins, bins::xmin, bins::xmax);
  THnSparseD hnAntidNetL("hAntidNetL", ";Subsample;Centrality (%);#Delta#eta;", 3, bins::bins, bins::xmin, bins::xmax);
  THnSparseD hnAntiLantid("hAntiLantid", ";Subsample;Centrality (%);#Delta#eta;", 3, bins::bins, bins::xmin, bins::xmax);
  THnSparseD hnAntipAntid("hAntipAntid", ";Subsample;Centrality (%);#Delta#eta;", 3, bins::bins, bins::xmin, bins::xmax);

  // loop over bins
  TH1D *proj_nev = nullptr, *proj_antiL = nullptr, *proj_L = nullptr, *proj_antid = nullptr, *proj_antip = nullptr;
  TH2D *proj_sqAntiL = nullptr, *proj_sqL = nullptr, *proj_sqAntip = nullptr, *proj_sqAntid = nullptr;
  TH2D *proj_LantiL = nullptr, *proj_Lantid = nullptr, *proj_antiLantid = nullptr, *proj_antipAntid = nullptr;
  for (int iS{1}; iS < bins::bins[0] + 1; ++iS) {
    int iSsmallMin = nEv->GetAxis(0)->FindBin(hnAntid_k1.GetAxis(0)->GetBinLowEdge(iS) + constants::epsilon);
    int iSsmallMax = nEv->GetAxis(0)->FindBin(hnAntid_k1.GetAxis(0)->GetBinUpEdge(iS) - constants::epsilon);
    std::cout << "(iSsmallMin, iSsmallMax) = (" << iSsmallMin << ", " << iSsmallMax << ")" << "...\n";
    setAxisRanges(nEv, {iSsmallMin}, {iSsmallMax});
    proj_nev = dynamic_cast<TH1D*>(nEv->Projection(1));
    std::cout << " --- S U B S A M P L E    N . " << iS << " --- \n";
    for (int iC{1}; iC < /* bins::bins[1] - 1 */ bins::bins[1] + 1; ++iC) { // bins::bins[1] - 1 to be used for centrality differential measurement
      auto effd = (TH1D*)file_eff->Get(Form("effD_%d", iC));
      auto effL = (TH1D*)file_eff->Get(Form("effL_%d", iC));
      auto effAntiL = (TH1D*)file_eff->Get(Form("effAntiL_%d", iC));
      auto effp = (TH1D*)file_eff->Get(Form("effP_%d", iC));
      // if (useRescaledMB) {
      //   effd = (TH1D*)file_eff_MB->Get("effD");
      //   effL = (TH1D*)file_eff_MB->Get("effL");
      //   effAntiL = (TH1D*)file_eff_MB->Get("effAntiL");
      //   effp = (TH1D*)file_eff_MB->Get("effP");
      //   auto ratioToMB_effd = (TH1D*)file_eff->Get(Form("ratioToMB_D_%d", iC));
      //   auto ratioToMB_effL = (TH1D*)file_eff->Get(Form("ratioToMB_L_%d", iC));
      //   auto ratioToMB_effAntiL = (TH1D*)file_eff->Get(Form("ratioToMB_AntiL_%d", iC));
      //   auto ratioToMB_effp = (TH1D*)file_eff->Get(Form("ratioToMB_P_%d", iC));
      //   effd->Multiply(ratioToMB_effd);
      //   effL->Multiply(ratioToMB_effL);
      //   effAntiL->Multiply(ratioToMB_effAntiL);
      //   effp->Multiply(ratioToMB_effp);
      // }
      effd->SetName("effd");
      effL->SetName("effL");
      effAntiL->SetName("effAntiL");
      effp->SetName("effp");
      int iCsmallMin = nEv->GetAxis(1)->FindBin(hAntid_k1->GetXaxis()->GetBinLowEdge(iC) + constants::epsilon);
      int iCsmallMax = nEv->GetAxis(1)->FindBin(hAntid_k1->GetXaxis()->GetBinUpEdge(iC) - constants::epsilon);
      double nev_tot = proj_nev->Integral(iCsmallMin, iCsmallMax); // number of events in large centrality bin
      std::cout << "nev_tot = " << nev_tot << "\n";
      std::cout << "Processing centrality " << iC << ", (iCsmallMin, iCsmallMax) = (" << iCsmallMin << ", " << iCsmallMax << ")" << "...\n";
      for (int iE{1}; iE < bins::bins[2] + 1; ++iE) {

        int idx[]{iS, iC, iE};
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
          double nev = proj_nev->GetBinContent(iC_small); // number of events in small centrality bin
          std::cout << "nev_small = " << nev << "\n";

          // select axis ranges
          setAxisRanges(nAntid, {iSsmallMin, iC_small, 1}, {iSsmallMax, iC_small, 1});
          setAxisRanges(nAntip, {iSsmallMin, iC_small, 1}, {iSsmallMax, iC_small, 1});
          setAxisRanges(nAntiL, {iSsmallMin, iC_small, 1}, {iSsmallMax, iC_small, 1});
          setAxisRanges(nL, {iSsmallMin, iC_small, 1}, {iSsmallMax, iC_small, 1});
          setAxisRanges(nSqAntid, {iSsmallMin, iC_small, 1}, {iSsmallMax, iC_small, 1});
          setAxisRanges(nSqAntip, {iSsmallMin, iC_small, 1}, {iSsmallMax, iC_small, 1});
          setAxisRanges(nSqAntiL, {iSsmallMin, iC_small, 1}, {iSsmallMax, iC_small, 1});
          setAxisRanges(nSqL, {iSsmallMin, iC_small, 1}, {iSsmallMax, iC_small, 1});
          setAxisRanges(nLAntiL, {iSsmallMin, iC_small, 1}, {iSsmallMax, iC_small, 1});
          setAxisRanges(nLAntid, {iSsmallMin, iC_small, 1}, {iSsmallMax, iC_small, 1});
          setAxisRanges(nAntiLAntid, {iSsmallMin, iC_small, 1}, {iSsmallMax, iC_small, 1});
          setAxisRanges(nAntipAntid, {iSsmallMin, iC_small, 1}, {iSsmallMax, iC_small, 1});

          double sum_antip = 0.;
          double sum_antip_2 = 0.;
          double sumSq_antip = 0.;
          double sum_antid = 0.;
          double sum_antid_2 = 0.;
          double sumSq_antid = 0.;
          double sum_antiL = 0.;
          double sum_antiL_2 = 0.;
          double sumSq_antiL = 0.;
          double sum_L = 0.;
          double sum_L_2 = 0.;
          double sumSq_L = 0.;
          double sum_LantiL = 0.;
          double sum_Lantid = 0.;
          double sum_antiLantid = 0.;
          double sum_antipAntid = 0.;

          // project on pT axes
          proj_antid = dynamic_cast<TH1D*>(nAntid->Projection(3));
          proj_sqAntid = dynamic_cast<TH2D*>(nSqAntid->Projection(4, 3));
          proj_antip = dynamic_cast<TH1D*>(nAntip->Projection(3));
          proj_sqAntip = dynamic_cast<TH2D*>(nSqAntip->Projection(4, 3));
          proj_antiL = dynamic_cast<TH1D*>(nAntiL->Projection(3));
          proj_sqAntiL = dynamic_cast<TH2D*>(nSqAntiL->Projection(4, 3));
          proj_L = dynamic_cast<TH1D*>(nL->Projection(3));
          proj_sqL = dynamic_cast<TH2D*>(nSqL->Projection(4, 3));
          proj_LantiL = dynamic_cast<TH2D*>(nLAntiL->Projection(4, 3));
          proj_Lantid = dynamic_cast<TH2D*>(nLAntid->Projection(4, 3));
          proj_antiLantid = dynamic_cast<TH2D*>(nAntiLAntid->Projection(4, 3));
          proj_antipAntid = dynamic_cast<TH2D*>(nAntipAntid->Projection(4, 3));

          for (int iPL{1}; iPL < proj_antiL->GetNbinsX() + 1; ++iPL) { // loop over pT bins (lambda)
            double n_L = proj_L->GetBinContent(iPL);
            double n_antiL = proj_antiL->GetBinContent(iPL);
            double effL_pt = kApplyEffCorrection ? effL->GetBinContent(iPL) * scale_eff : 1.;
            double effAntiL_pt = kApplyEffCorrection ? effAntiL->GetBinContent(iPL) * scale_eff : 1.; // TODO: apply charge-conjugate efficiencies
            sum_L += (n_L / effL_pt);
            std::cout << n_L / effL_pt << std::endl;
            sum_antiL += (n_antiL / effAntiL_pt);
            sum_L_2 += (n_L / std::pow(effL_pt, 2.));
            sum_antiL_2 += (n_antiL / std::pow(effAntiL_pt, 2.));
            for (int iPL2{1}; iPL2 < proj_antiL->GetNbinsX() + 1; ++iPL2) { // loop over pT bins (lambda, 2)
              double nsq_L = proj_sqL->GetBinContent(iPL, iPL2);
              double nsq_antiL = proj_sqAntiL->GetBinContent(iPL, iPL2);
              double effL_pt2 = kApplyEffCorrection ? effL->GetBinContent(iPL2) * scale_eff : 1.; // TODO: apply charge-conjugate efficiencies
              double effAntiL_pt2 = kApplyEffCorrection ? effAntiL->GetBinContent(iPL2) * scale_eff : 1.; // TODO: apply charge-conjugate efficiencies
              double n_LantiL = proj_LantiL->GetBinContent(iPL, iPL2);
              sumSq_L += (nsq_L / effL_pt / effL_pt2);
              sumSq_antiL += (nsq_antiL / effAntiL_pt / effAntiL_pt2);
              sum_LantiL += (n_LantiL / effL_pt / effAntiL_pt2);
            }
            for (int iPD{1}; iPD < proj_antid->GetNbinsX() + 1; ++iPD) { // loop over pT bins (deuteron)
              double n_Lantid = proj_Lantid->GetBinContent(iPL, iPD);
              double effd_pt = kApplyEffCorrection ? effd->GetBinContent(iPD) * scale_eff : 1.;
              double n_antiLantid = proj_antiLantid->GetBinContent(iPL, iPD);
              // std::cout << "n_antiLantid = " << n_antiLantid << std::endl;
              sum_Lantid += (n_Lantid / effL_pt / effd_pt);
              sum_antiLantid += (n_antiLantid / effAntiL_pt / effd_pt);
            }
          }
          for (int iPD{1}; iPD < proj_antid->GetNbinsX() + 1; ++iPD) { // loop over pT bins (deuteron)
            double n_antid = proj_antid->GetBinContent(iPD);
            double effd_pt = kApplyEffCorrection ? effd->GetBinContent(iPD) * scale_eff : 1.;
            sum_antid += (n_antid / effd_pt);
            sum_antid_2 += (n_antid / std::pow(effd_pt, 2.));
            for (int iPD2{1}; iPD2 < proj_antid->GetNbinsX() + 1; ++iPD2) { // loop over pT bins (deuteron, 2)
              double nsq_antid = proj_sqAntid->GetBinContent(iPD, iPD2);
              double effd_pt2 = kApplyEffCorrection ? effd->GetBinContent(iPD2) * scale_eff : 1.;
              sumSq_antid += (nsq_antid / effd_pt / effd_pt2);
            }
          }
          for (int iPP{1}; iPP < proj_antip->GetNbinsX() + 1; ++iPP) { // loop over pT bins (proton)
            double n_antip = proj_antip->GetBinContent(iPP);
            double effp_pt = kApplyEffCorrection ? effp->GetBinContent(iPP) * scale_eff : 1.;
            sum_antip += (n_antip / effp_pt);
            sum_antip_2 += (n_antip / std::pow(effp_pt, 2.));
            // std::cout << "eff_pt = " << effp_pt << std::endl;
            for (int iPP2{1}; iPP2 < proj_antip->GetNbinsX() + 1; ++iPP2) { // loop over pT bins (proton, 2)
              double nsq_antip = proj_sqAntip->GetBinContent(iPP, iPP2);
              double effp_pt2 = kApplyEffCorrection ? effp->GetBinContent(iPP2) * scale_eff : 1.;
              sumSq_antip += (nsq_antip / effp_pt2 / effp_pt);
            }
            for (int iPD2{1}; iPD2 < proj_antid->GetNbinsX() + 1; ++iPD2) { // loop over pT bins (deuteron, 2)
              double n_antipAntid = proj_antipAntid->GetBinContent(iPP, iPD2);
              double effd_pt2 = kApplyEffCorrection ? effd->GetBinContent(iPD2) * scale_eff : 1.;
              sum_antipAntid += (n_antipAntid / effd_pt2 / effp_pt);
            }
          }

          if (nev > 1.e-5) {
            k1_L += (sum_L);
            k2_L += (sumSq_L - std::pow(sum_L, 2.) / nev + sum_L - sum_L_2);
            k1_antiL += (sum_antiL);
            k2_antiL += (sumSq_antiL - std::pow(sum_antiL, 2.) / nev + sum_antiL - sum_antiL_2);
            k1_antid += (sum_antid);
            k2_antid += (sumSq_antid - std::pow(sum_antid, 2.) / nev + sum_antid - sum_antid_2);
            k1_antip += (sum_antip);
            k2_antip += (sumSq_antip - std::pow(sum_antip, 2.) / nev + sum_antip - sum_antip_2);
            k11_LantiL += (sum_LantiL - sum_L * sum_antiL / nev);
            k11_Lantid += (sum_Lantid - sum_L * sum_antid / nev);
            k11_antiLantid += (sum_antiLantid - sum_antiL * sum_antid / nev);
            k11_antipAntid += (sum_antipAntid - sum_antip * sum_antid / nev);

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

          std::cout << k1_L << ", k11_Lantid = " << k11_Lantid << ", k11_antiLantid = " << k11_antiLantid << "\n";

          // std::cout << "k11_antiLantid = " << k11_antiLantid << ", k2_antiL * k2_antid = " << k2_antiL * k2_antid << std::endl;
          // std::cout << "k11_Lantid = " << k11_Lantid << ", k2_L * k2_antid = " << k2_L * k2_antid << std::endl;
          // std::cout << "k11_antipAntid = " << k11_antipAntid << std::endl;

          // lambda
          hnL_k2.SetBinContent(idx, k2_L);
          hnL_k1.SetBinContent(idx, k1_L);
          if (k1_L > 0) {
            hnL_k2k1.SetBinContent(idx, k2_L / k1_L);
            hnL_k2k1check.SetBinContent(idx, (k2_L - k1_L) / std::pow(k1_L, 2.));
          }
          else {
            hnL_k2k1.SetBinContent(idx, 0.);
            hnL_k2k1check.SetBinContent(idx, 0.);
          }

          // anti-lambda
          hnAntiL_k2.SetBinContent(idx, k2_antiL);
          hnAntiL_k1.SetBinContent(idx, k1_antiL);
          if (k1_antiL > 0) {
            hnAntiL_k2k1.SetBinContent(idx, k2_antiL / k1_antiL);
            hnAntiL_k2k1check.SetBinContent(idx, (k2_antiL - k1_antiL) / std::pow(k1_antiL, 2.));
          }
          else {
            hnAntiL_k2k1.SetBinContent(idx, 0.);
            hnAntiL_k2k1check.SetBinContent(idx, 0.);
          }

          // anti-proton
          hnAntip_k2.SetBinContent(idx, k2_antip);
          hnAntip_k1.SetBinContent(idx, k1_antip);
          if (k1_antip > 0) {
            hnAntip_k2k1.SetBinContent(idx, k2_antip / k1_antip);
            hnAntip_k2k1check.SetBinContent(idx, (k2_antip - k1_antip) / std::pow(k1_antip, 2.));
          }
          else {
            hnAntip_k2k1.SetBinContent(idx, 0.);
            hnAntip_k2k1check.SetBinContent(idx, 0.);
          }

          // anti-deuteron
          hnAntid_k2.SetBinContent(idx, k2_antid);
          hnAntid_k1.SetBinContent(idx, k1_antid);
          if (k1_antid > 0) {
            hnAntid_k2k1.SetBinContent(idx, k2_antid / k1_antid);
            hnAntid_k2k1check.SetBinContent(idx, (k2_antid - k1_antid) / std::pow(k1_antid, 2.));
          }
          else {
            hnAntid_k2k1.SetBinContent(idx, 0.);
            hnAntid_k2k1check.SetBinContent(idx, 0.);
          }

          // net-lambda
          auto k2_netL = k2_L + k2_antiL - 2 * (k11_LantiL);
          if (k1_L + k1_antiL > 1.e-5) {
            hnNetL.SetBinContent(idx, ( k2_netL ) / (k1_L + k1_antiL));
          }
          else {
            hnNetL.SetBinContent(idx, 0.);
          }

          if (k1_L + k1_antiL > 1.e-5) {
            hnNetL_k2.SetBinContent(idx, k2_netL);
          }
          else {
            hnNetL_k2.SetBinContent(idx, 0.);
          }

          // net-lambda - deuteron
          if (k2_netL * k2_antid > 1.e-5) {
            hnAntidNetL.SetBinContent(idx, ( k11_Lantid - k11_antiLantid ) / std::sqrt( k2_antid * k2_netL ));
          }
          else {
            hnAntidNetL.SetBinContent(idx, 0.);
          }

          // Lambda-antideuteron
          if (k2_L * k2_antid > 0.) {
            hnLantid.SetBinContent(idx, ( k11_Lantid ) / std::sqrt( k2_L * k2_antid ));
            // std::cout << "L-antid" << ( k11_Lantid - k1_L * k1_antid ) / std::sqrt( k2_L * k2_antid ) << std::endl;
          }
          else {
            hnLantid.SetBinContent(idx, 0.);
          }

          // antiLambda-antideuteron
          if (k2_antiL * k2_antid > 0.) {
            hnAntiLantid.SetBinContent(idx, ( k11_antiLantid ) / std::sqrt( k2_antiL * k2_antid ));
            // std::cout << "antiL-antid" << ( k11_antiLantid - k1_antiL * k1_antid ) / std::sqrt( k2_antiL * k2_antid ) << std::endl;
          }
          else {
            hnAntiLantid.SetBinContent(idx, 0.);
          }

          // antiproton-antideuteron
          if (k2_antip * k2_antid > 1.e-9) {
            hnAntipAntid.SetBinContent(idx, ( k11_antipAntid ) / std::sqrt( k2_antid * k2_antip ));
            // std::cout << "antip-antid" << ( k11_antipAntid - k1_antip * k1_antid ) / std::sqrt( k2_antid * k2_antip ) << std::endl;
          }
          else {
            hnAntipAntid.SetBinContent(idx, 0.);
          }
        }
        else {
          hnL_k2k1.SetBinContent(idx, 0.);
          hnL_k2k1check.SetBinContent(idx, 0.);
          hnL_k1.SetBinContent(idx, 0.);
          hnL_k2.SetBinContent(idx, 0.);

          hnAntiL_k2k1.SetBinContent(idx, 0.);
          hnAntiL_k2k1check.SetBinContent(idx, 0.);
          hnAntiL_k1.SetBinContent(idx, 0.);
          hnAntiL_k2.SetBinContent(idx, 0.);

          hnAntip_k2k1.SetBinContent(idx, 0.);
          hnAntip_k2k1check.SetBinContent(idx, 0.);
          hnAntip_k1.SetBinContent(idx, 0.);
          hnAntip_k2.SetBinContent(idx, 0.);

          hnAntid_k2k1.SetBinContent(idx, 0.);
          hnAntid_k2k1check.SetBinContent(idx, 0.);
          hnAntid_k1.SetBinContent(idx, 0.);
          hnAntid_k2.SetBinContent(idx, 0.);

          hnNetL.SetBinContent(idx, 0.);
          hnNetL_k2.SetBinContent(idx, 0.);
          hnAntidNetL.SetBinContent(idx, 0.);
          hnLantid.SetBinContent(idx, 0.);
          hnAntiLantid.SetBinContent(idx, 0.);
          hnAntipAntid.SetBinContent(idx, 0.);
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
  subsample(&hnNetL_k2, hNetL_k2);
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
  hNetL_k2->Write();
  hAntidNetL->Write();
  hLantid->Write();
  hAntiLantid->Write();
  hAntipAntid->Write();

  of->Close();
}
