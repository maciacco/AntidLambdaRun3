#include "Constants.h"

int colors[] = {TColor::GetColor("#ff3300"), TColor::GetColor("#ec6e0a"), TColor::GetColor("#daaa14"), TColor::GetColor("#c7e51e"), TColor::GetColor("#85dd69"), TColor::GetColor("#42d6b4"), TColor::GetColor("#00ceff"), TColor::GetColor("#009adf"), TColor::GetColor("#0067c0"), TColor::GetColor("#595959"), TColor::GetColor("#0033a1"), kBlack};
const char* var[]{"cos#theta_{p}", "DCA_{tracks}", "DCA_{V^{0}PV}", "#it{M}", "#it{n}_{TPCCls}", "n#sigma_{TOF}", "n#sigma_{TPC}", "DCA", "V_{z}"};

bool combine = false;
bool doBarlow = false;

void SystematicsNetL(){
  gStyle->SetOptStat(0);
  TFile* files[2][17];
  TH1D* rho_tmp[2][17];
  TH1D* k2k1_tmp[2][17];
  TH1D* rho[17];
  TH1D* k2k1[17];
  TH1D* sys[10];
  TCanvas csys("csys", "csys", 600, 600);
  TCanvas cvar("cvar", "cvar", 600, 600);
  csys.SetTopMargin(0.03);
  csys.SetRightMargin(0.03);
  csys.SetLeftMargin(0.16);
  csys.SetBottomMargin(0.14);
  cvar.SetTopMargin(0.03);
  cvar.SetRightMargin(0.03);
  cvar.SetLeftMargin(0.16);
  cvar.SetBottomMargin(0.14);
  TLegend leg(0.2, 0.7, 0.7, 0.92);
  leg.SetNColumns(2);
  leg.SetTextFont(44);
  leg.SetTextSize(19);

  gStyle->SetPalette(kBlueYellow);
  TLegend legVar(0.2, 0.76, 0.8, 0.92);
  legVar.SetNColumns(2);
  legVar.SetTextFont(44);
  legVar.SetTextSize(17);
  legVar.SetNColumns(3);
  TFile out("outSys_newPtBins_netL.root", "recreate");
  std::string outFile[]{"out_run2_LHC15o_grid_sys_hyperloop2.root", /* "out_run2_LHC18qr_grid_sys_hyperloop2.root" *//*  "out_run2_LHC18qr_lambda_small_acc.root" */"out_run2_LHC18qr_grid_sys_hyperloop2.root"};
  if (combine) {
    for (int iP{0}; iP < 2; ++iP) {
      for (int i{0}; i < 17; ++i) {
        auto cutSet = cutSets[i];
        std::cout << cutSet << std::endl;
        files[iP][i] = TFile::Open(Form("%s_%s", cutSet, outFile[iP].c_str()));
        // rho[i] = (TH1D*)files[iP][i]->Get("hAntidNetL");
        // k2k1[i] = (TH1D*)files[iP][i]->Get("hNetL_k2k1");
        // rho[i]->SetName(Form("hAntidNetL%s", cutSet));
        // k2k1[i]->SetName(Form("hNetL_k2k1%s", cutSet));
        rho_tmp[iP][i] = (TH1D*)files[iP][i]->Get("hAntidNetL");
        k2k1_tmp[iP][i] = (TH1D*)files[iP][i]->Get("hNetL_k2k1");
        rho_tmp[iP][i]->SetName(Form("hAntidNetL%s", cutSet));
        k2k1_tmp[iP][i]->SetName(Form("hNetL_k2k1%s", cutSet));
      }
    }
    for (int i{0}; i < 17; ++i) {
      rho[i] = new TH1D(*rho_tmp[0][i]);
      k2k1[i] = new TH1D(*k2k1_tmp[0][i]);
      for (int iB{1}; iB <= 8; ++ iB) {
        double rho1 = rho_tmp[0][i]->GetBinContent(iB);
        double rho2 = rho_tmp[1][i]->GetBinContent(iB);
        double rhoe1 = rho_tmp[0][i]->GetBinError(iB);
        double rhoe2 = rho_tmp[1][i]->GetBinError(iB);
        rho[i]->SetBinContent(iB, (rho1 / std::pow(rhoe1, 2.) + rho2 / std::pow(rhoe2, 2.)) / (1. / std::pow(rhoe1, 2.) + 1. / std::pow(rhoe2, 2.)));
        rho[i]->SetBinError(iB, std::sqrt(1. / (1. / std::pow(rhoe1, 2.) + 1. / std::pow(rhoe2, 2.))));
        double k2k11 = k2k1_tmp[0][i]->GetBinContent(iB);
        double k2k12 = k2k1_tmp[1][i]->GetBinContent(iB);
        double k2k1e1 = k2k1_tmp[0][i]->GetBinError(iB);
        double k2k1e2 = k2k1_tmp[1][i]->GetBinError(iB);
        k2k1[i]->SetBinContent(iB, (k2k11 / std::pow(k2k1e1, 2.) + k2k12 / std::pow(k2k1e2, 2.)) / (1. / std::pow(k2k1e1, 2.) + 1. / std::pow(k2k1e2, 2.)));
        k2k1[i]->SetBinError(iB, std::sqrt(1. / (1. / std::pow(k2k1e1, 2.) + 1. / std::pow(k2k1e2, 2.))));
      }
      out.cd();
      rho[i]->Write();
    }
  } else {
    int iP = 1;
    for (int i{1}; i < 17; ++i) {
      if (i == 11) continue;
      auto cutSet = cutSets[i];
      std::cout << cutSet << std::endl;
      files[iP][i] = TFile::Open(Form("%s_%s", cutSet, outFile[iP].c_str()));
      rho[i] = (TH1D*)files[iP][i]->Get("hAntidNetL");
      k2k1[i] = (TH1D*)files[iP][i]->Get("hNetL_k2k1");
      rho[i]->SetName(Form("hAntidNetL%s", cutSet));
      k2k1[i]->SetName(Form("hNetL_k2k1%s", cutSet));

      out.cd();
      //rho[i]->Fit("pol0");

      if (i > 8 && i < 15) continue;
      cvar.cd();
      k2k1[i]->GetYaxis()->SetRangeUser(0.9, 1.01);
      k2k1[i]->GetXaxis()->SetRangeUser(0, 80);
      k2k1[i]->SetLineWidth(2);
      k2k1[i]->SetMarkerStyle(71);
      k2k1[i]->SetMarkerSize(1.5);
      // if (i == 0) k2k1[i]->SetLineColor(kBlack);
      k2k1[i]->Draw(i == 1 ? "PLC PMC" : "same PLC PMC");
      k2k1[i]->Write();
      legVar.AddEntry(k2k1[i], cutSetsName[i]);
    }

    for (int i{0}; i < 1; ++i) {
      auto cutSet = cutSets[i];
      std::cout << cutSet << std::endl;
      files[iP][i] = TFile::Open(Form("%s_%s", cutSet, outFile[iP].c_str()));
      rho[i] = (TH1D*)files[iP][i]->Get("hAntidNetL");
      k2k1[i] = (TH1D*)files[iP][i]->Get("hNetL_k2k1");
      rho[i]->SetName(Form("hAntidNetL%s", cutSet));
      k2k1[i]->SetName(Form("hNetL_k2k1%s", cutSet));

      out.cd();
      //rho[i]->Fit("pol0");
      cvar.cd();
      k2k1[i]->GetYaxis()->SetRangeUser(0.8, 1.2);
      k2k1[i]->GetXaxis()->SetRangeUser(0, 80);
      k2k1[i]->SetLineWidth(3);
      // k2k1[i]->SetMarkerStyle(20);
      // k2k1[i]->SetMarkerSize(1.5);
      k2k1[i]->SetLineColor(kRed);
      k2k1[i]->SetMarkerColor(kRed);
      k2k1[i]->Draw("same");
      k2k1[i]->Write();
      legVar.AddEntry(k2k1[i], "Nominal");
    }
  }
  legVar.Draw();
  cvar.Print("cvar.pdf");
  for (int i{0}; i < 9; ++i) {
    if (i > 3 && i < 8) continue;
    sys[i] = new TH1D(Form("sys_%d", i), ";Centrality (%);Sys. uncertainty", 8, 0, 80);
    sys[i]->SetName(Form("sys_%d", i));
    if (i < 7 && !(/* i == 1 || i == 2 || */ i == 5)) {
      for (int iB{1}; iB < 9; ++iB) {
        auto barlowLow = std::abs(k2k1[2 * i + 1]->GetBinContent(iB) - k2k1[0]->GetBinContent(iB)) / std::sqrt(std::abs(std::pow(k2k1[2 * i + 1]->GetBinError(iB), 2.) - std::pow(k2k1[0]->GetBinError(iB), 2.)));
        auto barlowUp = std::abs(k2k1[2 * (i + 1)]->GetBinContent(iB) - k2k1[0]->GetBinContent(iB)) / std::sqrt(std::abs(std::pow(k2k1[2 * (i + 1)]->GetBinError(iB), 2.) - std::pow(k2k1[0]->GetBinError(iB), 2.)));
        if (barlowUp < 2 && barlowLow < 2 && doBarlow) continue;
        //std::cout << k2k1[2 * i + 1]->GetBinContent(iB) - k2k1[2 * (i + 1)]->GetBinContent(iB) << std::endl;
        sys[i]->SetBinContent(iB, std::abs(k2k1[2 * i + 1]->GetBinContent(iB) - k2k1[2 * (i + 1)]->GetBinContent(iB)) * 0.5);
        sys[i]->SetBinError(iB, 0.);
      }
    }
    else if (i > 6){
      std::cout << i << std::endl;
      for (int iB{1}; iB < 9; ++iB) {
        auto barlow = std::abs(k2k1[8 + i]->GetBinContent(iB) - k2k1[0]->GetBinContent(iB)) / std::sqrt(std::abs(std::pow(k2k1[8 + i]->GetBinError(iB), 2.) - std::pow(k2k1[0]->GetBinError(iB), 2.)));
        if (barlow < 2 && doBarlow) continue;
        sys[i]->SetBinContent(iB, std::abs(k2k1[8 + i]->GetBinContent(iB) - k2k1[0]->GetBinContent(iB)));
        sys[i]->SetBinError(iB, 0.);
      }
    }
    else if (/* i == 1 || i == 2 || */ i == 5){
      for (int iB{1}; iB < 9; ++iB) {
        auto barlow = std::abs(k2k1[2 * (i + 1)]->GetBinContent(iB) - k2k1[0]->GetBinContent(iB)) / std::sqrt(std::abs(std::pow(k2k1[2 * (i + 1)]->GetBinError(iB), 2.) - std::pow(k2k1[0]->GetBinError(iB), 2.)));
        if (barlow < 2 && doBarlow) continue;
        sys[i]->SetBinContent(iB, std::abs(k2k1[2 * (i + 1)]->GetBinContent(iB) - k2k1[0]->GetBinContent(iB)));
        sys[i]->SetBinError(iB, 0.);
      }
    }
    out.cd();
    sys[i]->Smooth();
    sys[i]->Write();
    sys[i]->SetLineWidth(3);
    sys[i]->SetLineColor(colors[i]);
    csys.cd();
    sys[i]->GetYaxis()->SetRangeUser(0, 0.02);
    sys[i]->GetYaxis()->SetTitleOffset(1.7);
    leg.AddEntry(sys[i], var[i]);
    sys[i]->Draw(i == 0 ? "histo" : "histosame");
  }

  sys[9] = new TH1D("sysTot", ";Centrality (%);#k2k1", 8, 0, 80);
  for (int iB{1}; iB < 9; ++iB) {
    double sum = 0.;
    for (int i{0}; i < 9; ++i) {
    if (i > 3 && i < 8) continue;
    if (i == 5) continue;
      sum += std::pow(sys[i]->GetBinContent(iB), 2.);
    }
    sys[9]->SetBinContent(iB, std::sqrt(sum));
  }
  out.cd();
  k2k1[0]->Write("CentralValue");
  sys[9]->Smooth();
  sys[9]->Write();
  csys.cd();
  csys.Write();
  sys[9]->SetLineWidth(3);
  sys[9]->SetLineColor(kBlack);
  sys[9]->Draw("same");
  leg.AddEntry(sys[9], "Total");
  leg.Draw("same");
  csys.Print("csys.pdf");
}