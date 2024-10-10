#include "Constants.h"

int colors[] = {TColor::GetColor("#ff3300"), TColor::GetColor("#ec6e0a"), TColor::GetColor("#daaa14"), TColor::GetColor("#c7e51e"), TColor::GetColor("#85dd69"), TColor::GetColor("#42d6b4"), TColor::GetColor("#00ceff"), TColor::GetColor("#009adf"), TColor::GetColor("#0067c0"), TColor::GetColor("#595959"), TColor::GetColor("#0033a1"), kBlack};
//const char* var[]{"cosPA", "dcaLambdaDaugh", "dcaLambdatopv", "massLambda", "tpccls", "pidtof", "pidtpc", "dcaTrtackToPv", "zvtx"};
//const char* var[]{"_sys1", "_sys2", "_sys3", "_sys4", "_sys5", "_sys6", "_sys7", "_sys8", "_sys9", "_sys10", "_sys11", "_sys12", "_sys13", "_sys14", "_sys15", "_sys16", "_sys17", "_sys18", "_sys19", "_sys20", "_sys21", "_sys22", "_sys23", "_sys24", "_sys25", "_sys26", "_sys27", "_sys28", "_sys29", "_sys30"};
const char* var[] = {"", "_sys01", "_sys02", "_sys03", "_sys04", "_sys05", "_sys06", "_sys07", "_sys08", "_sys09", "_sys10", "_sys11", "_sys12", "_sys13", "_sys14", "_sys15", "_sys16", "_sys17", "_sys18", "_sys19", "_sys20", "_sys21", "_sys22", "_sys23", "_sys24", "_sys25", "_sys26", "_sys27", "_sys28", "_sys29", "_sys30", "_sys31", "_sys32", "_sys33", "_sys34", "_sys35", "_sys36", "_sys37", "_sys38", "_sys39", "_sys40", "_sys41", "_sys42", "_sys43", "_sys44", "_sys45", "_sys46", "_sys47", "_sys48", "_sys49", "_sys50"};


bool combine = false;
bool doBarlow = false;

void SystematicsMultitrialNetL(){
  gStyle->SetOptStat(0);
  TFile* files[51];
  TH1D* rho_tmp[51];
  TH1D* k2k1_tmp[51];
  TH1D* rho[51];
  TH1D* rhoRatio[51];
  TH1D* k2k1[51];
  TH1D* sys[9];
  TGraphErrors gout;
  TH1D hout("hout", ";Centrality (%);Sys", 8, 0, 80);
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
  leg.SetTextSize(15);

  gStyle->SetPalette(kBlueYellow);
  TLegend legVar(0.2, 0.16, 0.7, 0.36);
  legVar.SetNColumns(2);
  legVar.SetTextFont(44);
  legVar.SetTextSize(17);
  legVar.SetNColumns(3);

  TFile out("outSys_multitrial_netL.root", "recreate");
  std::string outFile{"out_run2_Multitrial.root"};//"out_run2_LHC18qr_grid_sys_hyperloop2.root"};

  std::string outFile0{"out_run2_Multitrial.root"};//"out_run2_LHC18qr_grid_sys_hyperloop2.root"};
  TFile *files0 = TFile::Open(Form("_%s", outFile0.c_str()));
  auto rho0 = (TH1D*)files0->Get("hNetL_k2k1");
  auto k2k10 = (TH1D*)files0->Get("hNetL_k2k1");
  rho0->SetName(Form("hNetL_k2k1%s", ""));
  k2k10->SetName(Form("hNetL_k2k1%s", ""));
  TH1D *rhoCpy = new TH1D("rhoCpy", ";Centrality (%);#rho_{#bar{d}#Delta#Lambda}", 8, 0, 80);
  TH1D *k2k1Cpy = new TH1D("k2k1Cpy", ";Centrality (%);#kappa_{2}/#kappa_{1}(#Delta#Lambda)}", 8, 0, 80);
  for (int i{1}; i < 9; ++i) {
    rhoCpy->SetBinContent(i, rho0->GetBinContent(i));
    rhoCpy->SetBinError(i, rho0->GetBinError(i));
    k2k1Cpy->SetBinContent(i, k2k10->GetBinContent(i));
    k2k1Cpy->SetBinError(i, k2k10->GetBinError(i));
  }
  rho0->GetYaxis()->SetRangeUser(0.8, 1.2);
  rho0->GetXaxis()->SetRangeUser(0, 80);
  rho0->SetLineWidth(2);
  rho0->SetMarkerStyle(71);
  rho0->SetMarkerSize(1.5);
  // if (i == 0) rho0->SetLineColor(kBlack);
  out.cd();
  rho0->Write();

  double corr_frac = 0.;

  for (int i{0}; i < 50; ++i) {
    auto cutSet = cutSets[i+1];
    std::cout << cutSet << std::endl;
    files[i] = TFile::Open(Form("%s_%s", cutSet, outFile.c_str()));
    rho[i] = (TH1D*)files[i]->Get("hNetL_k2k1");
    k2k1[i] = (TH1D*)files[i]->Get("hNetL_k2k1");
    rho[i]->SetName(Form("hNetL_k2k1%s", cutSet));
    k2k1[i]->SetName(Form("hNetL_k2k1%s", cutSet));

    out.cd();
    //rho[i]->Fit("pol0");
    cvar.cd();
    rho[i]->GetYaxis()->SetRangeUser(0.8, 1.2);
    rho[i]->GetXaxis()->SetRangeUser(0, 80);
    rho[i]->SetLineWidth(2);
    rho[i]->SetMarkerStyle(71);
    rho[i]->SetMarkerSize(1.5);
    // if (i == 0) rho[i]->SetLineColor(kBlack);
    rho[i]->Draw(i == 1 ? "PLC PMC" : "same PLC PMC");
    rho[i]->Write();
    rhoRatio[i] = new TH1D(*rho[i]);
    rhoRatio[i]->SetName(Form("rhoRatio_%s", cutSet));
    rhoRatio[i]->Divide(rhoCpy);
    double mean = 0;
    double sq = 0;
    for (int iP{1}; iP < 9; ++iP){
      rhoRatio[i]->SetBinContent(iP, rho[i]->GetBinContent(iP) - rhoCpy->GetBinContent(iP));
      rhoRatio[i]->SetBinError(iP, std::sqrt(std::abs(std::pow(rho[i]->GetBinError(iP), 2.) /* - std::pow(rhoCpy->GetBinError(iP), 2.) */)));
      mean += rho[i]->GetBinContent(iP) - rhoCpy->GetBinContent(iP);
      sq += std::pow(rho[i]->GetBinContent(iP) - rhoCpy->GetBinContent(iP), 2.);
    }
    mean /= 8.;
    sq /= 8.;
    sq = std::sqrt((sq - std::pow(mean, 2.))/8.);
    std::cout << mean << " +/- " << sq << std::endl;
    rhoRatio[i]->GetYaxis()->SetRangeUser(0.8, 1.2);
    rhoRatio[i]->Fit("pol0", "Q");
    for (int iP{1}; iP < 9; ++iP) {
      rhoRatio[i]->SetBinContent(iP, (rhoRatio[i]->GetBinContent(iP)) - rhoRatio[i]->GetFunction("pol0")->GetParameter(0));
      rhoRatio[i]->SetBinError(iP, std::sqrt(std::abs(std::pow(rho[i]->GetBinError(iP), 2.) - std::pow(rhoRatio[i]->GetFunction("pol0")->GetParError(0), 2.))));
    }
    // rhoRatio[i]->Fit("pol0", "Q");
    std::cout << rhoRatio[i]->GetFunction("pol0")->GetParameter(0) << std::endl;
    double test = /* mean / sq; // */rhoRatio[i]->GetFunction("pol0")->GetParameter(0) / rhoRatio[i]->GetFunction("pol0")->GetParError(0);
    if (std::abs(test) > 3.) {
      std::cout << i << std::endl;
      corr_frac += 1.;
    }
    rhoRatio[i]->Write();
    //legVar.AddEntry(rho[i], cutSets[i]);
  }

  corr_frac /= 50.;
  std::cout << "correlated fraction = " << corr_frac << std::endl;

  for (int i{0}; i < 8; ++i) {
    sys[i] = new TH1D(Form("hsys_%d", i), ";Var;Entries", 400, -0.05, 0.05);
    for (int ii{0}; ii < 50; ++ii) {
      sys[i]->Fill(rhoRatio[ii]->GetBinContent(i + 1));
    }
  }

  for (int i{8}; i < 9; ++i) {
    sys[i] = new TH1D(Form("hsys_%d", i), ";Var;Entries", 800, -0.05, 0.05);
    for (int ii{0}; ii < 50; ++ii) {
      sys[i]->Fill(rhoRatio[ii]->GetFunction("pol0")->GetParameter(0));
    }
    sys[i]->Write();
  }

  out.cd();
  for (int i{0}; i < 8; ++i) {
    auto mm = sys[i]->GetMean();
    double dist = 100000;
    int j_var = 0;
    for (int ii{0}; ii < 50; ++ ii) {
      auto mm_tmp = rho[i]->GetBinContent(i + 1);
      auto tmp_dist = std::abs(mm - mm_tmp);
      if (tmp_dist < dist) {
        dist = tmp_dist;
        j_var = ii;
      }
    }
    gout.AddPoint(rho0->GetBinCenter(i+1), rho[j_var]->GetBinContent(i+1));
    gout.SetPointError(gout.GetN() - 1, 1., rho[j_var]->GetBinError(i+1));
    //gout.AddPoint(rho0->GetBinCenter(i+1), rho0->GetBinContent(i+1));
    //gout.SetPointError(gout.GetN() - 1, 1., sys[i]->GetStdDev());
    hout.SetBinContent(i + 1, sys[i]->GetStdDev());
    sys[i]->Write();
  }
  hout.GetXaxis()->SetRangeUser(-10000, 90000);
  // hout.Smooth(1, "R");
  hout.Write();
  gout.Write("gout");
  //leg.AddEntry(sys[9], "Total");
  //leg.Draw("same");
  //csys.Print("csys.pdf");
}
