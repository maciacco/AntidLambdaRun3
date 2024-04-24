#include <TCanvas.h>

void SetOptCanv(TCanvas &cc){
  cc.SetLeftMargin(0.16);
  cc.SetRightMargin(0.05);
  cc.SetTopMargin(0.05);
}

constexpr const char* period = "15o";

void FinalPlotsWithModels(){
  gStyle->SetOptStat(0);
  TFile* f = TFile::Open(Form("out_run2_LHC%s_20240420.root", period));
  TFile* f16 = TFile::Open("out_full_analysis_antid_antip_ptL_1_4_16.root");
  TFile* f30 = TFile::Open("out_full_analysis_antid_antip_ptL_1_4_30.root");

  // antid-netL
  TH1D* hAntidNetL = (TH1D*)f->Get("hAntidNetL");
  auto hCpy(*hAntidNetL);
  hCpy.Reset();
  hCpy.SetLineStyle(kDashed);
  hCpy.SetLineColor(kBlack);
  hCpy.GetYaxis()->SetTitleOffset(1.6);
  hCpy.GetYaxis()->SetRangeUser(-0.007, 0.013);
  hCpy.GetXaxis()->SetRangeUser(0., 80.);
  TGraphErrors* model_16_rho = (TGraphErrors*)f16->Get("Grrho");
  model_16_rho->SetName("model_16_rho");
  TGraphErrors* model_30_rho = (TGraphErrors*)f30->Get("Grrho");
  model_30_rho->SetName("model_30_rho");
  hAntidNetL->SetLineWidth(2);
  hAntidNetL->SetLineColor(kRed);
  TCanvas cNetLantid("rhoAntidNetL", "rhoAntidNetL", 700, 600);
  SetOptCanv(cNetLantid);
  cNetLantid.cd();
  //hAntidNetL->Fit("pol0");
  hCpy.Draw();
  model_30_rho->SetFillColor(kOrange - 2);
  model_16_rho->SetFillColor(kAzure + 7);
  model_30_rho->Draw("e3same");
  model_16_rho->Draw("e3same");
  hAntidNetL->Draw("pesame");
  TLatex t;
  t.SetTextFont(44);
  t.SetTextSize(23);
  t.SetNDC();
  t.DrawLatex(0.3, 0.33, "ALICE Pb#minusPb #sqrt{#it{s}_{NN}} = 5.02 TeV");
  t.DrawLatex(0.3, 0.26, "0.8 < #it{p}_{T}(#bar{d}) < 1.8 GeV/#it{c}");
  t.DrawLatex(0.3, 0.19, "1 < #it{p}_{T}(#Lambda) < 4 GeV/#it{c}");
  TLegend l(0.3, 0.8, 0.5, 0.9);
  l.SetTextFont(44);
  l.SetTextSize(23);
  l.AddEntry(model_16_rho, "Thermal-FIST, #it{V}_{c} = 1.6 d#it{V}/d#it{y}", "f");
  l.AddEntry(model_30_rho, "Thermal-FIST, #it{V}_{c} = 3.0 d#it{V}/d#it{y}", "f");
  l.Draw("same");
  cNetLantid.Print(Form("cNetLantidwModels_%s.pdf", period));

  // netL
  TH1D* hNetL = (TH1D*)f->Get("hNetL_k2k1");
  auto hCpyNetL(*hNetL);
  hCpyNetL.Reset();
  hCpyNetL.GetYaxis()->SetTitleOffset(1.6);
  hCpyNetL.GetYaxis()->SetRangeUser(0.79, 1.06);
  hCpyNetL.GetXaxis()->SetRangeUser(0., 80.);
  TGraphErrors* model_16_c2c1 = (TGraphErrors*)f16->Get("Grc2byc1xi");
  model_16_c2c1->SetName("model_16_c2c1");
  TGraphErrors* model_30_c2c1 = (TGraphErrors*)f30->Get("Grc2byc1xi");
  model_30_c2c1->SetName("model_30_c2c1");
  hNetL->SetLineWidth(2);
  hNetL->SetLineColor(kRed);
  TCanvas cNetL("netL", "netL", 700, 600);
  SetOptCanv(cNetL);
  cNetL.cd();
  //hNetL->Fit("pol0");
  hCpyNetL.Draw();
  model_30_c2c1->SetFillColor(kOrange - 2);
  model_16_c2c1->SetFillColor(kAzure + 7);
  model_30_c2c1->Draw("e3same");
  model_16_c2c1->Draw("e3same");
  hNetL->Draw("pesame");
  TLatex tNetL;
  tNetL.SetTextFont(44);
  tNetL.SetTextSize(23);
  tNetL.SetNDC();
  tNetL.DrawLatex(0.3, 0.33, "ALICE Pb#minusPb #sqrt{#it{s}_{NN}} = 5.02 TeV");
  //tNetL.DrawLatex(0.3, 0.26, "0.8 < #it{p}_{T}(#bar{d}) < 1.8 GeV/#it{c}");
  tNetL.DrawLatex(0.3, 0.19, "1 < #it{p}_{T}(#Lambda) < 4 GeV/#it{c}");
  TLegend lNetL(0.3, 0.8, 0.5, 0.9);
  lNetL.SetTextFont(44);
  lNetL.SetTextSize(23);
  lNetL.AddEntry(model_16_c2c1, "Thermal-FIST, #it{V}_{c} = 1.6 d#it{V}/d#it{y}", "f");
  lNetL.AddEntry(model_30_c2c1, "Thermal-FIST, #it{V}_{c} = 3.0 d#it{V}/d#it{y}", "f");
  lNetL.Draw("same");
  TLine lH(0., 1., 80., 1.);
  lH.SetLineStyle(kDashed);
  lH.SetLineColor(kBlack);
  lH.Draw("same");
  cNetL.Print(Form("cNetLwModels_%s.pdf", period));
}