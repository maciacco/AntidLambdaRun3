#include <TCanvas.h>

void SetOptCanv(TCanvas &cc){
  cc.SetLeftMargin(0.16);
  cc.SetRightMargin(0.05);
  cc.SetTopMargin(0.05);
}

constexpr const char* period = "15o";

void Plot(){
  TFile* f = TFile::Open(Form("out_run2_LHC%s_20240420_noEff.root", period));
  gStyle->SetOptStat(0);

  // correlations
  TH1D* hsame = (TH1D*)f->Get("hAntiLantid");
  TH1D* hopp = (TH1D*)f->Get("hLantid");
  hsame->SetLineWidth(2);
  hsame->SetLineColor(kBlue);
  hopp->SetLineWidth(2);
  hopp->SetLineColor(kRed);
  TCanvas c("rhoAntidLambda", "rhoAntidLambda", 700, 600);
  SetOptCanv(c);
  c.cd();
  hsame->GetYaxis()->SetTitle("#rho_{d #Lambda}");
  hsame->Draw("pe");
  hopp->Draw("pesame");
  hsame->Fit("pol0");
  hopp->Fit("pol0");
  hsame->GetFunction("pol0")->SetLineColor(kBlue);
  hsame->GetFunction("pol0")->SetLineColor(kRed);
  hsame->GetYaxis()->SetRangeUser(-0.0015, 0.0015);
  hsame->GetXaxis()->SetRangeUser(0., 80.);
  c.Print(Form("cSameOpp_%s.pdf", period));

  // antid-antip
  TH1D* hAntipAntid = (TH1D*)f->Get("hAntipAntid");
  hAntipAntid->SetLineWidth(2);
  hAntipAntid->SetLineColor(kBlue);
  TCanvas cAntipAntid("rhoAntipAntid");
  cAntipAntid.cd();
  hAntipAntid->Fit("pol0");
  hAntipAntid->Draw("pe");
  hAntipAntid->GetXaxis()->SetRangeUser(0., 80.);
  cAntipAntid.Print(Form("cAntipAntid_%s.pdf", period));

  // antid-netL
  TH1D* hAntidNetL = (TH1D*)f->Get("hAntidNetL");
  hAntidNetL->SetLineWidth(2);
  hAntidNetL->SetLineColor(kBlue);
  TCanvas cNetLantid("rhoAntidNetL", "rhoAntidNetL", 700, 600);
  SetOptCanv(cNetLantid);
  cNetLantid.cd();
  hAntidNetL->GetYaxis()->SetTitleOffset(1.6);
  hAntidNetL->GetYaxis()->SetRangeUser(-0.001, 0.002);
  hAntidNetL->GetXaxis()->SetRangeUser(0., 80.);
  //hAntidNetL->Fit("pol0");
  hAntidNetL->Draw("pe");
  cNetLantid.Print(Form("cNetLantid_%s.pdf", period));

  // second-to-first-order cumulant ratio
  TH1D* hL = (TH1D*)f->Get("hL_k2k1");
  TCanvas cL("cLk2k1");
  hL->SetBinContent(10, 0.);
  hL->SetLineWidth(2);
  cL.cd();
  hL->GetYaxis()->SetTitle("#kappa_{2}/#kappa_{1}(#Lambda)");
  hL->Draw("pe");
  hL->GetYaxis()->SetRangeUser(0.99, 1.01);
  hL->GetXaxis()->SetRangeUser(0., 80.);
  cL.Print(Form("cLk2k1_%s.pdf", period));

  TH1D* hAntiL = (TH1D*)f->Get("hAntiL_k2k1");
  TCanvas cAntiL("cAntiLk2k1");
  hAntiL->SetBinContent(10, 0.);
  hAntiL->SetLineWidth(2);
  cAntiL.cd();
  hAntiL->GetYaxis()->SetTitle("#kappa_{2}/#kappa_{1}(#bar{#Lambda})");
  hAntiL->Draw("pe");
  hAntiL->GetYaxis()->SetRangeUser(0.99, 1.01);
  hAntiL->GetXaxis()->SetRangeUser(0., 80.);
  cAntiL.Print(Form("cAntiLk2k1_%s.pdf", period));

  TH1D* hNetL = (TH1D*)f->Get("hNetL_k2k1");
  TCanvas cNetL("hNetLk2k1", "hNetLk2k1", 700, 600);
  SetOptCanv(cNetL);
  hNetL->SetBinContent(10, 0.);
  hNetL->SetLineWidth(2);
  hNetL->GetYaxis()->SetTitleOffset(1.6);
  cNetL.cd();
  hNetL->GetYaxis()->SetTitle("#kappa_{2}/#kappa_{1}(#Delta#Lambda)");
  hNetL->Draw("pe");
  hNetL->GetYaxis()->SetRangeUser(0.99, 1.001);
  hNetL->GetXaxis()->SetRangeUser(0., 80.);
  cNetL.Print(Form("cNetLk2k1_%s.pdf", period));

  TH1D* hAntip = (TH1D*)f->Get("hAntip_k2k1");
  TCanvas cAntip("cAntipk2k1");
  hAntip->SetBinContent(10, 0.);
  hAntip->SetLineWidth(2);
  cAntip.cd();
  hAntip->GetYaxis()->SetTitle("#kappa_{2}/#kappa_{1}(#bar{p})");
  hAntip->Draw("pe");
  hAntip->GetYaxis()->SetRangeUser(0.99, 1.01);
  hAntip->GetXaxis()->SetRangeUser(0., 80.);
  cAntip.Print(Form("cAntipk2k1_%s.pdf", period));

  TH1D* hAntid = (TH1D*)f->Get("hAntid_k2k1");
  TCanvas cAntid("cAntidk2k1", "cAntidk2k1", 700, 600);
  SetOptCanv(cAntid);
  hAntid->GetYaxis()->SetTitleOffset(1.6);
  hAntid->SetBinContent(10, 0.);
  hAntid->SetLineWidth(2);
  cAntid.cd();
  hAntid->GetYaxis()->SetTitle("#kappa_{2}/#kappa_{1}(#bar{d})");
  hAntid->GetXaxis()->SetRangeUser(0., 80.);
  hAntid->Draw("pe");
  hAntid->GetYaxis()->SetRangeUser(0.99, 1.01);
  cAntid.Print(Form("cAntidk2k1_%s.pdf", period));
}
