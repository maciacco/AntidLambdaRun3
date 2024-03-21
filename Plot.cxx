void Plot(){
  TFile* f = TFile::Open("out_run2_evecuts.root");

  // correlations
  TH1D* hsame = (TH1D*)f->Get("hAntiLantid");
  TH1D* hopp = (TH1D*)f->Get("hLantid");
  hsame->SetLineWidth(2);
  hsame->SetLineColor(kBlue);
  hopp->SetLineWidth(2);
  hopp->SetLineColor(kRed);
  TCanvas c("rhoAntidLambda");
  c.cd();
  hsame->GetYaxis()->SetTitle("#rho_{d #Lambda}");
  hsame->Draw("pe");
  hopp->Draw("pesame");
  hsame->Fit("pol0");
  hopp->Fit("pol0");
  hsame->GetFunction("pol0")->SetLineColor(kBlue);
  hsame->GetFunction("pol0")->SetLineColor(kRed);
  hsame->GetYaxis()->SetRangeUser(-0.0015, 0.0015);
  c.Print("cSameOpp.pdf");

  // antid-antip
  TH1D* hAntipAntid = (TH1D*)f->Get("hAntipAntid");
  hAntipAntid->SetLineWidth(2);
  hAntipAntid->SetLineColor(kBlue);
  TCanvas cAntipAntid("rhoAntipAntid");
  cAntipAntid.cd();
  hAntipAntid->Fit("pol0");
  hAntipAntid->Draw("pe");
  cAntipAntid.Print("cAntipAntid.pdf");

  // antid-netL
  TH1D* hAntidNetL = (TH1D*)f->Get("hAntidNetL");
  hAntidNetL->SetLineWidth(2);
  hAntidNetL->SetLineColor(kBlue);
  TCanvas cNetLantid("rhoAntidNetL");
  cNetLantid.cd();
  hAntidNetL->Fit("pol0");
  hAntidNetL->Draw("pe");
  cNetLantid.Print("cNetLantid.pdf");

  // second-to-first-order cumulant ratio
  TH1D* hL = (TH1D*)f->Get("hL_k2k1");
  TCanvas cL("cLk2k1");
  hL->SetBinContent(10, 0.);
  hL->SetLineWidth(2);
  cL.cd();
  hL->GetYaxis()->SetTitle("#kappa_{2}/#kappa_{1}(#Lambda)");
  hL->Draw("pe");
  hL->GetYaxis()->SetRangeUser(0.99, 1.01);
  cL.Print("cLk2k1.pdf");

  TH1D* hAntiL = (TH1D*)f->Get("hAntiL_k2k1");
  TCanvas cAntiL("cAntiLk2k1");
  hAntiL->SetBinContent(10, 0.);
  hAntiL->SetLineWidth(2);
  cAntiL.cd();
  hAntiL->GetYaxis()->SetTitle("#kappa_{2}/#kappa_{1}(#bar{#Lambda})");
  hAntiL->Draw("pe");
  hAntiL->GetYaxis()->SetRangeUser(0.99, 1.01);
  cAntiL.Print("cAntiLk2k1.pdf");

  TH1D* hNetL = (TH1D*)f->Get("hNetL_k2k1");
  TCanvas cNetL("hNetLk2k1");
  hNetL->SetBinContent(10, 0.);
  hNetL->SetLineWidth(2);
  cNetL.cd();
  hNetL->GetYaxis()->SetTitle("#kappa_{2}/#kappa_{1}(#Delta#Lambda)");
  hNetL->Draw("pe");
  hNetL->GetYaxis()->SetRangeUser(0.99, 1.01);
  cNetL.Print("cNetLk2k1.pdf");

  TH1D* hAntip = (TH1D*)f->Get("hAntip_k2k1");
  TCanvas cAntip("cAntipk2k1");
  hAntip->SetBinContent(10, 0.);
  hAntip->SetLineWidth(2);
  cAntip.cd();
  hAntip->GetYaxis()->SetTitle("#kappa_{2}/#kappa_{1}(#bar{p})");
  hAntip->Draw("pe");
  hAntip->GetYaxis()->SetRangeUser(0.99, 1.01);
  cAntip.Print("cAntipk2k1.pdf");

  TH1D* hAntid = (TH1D*)f->Get("hAntid_k2k1");
  TCanvas cAntid("cAntidk2k1");
  hAntid->SetBinContent(10, 0.);
  hAntid->SetLineWidth(2);
  cAntid.cd();
  hAntid->GetYaxis()->SetTitle("#kappa_{2}/#kappa_{1}(#bar{d})");
  hAntid->Draw("pe");
  hAntid->GetYaxis()->SetRangeUser(0.99, 1.01);
  cAntid.Print("cAntidk2k1.pdf");
}
