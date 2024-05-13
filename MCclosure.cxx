void MCclosure(){
  gStyle->SetOptStat(0);

  TFile *fgen = TFile::Open("_out_run2_LHC21k2abc_gen.root");
  TFile *frec = TFile::Open("_out_run2_LHC21k2abc_rec.root");
  TFile *fgenL = TFile::Open("_out_run2_LHC21d6abc_grid_gen.root");
  TFile *frecL = TFile::Open("_out_run2_LHC21d6abc_grid.root");
  TFile *fout = TFile::Open("_out_mcclosure.root", "recreate");

  TCanvas cRho("cRhoMCclosure", "cRhoMCclosure", 600, 600);
  TH1D* rhoGen = (TH1D*)fgen->Get("hAntidNetL");
  TH1D* rhoRec = (TH1D*)frec->Get("hAntidNetL");
  rhoRec->Add(rhoGen, -1);
  rhoRec->SetLineWidth(2);
  rhoRec->GetYaxis()->SetRangeUser(-0.1, 0.1);
  rhoRec->GetYaxis()->SetTitle("Rec - Gen (#rho_{#bar{d}#Delta#Lambda})");
  rhoRec->GetYaxis()->SetTitleOffset(1.6);
  rhoRec->GetXaxis()->SetRangeUser(0, 80);
  rhoRec->Fit("pol0");

  TCanvas cAntid("cAntidk2k1MCclosure", "cAntidk2k1MCclosure", 600, 600);
  TH1D* antidk2k1Gen = (TH1D*)fgen->Get("hAntid_k2k1");
  TH1D* antidk2k1Rec = (TH1D*)frec->Get("hAntid_k2k1");
  antidk2k1Rec->Divide(antidk2k1Gen);
  antidk2k1Rec->SetLineWidth(2);
  antidk2k1Rec->GetYaxis()->SetRangeUser(0.9, 1.1);
  antidk2k1Rec->GetYaxis()->SetTitle("Rec / Gen (#kappa_{2}/#kappa_{1}(#bar{d}))");
  antidk2k1Rec->GetYaxis()->SetTitleOffset(1.6);
  antidk2k1Rec->GetXaxis()->SetRangeUser(0, 80);
  antidk2k1Rec->Fit("pol0");

  TCanvas cNetL("cNetLk2k1MCclosure", "cNetLk2k1MCclosure", 600, 600);
  TH1D* netLk2k1Gen = (TH1D*)fgenL->Get("hNetL_k2k1");
  TH1D* netLk2k1Rec = (TH1D*)frecL->Get("hNetL_k2k1");
  netLk2k1Rec->Divide(netLk2k1Gen);
  netLk2k1Rec->SetLineWidth(2);
  netLk2k1Rec->GetYaxis()->SetRangeUser(0.9, 1.1);
  netLk2k1Rec->GetYaxis()->SetTitle("Rec / Gen (#kappa_{2}/#kappa_{1}(#Delta#Lambda))");
  netLk2k1Rec->GetYaxis()->SetTitleOffset(1.6);
  netLk2k1Rec->GetXaxis()->SetRangeUser(0, 80);
  netLk2k1Rec->Fit("pol0");

  cRho.cd();
  cRho.SetLeftMargin(0.16);
  cRho.SetRightMargin(0.05);
  cRho.SetTopMargin(0.05);
  rhoRec->Draw();

  cAntid.cd();
  cAntid.SetLeftMargin(0.16);
  cAntid.SetRightMargin(0.05);
  cAntid.SetTopMargin(0.05);
  antidk2k1Rec->Draw();

  cNetL.cd();
  cNetL.SetLeftMargin(0.16);
  cNetL.SetRightMargin(0.05);
  cNetL.SetTopMargin(0.05);
  netLk2k1Rec->Draw();

  fout->cd();
  rhoRec->Write();
  antidk2k1Rec->Write();
  netLk2k1Rec->Write();

  cRho.Print("cMcClosureRho.pdf");
  cAntid.Print("cMcClosureAntidk2k1.pdf");
  cNetL.Print("cMcClosureNetLk2k1.pdf");
}