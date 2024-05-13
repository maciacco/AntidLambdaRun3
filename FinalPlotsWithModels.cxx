#include <TCanvas.h>

void SetOptCanv(TCanvas &cc){
  cc.SetLeftMargin(0.16);
  cc.SetRightMargin(0.05);
  cc.SetTopMargin(0.05);
}

constexpr const char* period = "18qr";

double chi2Matrices(TH1* const& hData, TGraphErrors* const& gSys, TGraphErrors* const& gModel, int size = 8);

void FinalPlotsWithModels(){
  gStyle->SetOptStat(0);
  TFile* f = TFile::Open(Form("_out_run2_LHC%s_grid_sys_hyperloop2.root", period));
  TFile* f16 = TFile::Open("out_full_analysis_antid_antip_highStats_ptL_1_4_16.root");
  TFile* f30 = TFile::Open("out_full_analysis_antid_antip_highStats_ptL_1_4_30.root");
  TFile* fSys = TFile::Open("outSys_newPtBins.root");

  // antid-netL
  TH1D* hAntidNetL = (TH1D*)fSys->Get("CentralValue");
  auto hsys = (TH1D*)fSys->Get("sysTot");
  TGraphErrors gSys;
  for (int i{1}; i < 9; ++i) {
    gSys.AddPoint(hAntidNetL->GetBinCenter(i), hAntidNetL->GetBinContent(i));
    gSys.SetPointError(gSys.GetN() - 1, 2., hsys->GetBinContent(i));
  }
  gSys.SetLineWidth(2);
  gSys.SetFillStyle(0);
  gSys.SetLineColor(kRed);
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
  gSys.Draw("e5same");
  hAntidNetL->Draw("pesame");
  TLatex t;
  t.SetTextFont(44);
  t.SetTextSize(23);
  t.SetNDC();
  t.DrawLatex(0.2, 0.33, "ALICE Preliminary, Pb#minusPb #sqrt{#it{s}_{NN}} = 5.02 TeV");
  t.DrawLatex(0.2, 0.26, "0.6 < #it{p}_{T}(#bar{d}) < 1.8 GeV/#it{c}");
  t.DrawLatex(0.2, 0.19, "1 < #it{p}_{T}(#Lambda) < 4 GeV/#it{c}");
  TLegend l(0.3, 0.8, 0.5, 0.9);
  l.SetTextFont(44);
  l.SetTextSize(23);
  l.AddEntry(model_16_rho, "Thermal-FIST, #it{V}_{c} = 1.6 d#it{V}/d#it{y}", "f");
  l.AddEntry(model_30_rho, "Thermal-FIST, #it{V}_{c} = 3.0 d#it{V}/d#it{y}", "f");
  l.Draw("same");
  auto chi2 = chi2Matrices(hAntidNetL, &gSys, model_16_rho);
  std::cout << chi2 / 8. << std::endl;
  chi2 = chi2Matrices(hAntidNetL, &gSys, model_30_rho);
  std::cout << chi2 / 8. << std::endl;
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

double chi2Matrices(TH1* const& hData, TGraphErrors* const& gSys, TGraphErrors* const& gModel, int size = 8){
  ROOT::Math::SVector<double, 8> data_model_diff;
  ROOT::Math::SMatrix<double, 8, 8> covariance;

  // fill data, model, and uncorrelated uncertainty matrices
  for (int iDiag{0}; iDiag < size; ++iDiag) {
    data_model_diff[iDiag] = hData->GetBinContent(iDiag + 1) - gModel->GetPointY(iDiag);
    covariance[iDiag][iDiag] = std::pow(hData->GetBinError(iDiag + 1), 2.)+ std::pow(gModel->GetErrorY(iDiag), 2.) /* + std::pow(gSys->GetErrorY(iDiag), 2.) */;
  }

  //fill correlated uncertainty matrix
  for (int iC{0}; iC < size; ++iC) {
    for (int iR{0}; iR < size; ++iR) {
      covariance[iR][iC] += gSys->GetErrorY(iR) * gSys->GetErrorY(iC);
    }
  }

  // covariance.Print(std::cout);

  // combine errors
  auto inverse_covariance = covariance.Invert();

  return ROOT::Math::Similarity<double>(covariance, data_model_diff);
}