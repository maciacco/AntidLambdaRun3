#include <TCanvas.h>

void SetOptCanv(TCanvas &cc){
  cc.SetLeftMargin(0.16);
  cc.SetRightMargin(0.03);
  cc.SetTopMargin(0.03);
}

constexpr const char* period = "18qr";

double chi2Matrices(TH1* const& hData, TGraphErrors* const& gSys, TGraphErrors* const& gModel, int size = 8);

void FinalPlotsWithModels(){
  gStyle->SetOptStat(0);
  TFile* f = TFile::Open(Form("_out_run2_LHC%s_grid_sys_hyperloop2.root", period));
  TFile* f16 = TFile::Open("out_full_analysis_antid_antip_highStats_ptL_1_4_16_2.root");
  TFile* f30 = TFile::Open("out_full_analysis_antid_antip_highStats_ptL_1_4_16_30_2.root");
  TFile* fSys = TFile::Open("outSys_newPtBins.root");
  TFile* fSysNetL = TFile::Open("outSys_newPtBins_netL.root");

  // antid-netL
  TH1D* hAntidNetL = (TH1D*)fSys->Get("CentralValue");
  auto hsys = (TH1D*)fSys->Get("sysTot");
  TGraphErrors gSys;
  for (int i{1}; i < 9; ++i) {
    gSys.AddPoint(hAntidNetL->GetBinCenter(i), hAntidNetL->GetBinContent(i));
    gSys.SetPointError(gSys.GetN() - 1, 1., hsys->GetBinContent(i));
  }
  gSys.SetLineWidth(3);
  gSys.SetFillStyle(0);
  gSys.SetLineColor(kRed);
  auto hCpy(*hAntidNetL);
  hCpy.Reset();
  hCpy.SetLineStyle(kDashed);
  hCpy.SetLineWidth(1);
  hCpy.SetLineColor(kBlack);
  hCpy.GetYaxis()->SetTitle("#it{#rho}_{#bar{d}#Delta#Lambda}");
  hCpy.GetYaxis()->SetTitleOffset(1.6);
  hCpy.GetYaxis()->SetNdivisions(505);
  hCpy.GetXaxis()->SetNdivisions(505);
  hCpy.GetYaxis()->SetRangeUser(-0.004, 0.012);
  hCpy.GetXaxis()->SetRangeUser(0., 80.);
  TGraphErrors* model_16_rho = (TGraphErrors*)f16->Get("Grrho");
  model_16_rho->SetName("model_16_rho");
  TGraphErrors* model_30_rho = (TGraphErrors*)f30->Get("Grrho");
  model_30_rho->SetName("model_30_rho");
  hAntidNetL->SetLineWidth(3);
  hAntidNetL->SetLineColor(kRed);
  hAntidNetL->SetMarkerColor(kRed);
  hAntidNetL->SetMarkerStyle(20);
  hAntidNetL->SetMarkerSize(1.3);
  TCanvas cNetLantid("rhoAntidNetL", "rhoAntidNetL", 600, 600);
  SetOptCanv(cNetLantid);
  cNetLantid.cd();
  // hAntidNetL->Fit("pol0");
  hCpy.Draw();
  model_30_rho->SetFillColor(kOrange);
  model_16_rho->SetFillColor(kAzure + 7);
  model_30_rho->SetLineColor(kOrange);
  model_16_rho->SetLineColor(kAzure + 7);
  model_30_rho->Draw("le3same");
  model_16_rho->Draw("le3same");
  gSys.Draw("e5same");
  hAntidNetL->Draw("pex0same");
  TLatex t;
  t.SetTextFont(44);
  t.SetTextSize(23);
  t.SetNDC();
  t.DrawLatex(0.2, 0.35 + 0.55, "ALICE Preliminary, Pb#minusPb #sqrt{#it{s}_{NN}} = 5.02 TeV");
  t.DrawLatex(0.2, 0.29 + 0.55, "|#it{#eta}| < 0.8");
  t.DrawLatex(0.2, 0.23 + 0.55, "0.6 < #it{p}_{T}(#bar{d}) < 1.8 GeV/#it{c}");
  t.DrawLatex(0.2, 0.17 + 0.55, "1 < #it{p}_{T}(#Lambda) < 4 GeV/#it{c}");
  TLegend l(0.2, 0.15, 0.8, 0.21);
  l.SetNColumns(2);
  l.SetColumnSeparation(0.2);
  l.SetTextFont(44);
  l.SetTextSize(23);
  t.DrawLatex(0.2, 0.25, "#splitline{Thermal-FIST CE SHM}{#it{T}_{ch}, #it{#gamma}_{s}, d#it{V}/d#it{y} from PRC 100 (2019) 054906}");
  l.AddEntry(model_16_rho, "#it{V}_{c} = 1.6 d#it{V}/d#it{y}", "fl");
  l.AddEntry(model_30_rho, "#it{V}_{c} = 3.0 d#it{V}/d#it{y}", "fl");
  l.Draw("same");
  auto chi2 = chi2Matrices(hAntidNetL, &gSys, model_16_rho);
  std::cout << chi2 / 8. << std::endl;
  std::cout << "nSigma = " << ROOT::Math::normal_quantile_c(TMath::Prob(chi2, 8)/2., 1.) << std::endl;
  chi2 = chi2Matrices(hAntidNetL, &gSys, model_30_rho);
  std::cout << chi2 / 8. << std::endl;
  std::cout << "nSigma = " << ROOT::Math::normal_quantile_c(TMath::Prob(chi2, 8)/2., 1.) << std::endl;
  cNetLantid.Print(Form("cNetLantidwModels_%s_modelsOnly.pdf", period));

  // netL
  TH1D* hNetL = (TH1D*)f->Get("hNetL_k2k1");
  auto hCpyNetL(*hNetL);
  auto hsysNetL = (TH1D*)fSysNetL->Get("sysTot");
  TGraphErrors gSysNetL;
  for (int i{1}; i < 9; ++i) {
    gSysNetL.AddPoint(hNetL->GetBinCenter(i), hNetL->GetBinContent(i));
    gSysNetL.SetPointError(gSysNetL.GetN() - 1, 1., hsysNetL->GetBinContent(i));
  }
  gSysNetL.SetLineWidth(3);
  gSysNetL.SetFillStyle(0);
  gSysNetL.SetLineColor(kRed);
  hCpyNetL.Reset();
  hCpyNetL.GetYaxis()->SetTitle("#it{#kappa}_{2}(#Lambda - #bar{#Lambda}) / #it{#kappa}_{1}(#Lambda + #bar{#Lambda})");
  hCpyNetL.GetYaxis()->SetTitleOffset(1.6);
  hCpyNetL.GetYaxis()->SetRangeUser(0.82, 1.04);
  hCpyNetL.GetYaxis()->SetNdivisions(505);
  hCpyNetL.GetXaxis()->SetNdivisions(505);
  hCpyNetL.GetXaxis()->SetRangeUser(0., 80.);
  TGraphErrors* model_16_c2c1 = (TGraphErrors*)f16->Get("Grc2byc1xi");
  model_16_c2c1->SetName("model_16_c2c1");
  TGraphErrors* model_30_c2c1 = (TGraphErrors*)f30->Get("Grc2byc1xi");
  model_30_c2c1->SetName("model_30_c2c1");
  hNetL->SetLineWidth(3);
  hNetL->SetLineColor(kRed);
  hNetL->SetMarkerColor(kRed);
  hNetL->SetMarkerStyle(20);
  hNetL->SetMarkerSize(1.3);
  TCanvas cNetL("netL", "netL", 600, 600);
  SetOptCanv(cNetL);
  cNetL.cd();
  //hNetL->Fit("pol0");
  hCpyNetL.Draw();
  model_30_c2c1->SetFillColor(kOrange);
  model_16_c2c1->SetFillColor(kAzure + 7);
  model_30_c2c1->SetLineColor(kOrange);
  model_16_c2c1->SetLineColor(kAzure + 7);
  model_30_c2c1->SetLineWidth(2);
  model_16_c2c1->SetLineWidth(2);
  model_30_c2c1->Draw("le3same");
  model_16_c2c1->Draw("le3same");
  hNetL->Draw("pex0same");
  gSysNetL.Draw("e5same");
  TLatex tNetL;
  tNetL.SetTextFont(44);
  tNetL.SetTextSize(23);
  tNetL.SetNDC();
  t.DrawLatex(0.2, 0.35 + 0.55, "ALICE Preliminary, Pb#minusPb #sqrt{#it{s}_{NN}} = 5.02 TeV");
  t.DrawLatex(0.2, 0.29 + 0.55, "|#it{#eta}| < 0.8, 1 < #it{p}_{T} < 4 GeV/#it{c}");
  TLegend lNetL(0.2, 0.15, 0.8, 0.21);
  lNetL.SetNColumns(2);
  lNetL.SetColumnSeparation(0.2);
  lNetL.SetTextFont(44);
  lNetL.SetTextSize(23);
  t.DrawLatex(0.2, 0.25, "#splitline{Thermal-FIST CE SHM}{#it{T}_{ch}, #it{#gamma}_{s}, d#it{V}/d#it{y} from PRC 100 (2019) 054906}");
  lNetL.AddEntry(model_16_c2c1, "#it{V}_{c} = 1.6 d#it{V}/d#it{y}", "fl");
  lNetL.AddEntry(model_30_c2c1, "#it{V}_{c} = 3.0 d#it{V}/d#it{y}", "fl");
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
    covariance[iDiag][iDiag] = std::pow(hData->GetBinError(iDiag + 1), 2.) + std::pow(gModel->GetErrorY(iDiag), 2.) /* + std::pow(gSys->GetErrorY(iDiag), 2.) */;
  }

  //fill correlated uncertainty matrix
  for (int iC{0}; iC < size; ++iC) {
    for (int iR{0}; iR < size; ++iR) {
      covariance[iR][iC] += (gSys->GetErrorY(iR) * gSys->GetErrorY(iC));
    }
  }

  // covariance.Print(std::cout);

  // combine errors
  auto inverse_covariance = covariance.Invert();

  return ROOT::Math::Similarity<double>(covariance, data_model_diff);
}
