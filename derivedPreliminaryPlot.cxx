void derivedPreliminaryPlot(){
  TFile fout("derived.root", "recreate");
  gStyle->SetOptStat(0);
  TCanvas c("c", "c", 600, 600);
  c.SetLeftMargin(0.16);
  c.SetRightMargin(0.03);
  c.SetTopMargin(0.03);
  TFile fmodel16("out_ptL_1_4_16.root");
  TFile fmodel30("out_ptL_1_4_30.root");
  auto gmodel16 = (TGraphErrors*)fmodel16.Get("Grc2byc1xi");
  auto gmodel30 = (TGraphErrors*)fmodel30.Get("Grc2byc1xi");

  auto grae = new TGraphAsymmErrors(8);
  grae->SetName("hK2NetLambdaSkellamRatioCentCorr_syst");
  grae->SetTitle("ratio to Skellam baseline");
  grae->SetPoint(0,5,0.9610578);
  grae->SetPointError(0,1,1,0.0296402,0.05931001);
  grae->SetPoint(1,15,0.9438433);
  grae->SetPointError(1,1,1,0.02740532,0.05578747);
  grae->SetPoint(2,25,0.9421404);
  grae->SetPointError(2,1,1,0.02728678,0.05510609);
  grae->SetPoint(3,35,0.9470533);
  grae->SetPointError(3,1,1,0.02750569,0.05580928);
  grae->SetPoint(4,45,0.9375322);
  grae->SetPointError(4,1,1,0.02737545,0.05544129);
  grae->SetPoint(5,55,0.9352931);
  grae->SetPointError(5,1,1,0.02791681,0.0546885);
  grae->SetPoint(6,65,0.940146);
  grae->SetPointError(6,1,1,0.03927298,0.05766634);
  grae->SetPoint(7,75,0.9457662);
  grae->SetPointError(7,1,1,0.0281375,0.06360881);

  auto gre = new TGraphErrors(8);
  gre->SetName("hK2NetLambdaSkellamRatioCentCorr_stat");
  gre->SetTitle("ratio to Skellam baseline");
  gre->SetMarkerStyle(20);
  gre->SetMarkerSize(1.5);
  gre->SetPoint(0,5,0.9610578);
  gre->SetPointError(0,0,0.004199425);
  gre->SetPoint(1,15,0.9438433);
  gre->SetPointError(1,0,0.004332259);
  gre->SetPoint(2,25,0.9421404);
  gre->SetPointError(2,0,0.00385588);
  gre->SetPoint(3,35,0.9470533);
  gre->SetPointError(3,0,0.003964213);
  gre->SetPoint(4,45,0.9375322);
  gre->SetPointError(4,0,0.004342889);
  gre->SetPoint(5,55,0.9352931);
  gre->SetPointError(5,0,0.007860323);
  gre->SetPoint(6,65,0.940146);
  gre->SetPointError(6,0,0.009775012);
  gre->SetPoint(7,75,0.9457662);
  gre->SetPointError(7,0,0.02483227);

  TH2F frame("frame", ";Centrality (%);#it{#kappa}_{2}(#Lambda - #bar{#Lambda}) / #it{#kappa}_{1}(#Lambda + #bar{#Lambda})", 1, 0, 80, 1, 0.86, 1.05);

  c.cd();
  frame.GetYaxis()->SetNdivisions(505);
  frame.GetXaxis()->SetNdivisions(505);
  frame.Draw();
  frame.GetYaxis()->SetTitleOffset(1.6);
  gre->SetLineWidth(3);
  gre->SetLineColor(kViolet - 3);
  gre->SetMarkerColor(kViolet - 3);
  gre->SetMarkerStyle(53);
  gre->SetMarkerSize(1.3);
  grae->SetLineWidth(0);
  grae->SetFillStyle(3154);
  //grae->SetLineColor(kViolet - 3);
  grae->SetFillColor(kViolet - 9);
  gmodel30->SetLineColor(kOrange - 2);
  gmodel16->SetLineColor(kAzure + 7);
  gmodel30->SetFillColor(kOrange - 2);
  gmodel16->SetFillColor(kAzure + 7);
  gmodel16->Draw("e3same");
  gmodel30->Draw("e3same");
  grae->Draw("e5same");
  gre->Draw("pesame");

  // new analysis
  TFile* f = TFile::Open("_out_run2_LHC18qr_lambda_small_acc.root");
  TFile* fSysNetL = TFile::Open("outSys_newPtBins_netL_smallAcc.root");
  TH1D* hNetL = (TH1D*)f->Get("hNetL_k2k1");
  auto hsysNetL = (TH1D*)fSysNetL->Get("sysTot");
  TGraphErrors gSysNetL;
  for (int i{1}; i < 9; ++i) {
    gSysNetL.AddPoint(hNetL->GetBinCenter(i), hNetL->GetBinContent(i));
    gSysNetL.SetPointError(gSysNetL.GetN() - 1, 1., hsysNetL->GetBinContent(i));
  }
  gSysNetL.SetLineWidth(3);
  gSysNetL.SetFillStyle(0);
  gSysNetL.SetLineColor(kRed);
  hNetL->SetLineWidth(3);
  hNetL->SetLineColor(kRed);
  hNetL->SetMarkerColor(kRed);
  hNetL->SetMarkerStyle(20);
  hNetL->SetMarkerSize(1.3);
  hNetL->Draw("pex0same");
  gSysNetL.Draw("e5same");
  TLine hl(0, 1, 80, 1);
  hl.SetLineStyle(kDashed);
  hl.Draw("same");

  TLegend data(0.25, 0.76, 0.5, 0.82);
  data.SetNColumns(2);
  data.SetColumnSeparation(0.9);
  data.SetTextFont(44);
  data.SetTextSize(23);
  data.AddEntry(gre, "Run 2 2015", "pe");
  data.AddEntry(hNetL, "Run 2 2018", "pe");
  data.Draw("same");

  TLatex t;
  t.SetTextFont(44);
  t.SetTextSize(23);
  t.SetNDC();
  t.DrawLatex(0.2, 0.35 + 0.55, "ALICE Preliminary, Pb#minusPb #sqrt{#it{s}_{NN}} = 5.02 TeV");
  t.DrawLatex(0.25, 0.29 + 0.55, "|#it{#eta}| < 0.5, 1 < #it{p}_{T} < 4 GeV/#it{c}");
  TLegend lNetL(0.2, 0.15, 0.8, 0.21);
  lNetL.SetNColumns(2);
  lNetL.SetColumnSeparation(0.2);
  lNetL.SetTextFont(44);
  lNetL.SetTextSize(23);
  t.DrawLatex(0.2, 0.25, "#splitline{Thermal-FIST CE SHM}{#it{T}_{ch}, #it{#gamma}_{s}, d#it{V}/d#it{y} from PRC 100 (2019) 054906}");
  lNetL.AddEntry(gmodel16, "#it{V}_{c} = 1.6 d#it{V}/d#it{y}", "fl");
  lNetL.AddEntry(gmodel30, "#it{V}_{c} = 3.0 d#it{V}/d#it{y}", "fl");
  lNetL.Draw("same");

  c.Print("derived.pdf");
  fout.cd();
  c.Write();
}
