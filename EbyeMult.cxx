void SetHistStyle(TH1 *h)
{
  h->GetYaxis()->SetTitleFont(44);
  h->GetYaxis()->SetTitleSize(26);
  h->GetYaxis()->SetTitleOffset(1.5);
  h->GetXaxis()->SetTitleFont(44);
  h->GetXaxis()->SetTitleSize(26);
  h->GetXaxis()->SetTitleOffset(.9);
  h->SetMarkerStyle(20);
  h->SetMarkerSize(1.2);
  h->SetMarkerColor(kBlack);
  h->SetLineColor(kBlack);
  h->SetLineWidth(2);
}

void shift_binning(const TH1* in, TH1* out, double shift = -0.5){
  int n_bins = in->GetNbinsX();
  double left = in->GetXaxis()->GetBinLowEdge(1);
  double right = in->GetXaxis()->GetBinUpEdge(n_bins);
  out->SetName(Form("%s_shift", in->GetName()));
  out->SetBins(n_bins, left + shift, right + shift);
  out->GetXaxis()->SetTitle(in->GetXaxis()->GetTitle());
  out->GetYaxis()->SetTitle(in->GetYaxis()->GetTitle());
  for (int iB = 0; iB < n_bins; ++iB){
    out->SetBinContent(iB + 1, in->GetBinContent(iB + 1));
    out->SetBinError(iB + 1, in->GetBinError(iB + 1));
  }
}

constexpr const char* PartLabel[]{"#bar{d}", "#bar{#Lambda}", "#Lambda"};
constexpr const char* PartText[]{"Antid", "AntiL", "L"};
constexpr double rangesX[3][2]{{0., 6.}, {0., 25.}, {0., 25.}};
constexpr const char* ptLabel[]{"0.8 #leq #it{p}_{T} < 1.6 GeV/#it{c}", "0.6 #leq #it{p}_{T} < 3.0 GeV/#it{c}", "0.6 #leq #it{p}_{T} < 3.0 GeV/#it{c}"};

void EbyeMult(){
  gStyle->SetOptStat(0);
  auto _file0 = TFile::Open("AnalysisResultsLHC18qrSysHyperloop2.root");
  auto fileOut = TFile::Open("multiplicity.root", "recreate");
  TH2* hRec[3];
  hRec[0] = (TH2F *)_file0->Get("antid-lambda-ebye/QA/nRecPerEvAntid");
  hRec[1] = (TH2F *)_file0->Get("antid-lambda-ebye/QA/nRecPerEvAntiL");
  hRec[2] = (TH2F *)_file0->Get("antid-lambda-ebye/QA/nRecPerEvL");
  for (int iP{0}; iP < 3; ++iP) {
    hRec[iP]->RebinX(10);
    for (int iC{1}; iC < 9; ++iC) {
        TCanvas cc(Form("cNRec_%s_%d", PartText[iP], iC), "", 600, 600);
        cc.SetLeftMargin(0.14);
        cc.SetRightMargin(0.02);
        cc.SetTopMargin(0.05);
        TH1F *proj_ = (TH1F*)hRec[iP]->ProjectionY(Form("nRec_%s_%d", PartText[iP], iC), iC, iC);
        TH1F *proj = new TH1F();
        proj_->GetYaxis()->SetTitle(Form("#it{P}(#it{N}_{%s})", PartLabel[iP]));
        proj_->Scale(1./proj_->GetEntries());
        shift_binning(proj_, proj);
        proj->GetXaxis()->SetRangeUser(rangesX[iP][0], rangesX[iP][1]);
        proj->GetYaxis()->SetRangeUser(2.e-9, 4.);
        TF1 pois("pois", "TMath::Poisson(x, [0])", 0, 1000.);
        pois.SetNpx(10000);
        pois.SetParLimits(0, 0, 1000.);
        proj->Fit("pois", "MNL+");
        SetHistStyle(proj);
        fileOut->cd();
        //proj->Write();
        cc.cd();
        cc.SetLogy();
        proj->Draw("pe");
        pois.Draw("same");
        TLatex txt;
        txt.SetNDC(true);
        txt.SetTextFont(44);
        txt.SetTextSize(25);
        txt.DrawLatex(0.55, 0.87, "Pb#minusPb, #sqrt{s_{NN}} = 5.02 TeV");
        txt.DrawLatex(0.55, 0.82, Form("Centrality %.f-%.f%%", hRec[iP]->GetXaxis()->GetBinLowEdge(iC), hRec[iP]->GetXaxis()->GetBinUpEdge(iC)));
        txt.DrawLatex(0.55, 0.77, "|#eta| < 0.8");
        txt.DrawLatex(0.55, 0.72, PartLabel[iP]);
        txt.DrawLatex(0.55, 0.67, ptLabel[iP]);
        TLegend ll(0.6, 0.6, 0.8, 0.65);
        ll.SetTextFont(44);
        ll.SetTextSize(25);
        ll.AddEntry(&pois, "Poisson fit", "l");
        ll.Draw("same");
        cc.Write();
        cc.Print(Form("plots/Mult%s_%d.pdf", PartText[iP], iC));
    }
  }
  fileOut->Close();
}