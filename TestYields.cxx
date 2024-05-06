Double_t myfunction(Double_t *x, Double_t *par)
{
  Float_t xx = x[0];
  if (xx > par[0] && xx < par[1])
  {
    return 0.;
  }
  Double_t f = TMath::Exp(par[2] + par[3] * xx);
  return f;
}

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

constexpr double invMassCut = 0.0013;
constexpr double nSigmaTPCCutLow = 2.;
constexpr double nSigmaTPCCut = 2.;
constexpr double massTOFCut = .085;
constexpr bool kLambda = true;
constexpr bool kDeuteron = true;

int colors[] = {TColor::GetColor("#ff3300"), TColor::GetColor("#ec6e0a"), TColor::GetColor("#daaa14"), TColor::GetColor("#c7e51e"), TColor::GetColor("#85dd69"), TColor::GetColor("#42d6b4"), TColor::GetColor("#00ceff"), TColor::GetColor("#009adf"), TColor::GetColor("#0067c0"), TColor::GetColor("#595959"), TColor::GetColor("#0033a1")};

void TestYields()
{
  gStyle->SetOptStat(0);
  auto _file0 = TFile::Open("AnalysisResults_LHC15o_20240501_cpaTests.root");
  auto fileOut = TFile::Open("testYields_cpa4_.root", "recreate");
  auto file_eff = TFile::Open("efficiency__cpa4_.root");
  TH1D k1("k1", ";Centrality (%);Yield", 10, 0, 100);
  if (kLambda)
  {
    auto hMass = (TH3F *)_file0->Get("antid-lambda-ebye_cpa4/QA/massLambda");
    hMass->RebinX(10);
    hMass->RebinY(2);
    TCanvas ctot("PurityLambda", "purityLambda", 600, 600);
    TH1D *res[9]{nullptr};
    TLegend leg(0.3, 0.2, 0.7, 0.5);
    for (int iC{1}; iC < 9; ++iC)
    {
      double bin[]{1., 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2., 2.2, 2.4, 2.6, 2.8, 3., 3.4, 4};
      res[iC - 1] = new TH1D(Form("res_%d", iC), ";#it{p}_{T} (GeV/#it{c});Purity", 17, bin);
      auto effL = (TH1D*)file_eff->Get(Form("effL_%d", iC));
      double totYield = 0.;
      double totYieldErr = 0.;
      for (int iP{1}; iP <= 17; ++iP)
      {
        auto eff_pt = effL->GetBinContent(iP);
        auto binLeft = hMass->GetYaxis()->FindBin(effL->GetXaxis()->GetBinLowEdge(iP) + 0.001);
        auto binRight = hMass->GetYaxis()->FindBin(effL->GetXaxis()->GetBinUpEdge(iP) - 0.001);
        auto proj = (TH1F *)hMass->ProjectionZ(Form("mass_%.1f_%.1f", hMass->GetYaxis()->GetBinLowEdge(iP), hMass->GetYaxis()->GetBinUpEdge(iP)), iC, iC, binLeft, binRight);
        proj->SetTitle(";#it{M}(p + #pi^{-}) (GeV/#it{c});Entries");
        // auto effAntiL = (TH1D*)file_eff->Get(Form("effAntiL_%d", iC));
        // effAntiL->SetName("effAntiL");
        TH1F projCpy(*proj);
        for (int iB{proj->FindBin(1.109)}; iB < proj->FindBin(1.123); ++iB)
        {
          projCpy.SetBinContent(iB, 0.);
        }
        TF1 ff("ff", "pol3", 1.1, 1.125);
        projCpy.Fit("ff", "QRM+", "", 1.15, 1.13);
        fileOut->cd();
        TF1 f("f", "pol3", 1.102, 1.132);
        f.SetParameter(0, ff.GetParameter(0));
        f.SetParameter(1, ff.GetParameter(1));
        f.SetParameter(2, ff.GetParameter(2));
        f.SetParameter(3, ff.GetParameter(3));
        auto bkg = f.Integral(1.1157 - invMassCut, 1.1157 + invMassCut) / proj->GetBinWidth(1);
        auto sig = proj->Integral(proj->FindBin(1.1157 - invMassCut), proj->FindBin(1.1157 + invMassCut));
        auto deltaPt = effL->GetXaxis()->GetBinUpEdge(iP) - effL->GetXaxis()->GetBinLowEdge(iP);
        std::cout << "Purity = " << (sig - bkg) << std::endl;
        res[iC - 1]->SetBinContent(iP, (sig - bkg) / eff_pt);
        res[iC - 1]->SetBinError(iP, effL->GetBinError(iP) * (sig - bkg) / eff_pt / eff_pt);
        totYield += (sig - bkg) / eff_pt;
        totYieldErr += pow(effL->GetBinError(iP) * (sig - bkg) / eff_pt / eff_pt, 2.);
        // proj->Write();
        // projCpy.Write();
        TCanvas c(Form("pT_%.1f_%.1f", hMass->GetYaxis()->GetBinLowEdge(iP), hMass->GetYaxis()->GetBinUpEdge(iP)), "c", 600, 600);
        c.SetLeftMargin(0.14);
        c.SetRightMargin(0.02);
        c.SetTopMargin(0.05);
        c.cd();
        TLatex txt;
        proj->GetXaxis()->SetRangeUser(1.096, 1.134);
        proj->GetYaxis()->SetRangeUser(0, proj->GetBinContent(proj->GetMaximumBin()) * 1.3);
        txt.SetNDC(true);
        txt.SetTextFont(44);
        txt.SetTextSize(21);
        SetHistStyle(proj);
        proj->Draw("pe");
        f.SetLineWidth(3);
        f.Draw("same");
        TLine ll(1.1157 - invMassCut, 0, 1.1157 - invMassCut, proj->GetBinContent(proj->GetMaximumBin()) * 0.6);
        TLine lr(1.1157 + invMassCut, 0, 1.1157 + invMassCut, proj->GetBinContent(proj->GetMaximumBin()) * 0.6);
        ll.SetLineStyle(kDashed);
        lr.SetLineStyle(kDashed);
        ll.Draw("same");
        lr.Draw("same");
        txt.DrawLatex(0.18, 0.67, Form("%.1f #leq #it{p}_{T} < %.1f GeV/#it{c}", hMass->GetYaxis()->GetBinLowEdge(iP), hMass->GetYaxis()->GetBinUpEdge(iP)));
        txt.DrawLatex(0.18, 0.77, "#Lambda+#bar{#Lambda}");
        txt.DrawLatex(0.18, 0.87, "Pb#minusPb, #sqrt{s_{NN}} = 5.02 TeV");
        txt.DrawLatex(0.18, 0.82, Form("Centrality %.f-%.f%%", hMass->GetXaxis()->GetBinLowEdge(iC), hMass->GetXaxis()->GetBinUpEdge(iC)));
        txt.DrawLatex(0.18, 0.72, "|#eta| < 0.8");
        txt.DrawLatex(0.62, 0.87, Form("#frac{S}{S + B} = %.3f #pm %.3f", (sig - bkg) / sig, (1. / sig - (sig - bkg) / sig / sig) * std::sqrt(sig)));
        fileOut->cd();
        c.Write();
        //c.Print(Form("plots/LambdaPurity_cent%d_pt_%.1f_%.1f.pdf", iC, hMass->GetYaxis()->GetBinLowEdge(iP), hMass->GetYaxis()->GetBinUpEdge(iP)));
      }
      fileOut->cd();
      res[iC - 1]->SetMarkerStyle(20);
      res[iC - 1]->SetMarkerSize(1.5);
      res[iC - 1]->SetMarkerColor(colors[iC - 1]);
      res[iC - 1]->SetLineColor(colors[iC - 1]);
      res[iC - 1]->SetLineWidth(2);
      res[iC - 1]->Write();
      leg.AddEntry(res[iC - 1], Form("%.0f-%.0f%%", hMass->GetXaxis()->GetBinLowEdge(iC), hMass->GetXaxis()->GetBinUpEdge(iC)));
      k1.SetBinContent(iC, totYield);
      k1.SetBinError(iC, std::sqrt(totYieldErr));
    }
    for (int iC{0}; iC < 8; ++iC)
    {
      ctot.cd();
      //res[iC]->GetYaxis()->SetRangeUser(0., 1.1);
      res[iC]->Draw(iC == 0 ? "pe" : "pesame");
    }
    ctot.cd();
    leg.SetTextFont(44);
    leg.SetTextSize(25);
    leg.SetHeader("#Lambda + #bar{#Lambda}");
    leg.SetNColumns(2);
    leg.Draw("same");
    fileOut->cd();
    ctot.Write();
  }
  k1.Write();
  fileOut->Close();
}
