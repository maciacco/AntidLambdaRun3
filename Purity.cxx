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

constexpr double invMassCut = 0.0025;
constexpr double nSigmaTPCCutLow = 2.;
constexpr double nSigmaTPCCut = 2.;
constexpr double massTOFCut = .085;
constexpr bool kLambda = true;
constexpr bool kDeuteron = true;

int colors[] = {TColor::GetColor("#ff3300"), TColor::GetColor("#ec6e0a"), TColor::GetColor("#daaa14"), TColor::GetColor("#c7e51e"), TColor::GetColor("#85dd69"), TColor::GetColor("#42d6b4"), TColor::GetColor("#00ceff"), TColor::GetColor("#009adf"), TColor::GetColor("#0067c0"), TColor::GetColor("#595959"), TColor::GetColor("#0033a1")};

void Purity()
{
  gStyle->SetOptStat(0);
  auto _file0 = TFile::Open("./AnalysisResultsLHC15oGridJob.root");
  auto fileOut = TFile::Open("purity.root", "recreate");
  if (kLambda)
  {
    auto hMass = (TH3F *)_file0->Get("antid-lambda-ebye/QA/massLambda");
    hMass->RebinX(10);
    hMass->RebinY(2);
    hMass->RebinZ(2);
    TCanvas ctot("PurityLambda", "purityLambda", 600, 600);
    TH1D *res[9]{nullptr};
    TLegend leg(0.3, 0.2, 0.7, 0.5);
    for (int iC{1}; iC < 10; ++iC)
    {
      res[iC - 1] = new TH1D(Form("res_%d", iC), ";#it{p}_{T} (GeV/#it{c});Purity", 20, 0., 4.);
      for (int iP{6}; iP <= 20; ++iP)
      {
        auto proj = (TH1F *)hMass->ProjectionZ(Form("mass_%.1f_%.1f", hMass->GetYaxis()->GetBinLowEdge(iP), hMass->GetYaxis()->GetBinUpEdge(iP)), iC, iC, iP, iP);
        proj->SetTitle(";#it{M}(p + #pi^{-}) (GeV/#it{c});Entries");
        TH1F projCpy(*proj);
        for (int iB{proj->FindBin(1.109)}; iB < proj->FindBin(1.123); ++iB)
        {
          projCpy.SetBinContent(iB, 0.);
        }
        TF1 ff("ff", "expo", 1.09, 1.14);
        projCpy.Fit("ff", "QRM+", "", 1.102, 1.132);
        fileOut->cd();
        TF1 f("f", "expo", 1.102, 1.132);
        f.SetParameter(0, ff.GetParameter(0));
        f.SetParameter(1, ff.GetParameter(1));
        // f.SetParameter(2, ff.GetParameter(2));
        auto bkg = f.Integral(1.1157 - invMassCut, 1.1157 + invMassCut) / proj->GetBinWidth(1);
        auto sig = proj->Integral(proj->FindBin(1.1157 - invMassCut), proj->FindBin(1.1157 + invMassCut));
        std::cout << "Purity = " << (sig - bkg) / sig << std::endl;
        res[iC - 1]->SetBinContent(iP, (sig - bkg) / sig);
        res[iC - 1]->SetBinError(iP, (1. / sig - (sig - bkg) / sig / sig) * std::sqrt(sig));
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
        c.Print(Form("plots/LambdaPurity_cent%d_pt_%.1f_%.1f.pdf", iC, hMass->GetYaxis()->GetBinLowEdge(iP), hMass->GetYaxis()->GetBinUpEdge(iP)));
      }
      fileOut->cd();
      res[iC - 1]->SetMarkerStyle(20);
      res[iC - 1]->SetMarkerSize(1.5);
      res[iC - 1]->SetMarkerColor(colors[iC - 1]);
      res[iC - 1]->SetLineColor(colors[iC - 1]);
      res[iC - 1]->SetLineWidth(2);
      res[iC - 1]->Write();
      leg.AddEntry(res[iC - 1], Form("%.0f-%.0f%%", hMass->GetXaxis()->GetBinLowEdge(iC), hMass->GetXaxis()->GetBinUpEdge(iC)));
    }
    for (int iC{0}; iC < 9; ++iC)
    {
      ctot.cd();
      res[iC]->GetYaxis()->SetRangeUser(0., 1.1);
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
    ctot.Print("purityLambda.pdf");
  }

  if (kDeuteron)
  {
    auto nSigmaTPC = (TH3F *)_file0->Get("antid-lambda-ebye/QA/tpcNsigmaGlo_d");
    auto massTOF = (TH3F *)_file0->Get("antid-lambda-ebye/QA/tofMass_d");
    nSigmaTPC->RebinX(10);
    //nSigmaTPC->RebinY(4);
    nSigmaTPC->RebinZ(2);
    massTOF->RebinX(10);
    //massTOF->RebinY(4);
    massTOF->RebinZ(2);
    TCanvas ctot("PurityAntideuteron", "purityAntideuteron", 600, 600);
    TH1D *res[10]{nullptr};
    TLegend leg(0.3, 0.2, 0.7, 0.5);
    for (int iC{1}; iC < 10; ++iC)
    {
      res[iC - 1] = new TH1D(Form("res_%d", iC), ";#it{p}_{T} (GeV/#it{c});Purity", 20, 0., 2.);
      for (int iP{7}; iP < 11; ++iP)
      {
        auto proj = (TH1F *)nSigmaTPC->ProjectionZ(Form("nSigmaTpc_%.1f_%.1f", nSigmaTPC->GetYaxis()->GetBinLowEdge(iP), nSigmaTPC->GetYaxis()->GetBinUpEdge(iP)), iC, iC, iP, iP);
        proj->SetTitle(";n#sigma_{TPC};Entries");
        TH1F projCpy(*proj);
        for (int iB{proj->FindBin(-4.)}; iB < proj->FindBin(5); ++iB)
        {
          projCpy.SetBinContent(iB, 0.);
        }
        TF1 ff("ff", "gaus", -5, 5);
        ff.SetParLimits(1, -10., -4);
        ff.SetParLimits(2, 0., 1.5);
        projCpy.Fit("ff", "QRM+", "", -5, 5);
        fileOut->cd();
        TF1 f("f", "gaus", -5, 5);
        f.SetParameter(0, ff.GetParameter(0));
        f.SetParameter(1, ff.GetParameter(1));
        f.SetParameter(2, ff.GetParameter(2));
        // f.SetParameter(2, ff.GetParameter(2));
        auto bkg = f.Integral(-nSigmaTPCCut, nSigmaTPCCut) / proj->GetBinWidth(1);
        auto sig = proj->Integral(proj->FindBin(-nSigmaTPCCut), proj->FindBin(nSigmaTPCCut));
        std::cout << "Purity = " << (sig - bkg) / sig << std::endl;
        res[iC - 1]->SetBinContent(iP, (sig - bkg) / sig);
        res[iC - 1]->SetBinError(iP, (1. / sig - (sig - bkg) / sig / sig) * std::sqrt(sig));
        // proj->Write();
        // projCpy.Write();
        TCanvas c(Form("pT_%.1f_%.1f", nSigmaTPC->GetYaxis()->GetBinLowEdge(iP), nSigmaTPC->GetYaxis()->GetBinUpEdge(iP)), "c", 600, 600);
        c.SetLeftMargin(0.14);
        c.SetRightMargin(0.02);
        c.SetTopMargin(0.05);
        c.cd();
        TLatex txt;
        proj->GetXaxis()->SetRangeUser(-5, 5);
        proj->GetYaxis()->SetRangeUser(0, proj->GetBinContent(proj->GetMaximumBin()) * 1.3);
        txt.SetNDC(true);
        txt.SetTextFont(44);
        txt.SetTextSize(21);
        SetHistStyle(proj);
        proj->Draw("pe");
        f.SetLineWidth(3);
        f.Draw("same");
        TLine ll(-nSigmaTPCCutLow, 0, -nSigmaTPCCutLow, proj->GetBinContent(proj->GetMaximumBin()) * 0.6);
        TLine lr(nSigmaTPCCut, 0, nSigmaTPCCut, proj->GetBinContent(proj->GetMaximumBin()) * 0.6);
        ll.SetLineStyle(kDashed);
        lr.SetLineStyle(kDashed);
        ll.Draw("same");
        lr.Draw("same");
        txt.DrawLatex(0.18, 0.67, Form("%.1f #leq #it{p}_{T} < %.1f GeV/#it{c}", nSigmaTPC->GetYaxis()->GetBinLowEdge(iP), nSigmaTPC->GetYaxis()->GetBinUpEdge(iP)));
        txt.DrawLatex(0.18, 0.77, "#bar{d}");
        txt.DrawLatex(0.18, 0.87, "Pb#minusPb, #sqrt{s_{NN}} = 5.02 TeV");
        txt.DrawLatex(0.18, 0.82, Form("Centrality %.f-%.f%%", nSigmaTPC->GetXaxis()->GetBinLowEdge(iC), nSigmaTPC->GetXaxis()->GetBinUpEdge(iC)));
        txt.DrawLatex(0.18, 0.72, "|#eta| < 0.8");
        txt.DrawLatex(0.62, 0.87, Form("#frac{S}{S + B} = %.3f #pm %.3f", (sig - bkg) / sig, (1. / sig - (sig - bkg) / sig / sig) * std::sqrt(sig)));
        fileOut->cd();
        c.Write();
        c.Print(Form("plots/AntideuteronPurity_cent%d_pt_%.1f_%.1f.pdf", iC, nSigmaTPC->GetYaxis()->GetBinLowEdge(iP), nSigmaTPC->GetYaxis()->GetBinUpEdge(iP)));
      }
      for (int iP{11}; iP <= 18; ++iP)
      {
        auto proj = (TH1F *)massTOF->ProjectionZ(Form("massTOF_%.1f_%.1f", massTOF->GetYaxis()->GetBinLowEdge(iP), massTOF->GetYaxis()->GetBinUpEdge(iP)), iC, iC, iP, iP);
        proj->SetTitle(";#it{M}_{TOF} (GeV/#it{c}^{2});Entries");
        TH1F projCpy(*proj);
        for (int iB{proj->FindBin(1.57)}; iB < proj->FindBin(2.1); ++iB)
        {
          projCpy.SetBinContent(iB, 0.);
        }
        TF1 ff("ff", myfunction, 1.4, 2.3, 4);
        ff.SetParLimits(2, 0., 1.e6);
        ff.SetParLimits(3, -1.e4, 0.);
        ff.FixParameter(0, proj->GetXaxis()->GetBinLowEdge(proj->FindBin(1.57)));
        ff.FixParameter(1, proj->GetXaxis()->GetBinUpEdge(proj->FindBin(2.1)));
        projCpy.Fit("ff", "QRML+", "", 1.4, 2.3);
        projCpy.Write();
        fileOut->cd();
        TF1 f("f", "expo", 1.4, 2.3);
        f.SetParameter(0, ff.GetParameter(2));
        f.SetParameter(1, ff.GetParameter(3));
        // f.SetParameter(2, ff.GetParameter(2));
        auto bkg = f.Integral(1.87561294257 - massTOFCut, 1.87561294257 + massTOFCut) / proj->GetBinWidth(1);
        auto sig = proj->Integral(proj->FindBin(1.87561294257 - massTOFCut), proj->FindBin(1.87561294257 + massTOFCut));
        std::cout << "Purity = " << (sig - bkg) / sig << std::endl;
        res[iC - 1]->SetBinContent(iP, (sig - bkg) / sig);
        res[iC - 1]->SetBinError(iP, (1. / sig - (sig - bkg) / sig / sig) * std::sqrt(sig));
        // proj->Write();
        // projCpy.Write();
        TCanvas c(Form("pT_%.1f_%.1f", massTOF->GetYaxis()->GetBinLowEdge(iP), massTOF->GetYaxis()->GetBinUpEdge(iP)), "c", 600, 600);
        c.SetLeftMargin(0.14);
        c.SetRightMargin(0.02);
        c.SetTopMargin(0.05);
        c.cd();
        TLatex txt;
        proj->GetYaxis()->SetRangeUser(0, proj->GetBinContent(proj->GetMaximumBin()) * 1.3);
        proj->GetXaxis()->SetRangeUser(1.5, 2.2);
        txt.SetNDC(true);
        txt.SetTextFont(44);
        txt.SetTextSize(21);
        SetHistStyle(proj);
        proj->Draw("pe");
        f.SetLineWidth(3);
        f.Draw("same");
        TLine ll(1.87561294257 - massTOFCut, 0, 1.87561294257 - massTOFCut, proj->GetBinContent(proj->GetMaximumBin()) * 0.6);
        TLine lr(1.87561294257 + massTOFCut, 0, 1.87561294257 + massTOFCut, proj->GetBinContent(proj->GetMaximumBin()) * 0.6);
        ll.SetLineStyle(kDashed);
        lr.SetLineStyle(kDashed);
        ll.Draw("same");
        lr.Draw("same");
        txt.DrawLatex(0.18, 0.67, Form("%.1f #leq #it{p}_{T} < %.1f GeV/#it{c}", massTOF->GetYaxis()->GetBinLowEdge(iP), massTOF->GetYaxis()->GetBinUpEdge(iP)));
        txt.DrawLatex(0.18, 0.77, "#bar{d}");
        txt.DrawLatex(0.18, 0.72, "|#eta| < 0.8");
        txt.DrawLatex(0.18, 0.87, "Pb#minusPb, #sqrt{s_{NN}} = 5.02 TeV");
        txt.DrawLatex(0.18, 0.82, Form("Centrality %.f-%.f%%", massTOF->GetXaxis()->GetBinLowEdge(iC), massTOF->GetXaxis()->GetBinUpEdge(iC)));
        txt.DrawLatex(0.62, 0.87, Form("#frac{S}{S + B} = %.3f #pm %.3f", (sig - bkg) / sig, (1. / sig - (sig - bkg) / sig / sig) * std::sqrt(sig)));
        fileOut->cd();
        c.Write();
        c.Print(Form("plots/AntideuteronPurity_cent%d_pt_%.1f_%.1f.pdf", iC, massTOF->GetYaxis()->GetBinLowEdge(iP), massTOF->GetYaxis()->GetBinUpEdge(iP)));
      }
      fileOut->cd();
      res[iC - 1]->SetMarkerStyle(20);
      res[iC - 1]->SetMarkerSize(1.5);
      res[iC - 1]->SetMarkerColor(colors[iC - 1]);
      res[iC - 1]->SetLineColor(colors[iC - 1]);
      res[iC - 1]->SetLineWidth(2);
      res[iC - 1]->Write();
      leg.AddEntry(res[iC - 1], Form("%.0f-%.0f%%", nSigmaTPC->GetXaxis()->GetBinLowEdge(iC), nSigmaTPC->GetXaxis()->GetBinUpEdge(iC)));
    }
    for (int iC{0}; iC < 9; ++iC)
    {
      ctot.cd();
      res[iC]->GetYaxis()->SetRangeUser(0., 1.1);
      res[iC]->Draw(iC == 0 ? "pe" : "pesame");
    }
    ctot.cd();
    leg.SetTextFont(44);
    leg.SetTextSize(25);
    leg.SetHeader("#bar{d}");
    leg.SetNColumns(2);
    leg.Draw("same");
    fileOut->cd();
    ctot.Write();
    ctot.Print("purityAntideuteron.pdf");
  }
  fileOut->Close();
}
