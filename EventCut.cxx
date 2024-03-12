void EventCut(){
  auto fin = TFile::Open("~/Downloads/AnalysisResultsTot.root");
  auto fo = TFile::Open("eventCut.root", "recreate");
  auto gloVsCent = (TH2F*)fin->Get("antid-lambda-ebye/QA/GlobalMultVsCent");
  auto pvVsCent = (TH2F*)fin->Get("antid-lambda-ebye/QA/PvMultVsCent");
  TH1D gloVsCentCalib("gloVsCentCalib", "Centrality T0C (%);#it{N}_{global}", 100, 0, 100);
  TH1D gloVsCentCalibErr("gloVsCentCalibErr", "Centrality T0C (%);#it{N}_{global}", 100, 0, 100);
  TH1D pvVsCentCalib("pvVsCentCalib", "Centrality T0C (%);#it{N}_{pv}", 100, 0, 100);
  TH1D pvVsCentCalibErr("pvVsCentCalibErr", "Centrality T0C (%);#it{N}_{pv}", 100, 0, 100);
  for (int iB{1}; iB < 89; ++iB) {
    TH1D* projGloVsCent = (TH1D*)gloVsCent->ProjectionY(Form("projGlo_%d", iB), iB, iB);
    TH1D* projPvVsCent = (TH1D*)pvVsCent->ProjectionY(Form("projPv_%d", iB), iB, iB);
    auto r = projGloVsCent->Fit("gaus", "SML+");
    if (r->CovMatrixStatus() != 3) continue;
    auto f = projGloVsCent->GetFunction("gaus");
    f->SetName("gaus1");
    auto mu1 = f->GetParameter(1);
    auto sig1 = f->GetParameter(2);
    projGloVsCent->Fit("gaus", "RML+", "", mu1 - .5 * sig1, mu1 + 1.5 * sig1);

    f = projGloVsCent->GetFunction("gaus");
    auto mu = f->GetParameter(1);
    auto muErr = f->GetParError(1);
    auto sig = f->GetParameter(2);
    auto sigErr = f->GetParError(2);
    gloVsCentCalib.SetBinContent(iB, mu);
    gloVsCentCalib.SetBinError(iB, muErr);
    gloVsCentCalibErr.SetBinContent(iB, sig);
    gloVsCentCalibErr.SetBinError(iB, sigErr);
    fo->cd();
    projGloVsCent->Write();
    projPvVsCent->Write();

    r = projPvVsCent->Fit("gaus", "SML+");
    if (r->CovMatrixStatus() != 3) continue;
    f = projPvVsCent->GetFunction("gaus");
    f->SetName("gaus1");
    mu1 = f->GetParameter(1);
    sig1 = f->GetParameter(2);
    projPvVsCent->Fit("gaus", "RML+", "", mu1 - .5 * sig1, mu1 + 1.5 * sig1);
    f = projPvVsCent->GetFunction("gaus");
    mu = f->GetParameter(1);
    muErr = f->GetParError(1);
    sig = f->GetParameter(2);
    sigErr = f->GetParError(2);
    pvVsCentCalib.SetBinContent(iB, mu);
    pvVsCentCalib.SetBinError(iB, muErr);
    pvVsCentCalibErr.SetBinContent(iB, sig);
    pvVsCentCalibErr.SetBinError(iB, sigErr);
    fo->cd();
    projPvVsCent->Write();
  }
  fo->cd();

  TF1 fgloMu("fgloMu", "pol3");
  TF1 fgloSig("fgloSig", "pol2");
  TF1 fpvMu("fpvMu", "pol3");
  TF1 fpvSig("fpvSig", "pol2");

  gloVsCentCalib.Fit("fgloMu");
  gloVsCentCalibErr.Fit("fgloSig");
  pvVsCentCalib.Fit("fpvMu");
  pvVsCentCalibErr.Fit("fpvSig");

  gloVsCentCalib.Write();
  gloVsCentCalibErr.Write();
  pvVsCentCalib.Write();
  pvVsCentCalibErr.Write();

  TF1 fCutUp("fCutUp", "pol3(0) + 3 * pol2(4)", 4);
  fCutUp.SetParameter(0, fgloMu.GetParameter(0));
  fCutUp.SetParameter(1, fgloMu.GetParameter(1));
  fCutUp.SetParameter(2, fgloMu.GetParameter(2));
  fCutUp.SetParameter(3, fgloMu.GetParameter(3));
  fCutUp.SetParameter(4, fgloSig.GetParameter(0));
  fCutUp.SetParameter(5, fgloSig.GetParameter(1));
  fCutUp.SetParameter(6, fgloSig.GetParameter(2));

  TF1 fCutDown("fCutDown", "pol3(0) - 3 * pol2(4)", 4);
  fCutDown.SetParameter(0, fgloMu.GetParameter(0));
  fCutDown.SetParameter(1, fgloMu.GetParameter(1));
  fCutDown.SetParameter(2, fgloMu.GetParameter(2));
  fCutDown.SetParameter(3, fgloMu.GetParameter(3));
  fCutDown.SetParameter(4, fgloSig.GetParameter(0));
  fCutDown.SetParameter(5, fgloSig.GetParameter(1));
  fCutDown.SetParameter(6, fgloSig.GetParameter(2));

  TF1 fCutPvUp("fCutPvUp", "pol3(0) + 3 * pol2(4)", 4);
  fCutPvUp.SetParameter(0, fpvMu.GetParameter(0));
  fCutPvUp.SetParameter(1, fpvMu.GetParameter(1));
  fCutPvUp.SetParameter(2, fpvMu.GetParameter(2));
  fCutPvUp.SetParameter(3, fpvMu.GetParameter(3));
  fCutPvUp.SetParameter(4, fpvSig.GetParameter(0));
  fCutPvUp.SetParameter(5, fpvSig.GetParameter(1));
  fCutPvUp.SetParameter(6, fpvSig.GetParameter(2));

  TF1 fCutPvDown("fCutPvDown", "pol3(0) - 3 * pol2(4)", 4);
  fCutPvDown.SetParameter(0, fpvMu.GetParameter(0));
  fCutPvDown.SetParameter(1, fpvMu.GetParameter(1));
  fCutPvDown.SetParameter(2, fpvMu.GetParameter(2));
  fCutPvDown.SetParameter(3, fpvMu.GetParameter(3));
  fCutPvDown.SetParameter(4, fpvSig.GetParameter(0));
  fCutPvDown.SetParameter(5, fpvSig.GetParameter(1));
  fCutPvDown.SetParameter(6, fpvSig.GetParameter(2));

  fCutUp.Write();
  fCutDown.Write();
  fCutPvUp.Write();
  fCutPvDown.Write();
  fo->Close();
  fin->Close();
}