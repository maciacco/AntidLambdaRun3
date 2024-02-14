void test(){
  auto fIn = TFile::Open("~/Downloads/AnalysisResults_0.root");
  auto fOut = TFile::Open("out_test_h.root", "recreate");
  auto ev = dynamic_cast<THnSparse*>(fIn->Get("antid-lambda-ebye/nEv"));
  auto sp = dynamic_cast<THnSparse*>(fIn->Get("antid-lambda-ebye/nL"));
  auto sq = dynamic_cast<THnSparse*>(fIn->Get("antid-lambda-ebye/nSqL"));
  sp->GetAxis(2)->SetRange(4, 4);
  sq->GetAxis(2)->SetRange(4, 4);
  auto hev = dynamic_cast<TH1D*>(ev->Projection(1));
  auto h = dynamic_cast<TH1D*>(sp->Projection(1));
  auto hq = dynamic_cast<TH1D*>(sq->Projection(1));
  auto hh = dynamic_cast<TH1D*>(fIn->Get("antid-lambda-ebye/q1L"));
  auto hhsq = dynamic_cast<TH1D*>(fIn->Get("antid-lambda-ebye/q1sqL"));
  for (int i{1}; i < hhsq->GetNbinsX() + 1; ++i) {
    auto n = h->GetBinContent(i);
    auto nsq = hq->GetBinContent(i);
    auto nev = hev->GetBinContent(i);
    hh->SetBinContent(i, nev > 0 ? n / nev : 0.);
    h->SetBinContent(i, nev > 0 && n > 0 ? (nsq / nev - std::pow(n / nev , 2) - n / nev) / std::pow(n / nev, 2.) : 0.);
    h->SetBinError(i, 0.);
    hh->SetBinError(i, 0.);
  }
  fOut->cd();
  h->Write();
  hev->Write();
  hh->Write();
  fOut->Close();
  fIn->Close();
}