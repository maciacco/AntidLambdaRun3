void TestEfficiency(){
  auto _file0 = TFile::Open("~/AnalysisResultsLHC21k2b.root");
  auto fileOut = TFile::Open("efficiencyLb.root", "recreate");
  auto hRec = (TH3F*)_file0->Get("antid-lambda-ebye/recAntiL");
  auto hRecM = (TH3F*)_file0->Get("antid-lambda-ebye/recL");
  auto hGen = (TH3F*)_file0->Get("antid-lambda-ebye/genAntiL");
  auto hGenM = (TH3F*)_file0->Get("antid-lambda-ebye/genL");
  hRec->Add(hRecM);
  hGen->Add(hGenM);
  for (int iC{1}; iC <= 1; ++iC) {
    auto projNum = (TH1F*)hRec->ProjectionY(Form("eff_%d", iC), 1, 100, 1, 1);
    TH1F projNumCpy(*projNum);
    projNumCpy.SetName(Form("projNum_%d", iC));
    auto projDen = (TH1F*)hGen->ProjectionY(Form("projDen_%d", iC), 1, 100, 1, 1);
    // projNum->Rebin(2);
    // projDen->Rebin(2);
    projNum->Divide(projNum, projDen, 1., 1., "B");
    fileOut->cd();
    projDen->Write();
    projNumCpy.Write();
    projNum->Write();
  }
  fileOut->Close();
}
