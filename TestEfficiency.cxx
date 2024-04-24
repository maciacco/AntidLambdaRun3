void TestEfficiency(){
  auto _file0 = TFile::Open("~/AnalysisResultsLHC21k2.root");
  auto fileOut = TFile::Open("efficiencyP.root", "recreate");
  auto hRec = (TH3F*)_file0->Get("antid-lambda-ebye/recAntip");
  auto hRecM = (TH3F*)_file0->Get("antid-lambda-ebye/recP");
  auto hGen = (TH3F*)_file0->Get("antid-lambda-ebye/genAntip");
  auto hGenM = (TH3F*)_file0->Get("antid-lambda-ebye/genP");
  //hRec->Add(hRecM);
  //hGen->Add(hGenM);
  // lambda
  //int centBins[10][2]{{0, 20}, {0, 20}, {20, 50}, {20, 50}, {20, 50}, {50, 90}, {50, 90}, {50, 90}, {50, 90}};
  // deuteron
  int centBins[10][2]{{0, 20}, {0, 20}, {20, 40}, {20, 40}, {40, 60}, {40, 60}, {60, 80}, {60, 80}, {80, 90}};
  for (int iC{1}; iC <= 9; ++iC) {
    auto projNum = (TH1F*)hRec->ProjectionY(Form("eff_%d", iC), centBins[iC - 1][0] + 1, centBins[iC - 1][1], 1, 1);
    TH1F projNumCpy(*projNum);
    projNumCpy.SetName(Form("projNum_%d", iC));
    auto projDen = (TH1F*)hGen->ProjectionY(Form("projDen_%d", iC), centBins[iC - 1][0] + 1, centBins[iC - 1][1], 1, 1);
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
