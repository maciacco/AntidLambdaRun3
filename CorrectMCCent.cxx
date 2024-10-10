#include "Constants.h"

void CorrectMCCent(const char* ifile = "TList_1715240052757"){
  TFile *_file0 = TFile::Open(Form("%s.root", ifile));
  TList* ccdb_object = (TList*)_file0->Get("ccdb_object");
  TH1D* v0mCalib = (TH1D*)ccdb_object->FindObject("hMultSelCalib_V0M");
  TGraph g;
  //std::cout << v0mCalib->GetBinContent(1) << std::endl;
  for (int iX{1}; iX <= v0mCalib->GetNbinsX(); ++iX) {
    int centMeasd = v0mCalib->GetBinContent(iX);
    double X = v0mCalib->GetXaxis()->GetBinLowEdge(iX);
    double newX = std::pow(((calib_param[0]+calib_param[1]*std::pow(X,calib_param[2]))-calib_param[3])/calib_param[4],1.0/calib_param[5]);
    int centTruth = v0mCalib->GetBinContent(v0mCalib->FindBin(newX));
    g.AddPoint(centMeasd, centTruth);
    if ((centTruth % 10) == 0) {std::cout << centMeasd << std::endl;}
  }
  TCanvas c;
  c.cd();
  g.Draw();
  c.Print("CorrectCalib.pdf");

  TFile out("CorrectCalib.root", "recreate");
  out.cd();
  g.Write();
  out.Close();
}