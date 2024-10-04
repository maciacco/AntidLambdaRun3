#include "Constants.h"

// ported from o2
constexpr double MassProton = 0.9382721;
constexpr double MassDeuteron = 1.87561294257;
Double_t alephBetheBloch(Double_t *x, Double_t *par)
{
   Float_t bg =x[0];
   Double_t beta = bg / std::sqrt(1. + bg * bg);
   Double_t aa = std::pow(beta, par[3]);
   Double_t bb = std::pow(1. / bg, par[4]);
   bb = std::log(par[2] + bb);
   Double_t f = (par[1] - aa - bb) * par[0] / aa;
   return f;
}

void BetheBlochCalib(const char *period = "calibMc", const double low = 0.7, const double up = 4, const double low2 = 0.45, const double lowd = 0.4, const double upd = 1.){
  auto inFile = TFile::Open(Form("%s/AnalysisResults_%s.root", kDataDir, period));
  auto outFile = TFile::Open(Form("%s/out_file_%s.root", kCalibDir, period), "recreate");
  outFile->mkdir("fits");
  TH2F* hTPC = dynamic_cast<TH2F*>(inFile->Get("ebye-maker/QA/tpcSignalPr"));
  TH2F* hTPCFull = dynamic_cast<TH2F*>(inFile->Get("ebye-maker/QA/tpcSignal"));
  TGraphErrors gFull;
  // hTPC->RebinX(5);
  int xmin = hTPC->GetXaxis()->FindBin(low);
  int xmax = hTPC->GetXaxis()->FindBin(up);
  int xmin2 = hTPCFull->GetXaxis()->FindBin(low2);
  int xmind = hTPCFull->GetXaxis()->FindBin(lowd);
  int xmaxd = hTPCFull->GetXaxis()->FindBin(upd);
  // std::cout << xmax - xmin << std::endl;
  TH1D* proj = new TH1D[xmax - xmin];
  TH1D hTPCPr("hTPCPr", ";#beta#gamma;d#it{E}/d#it{x}_{TPC} (a.u.)", xmax - 1, 0., up / MassProton);
  TH1D hTPCDe("hTPCDe", ";#beta#gamma;d#it{E}/d#it{x}_{TPC} (a.u.)", xmax - 1, 0., up / MassDeuteron);
  for (int iP{xmin}; iP < xmax; ++iP) {
    // if (hTPC->GetXaxis()->GetBinCenter(iP) > 3. && hTPC->GetXaxis()->GetBinCenter(iP) < 4.5) continue;
    TH1D* proj_tmp = dynamic_cast<TH1D*>(hTPC->ProjectionY(Form("h_%d", iP), iP, iP));
    proj_tmp->Copy(proj[iP - xmin]);
    proj[iP - xmin].GetXaxis()->SetRangeUser(proj[iP - xmin].GetMean() - proj[iP - xmin].GetStdDev(), proj[iP - xmin].GetMean() + proj[iP - xmin].GetStdDev());
    double mean = proj[iP - xmin].GetBinCenter(proj[iP - xmin].GetMaximumBin());
    proj[iP - xmin].GetXaxis()->SetRangeUser(mean - 0.12 * mean, mean + 0.12 * mean);
    double rms = proj[iP - xmin].GetRMS();
    proj[iP - xmin].GetXaxis()->SetRangeUser(0, 1.e4);
    proj[iP - xmin].Fit("gaus", "MRQL+", "", mean - 2.5 * rms, mean + 2.5 * rms);
    auto fitFunc = proj[iP - xmin].GetFunction("gaus");
    outFile->cd("fits");
    proj[iP - xmin].Write();
    hTPCPr.SetBinContent(iP, fitFunc->GetParameter(1));
    hTPCPr.SetBinError(iP, fitFunc->GetParameter(2));
    gFull.AddPoint(hTPCPr.GetBinCenter(iP), hTPCPr.GetBinContent(iP));
    gFull.SetPointError(gFull.GetN() - 1, 0., hTPCPr.GetBinError(iP));
  }
  TF1 alephBB("alephBB", alephBetheBloch, 0.1, 10., 5);
  hTPCPr.Fit("alephBB", "RM+", "", low, up);
  TF1 alephBBcpy("alephBBcpy", alephBetheBloch, 0.1, 10., 5);
  alephBBcpy.SetParameter(0, alephBB.GetParameter(0));
  alephBBcpy.SetParameter(1, alephBB.GetParameter(1));
  alephBBcpy.SetParameter(2, alephBB.GetParameter(2));
  alephBBcpy.SetParameter(3, alephBB.GetParameter(3));
  alephBBcpy.SetParameter(4, alephBB.GetParameter(4));

  TH1D* proj2 = new TH1D[xmin - xmin2];
  for (int iP{xmin2}; iP < xmin; ++iP) {
    // if (hTPC->GetXaxis()->GetBinCenter(iP) > 3. && hTPC->GetXaxis()->GetBinCenter(iP) < 4.5) continue;
    TH1D* proj_tmp = dynamic_cast<TH1D*>(hTPCFull->ProjectionY(Form("h_%d", iP), iP, iP));
    proj_tmp->Copy(proj2[iP - xmin2]);
    double val = alephBBcpy.Eval(hTPCPr.GetXaxis()->GetBinCenter(iP));
    proj2[iP - xmin2].GetXaxis()->SetRangeUser(val - 20., val + 20.);
    // std::cout << val << std::endl;
    double mean = proj2[iP - xmin2].GetBinCenter(proj2[iP - xmin2].GetMaximumBin());
    proj2[iP - xmin2].GetXaxis()->SetRangeUser(mean - 20., mean + 20.);
    // std::cout << mean << std::endl;
    double rms = proj2[iP - xmin2].GetRMS();
    proj2[iP - xmin2].GetXaxis()->SetRangeUser(0, 1.e4);
    proj2[iP - xmin2].Fit("gaus", "MRQL+", "", mean - 2.5 * rms, mean + 2.5 * rms);
    auto fitFunc = proj2[iP - xmin2].GetFunction("gaus");
    outFile->cd("fits");
    proj2[iP - xmin2].Write();
    hTPCPr.SetBinContent(iP, fitFunc->GetParameter(1));
    hTPCPr.SetBinError(iP, fitFunc->GetParameter(2));
    gFull.AddPoint(hTPCPr.GetBinCenter(iP), hTPCPr.GetBinContent(iP));
    gFull.SetPointError(gFull.GetN() - 1, 0., hTPCPr.GetBinError(iP));
  }
  hTPCPr.Fit("alephBB", "RM+", "", low2, up);


  TH1D* projd = new TH1D[xmaxd - xmind];
  for (int iP{xmaxd - 1}; iP >= xmind; --iP) {
    // if (hTPC->GetXaxis()->GetBinCenter(iP) > 3. && hTPC->GetXaxis()->GetBinCenter(iP) < 4.5) continue;
    TH1D* proj_tmp = dynamic_cast<TH1D*>(hTPCFull->ProjectionY(Form("h_%d", iP), iP, iP));
    proj_tmp->Copy(projd[iP - xmind]);
    double val = alephBB.Eval(hTPCDe.GetXaxis()->GetBinCenter(iP));
    projd[iP - xmind].GetXaxis()->SetRangeUser(val - 20, val + 100);
    // std::cout << val << std::endl;
    double mean = projd[iP - xmind].GetBinCenter(projd[iP - xmind].GetMaximumBin());
    projd[iP - xmind].GetXaxis()->SetRangeUser(mean - 20., mean + 100.);
    // std::cout << mean << std::endl;
    double rms = projd[iP - xmind].GetRMS();
    projd[iP - xmind].GetXaxis()->SetRangeUser(0, 1.e4);
    projd[iP - xmind].Fit("gaus", "MRQL+", "", mean - 2.5 * rms, mean + 2.5 * rms);
    auto fitFunc = projd[iP - xmind].GetFunction("gaus");
    outFile->cd("fits");
    projd[iP - xmind].Write();
    hTPCDe.SetBinContent(iP, fitFunc->GetParameter(1));
    hTPCDe.SetBinError(iP, fitFunc->GetParameter(2));
    gFull.AddPoint(hTPCDe.GetBinCenter(iP), hTPCDe.GetBinContent(iP));
    gFull.SetPointError(gFull.GetN() - 1, 0., hTPCDe.GetBinError(iP));
    gFull.Fit("alephBB", "NRM+", "", hTPCDe.GetBinCenter(iP), 5.);
  }
  // hTPCDe.Fit("alephBB", "RM+", "", low2, up);

  TF1 alephBB2("alephBB2", alephBetheBloch, 0.1, 10., 5);
  gFull.Fit("alephBB2", "RM+", "", 0.2, 5.);
  outFile->cd();
  hTPCPr.Write();
  hTPCDe.Write();
  gFull.Write();
  alephBB.Write();
  outFile->Close();
  inFile->Close();
}