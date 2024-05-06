#include "Constants.h"
#include <TFile.h>
#include <THnSparse.h>
#include <TH1D.h>
#include <TColor.h>
#include <TStyle.h>
#include <TH3F.h>
#include <TLegend.h>
#include <TCanvas.h>
int colors[] = {TColor::GetColor("#ff3300"), TColor::GetColor("#ec6e0a"), TColor::GetColor("#daaa14"), TColor::GetColor("#c7e51e"), TColor::GetColor("#85dd69"), TColor::GetColor("#42d6b4"), TColor::GetColor("#00ceff"), TColor::GetColor("#009adf"), TColor::GetColor("#0067c0"), TColor::GetColor("#595959"), TColor::GetColor("#0033a1")};

void TestEfficiency(const int cut = 0){
  const char* cutSet = cutSets[cut];
  gStyle->SetOptStat(0);
  auto _file0 = TFile::Open("AnalysisResultsLHC21k2abcGrid.root");
  auto _file1 = TFile::Open("AnalysisResultsLHC21d6abcGrid.root");
  auto fileOut = TFile::Open(Form("efficiency_%s_grid_.root", cutSet), "recreate");

  // lambda
  auto hSpRec = (THnSparse*)_file1->Get(Form("antid-lambda-ebye%s/nAntiL", cutSet));
  auto hSpRecM = (THnSparse*)_file1->Get(Form("antid-lambda-ebye%s/nL", cutSet));
  auto hSpGen = (THnSparse*)_file1->Get(Form("antid-lambda-ebye%s/nGenAntiL", cutSet));
  auto hSpGenM = (THnSparse*)_file1->Get(Form("antid-lambda-ebye%s/nGenL", cutSet));
  // auto hRec = (TH3F*)_file0->Get(Form("antid-lambda-ebye%s/recAntiL", cutSet));
  // auto hRecM = (TH3F*)_file0->Get(Form("antid-lambda-ebye%s/recL", cutSet));
  // auto hGen = (TH3F*)_file0->Get(Form("antid-lambda-ebye%s/genAntiL", cutSet));
  // auto hGenM = (TH3F*)_file0->Get(Form("antid-lambda-ebye%s/genL", cutSet));
  // hRec->Add(hRecM);
  // hGen->Add(hGenM);
  TH1F *projNum[9]{nullptr};
  TLegend lL(0.2, 0.65, 0.5, 0.85);
  lL.SetNColumns(2);
  lL.SetTextFont(44);
  lL.SetTextSize(23);
  TCanvas cL("cL", "cL", 600, 700);
  int centBinsL[10][2]{{0, 10}, {10, 20}, {20, 30}, {30, 40}, {40, 50}, {50, 60}, {60, 70}, {70, 80}, {80, 90}, {90, 100}};//{{0, 20}, {0, 20}, {20, 50}, {20, 50}, {20, 50}, {50, 90}, {50, 90}, {50, 90}, {50, 90}};
  for (int iC{1}; iC <= 8; ++iC) {
    hSpRecM->GetAxis(1)->SetRange(centBinsL[iC - 1][0] + 1, centBinsL[iC - 1][1]);
    hSpGenM->GetAxis(1)->SetRange(centBinsL[iC - 1][0] + 1, centBinsL[iC - 1][1]);
    projNum[iC - 1] = (TH1F*)hSpRecM->Projection(3);
    //projNum[iC - 1] = (TH1F*)hRecM->ProjectionY(Form("eff_%d", iC), centBinsL[iC - 1][0] + 1, centBinsL[iC - 1][1], 1, 1);
    TH1F projNumCpy(*projNum[iC - 1]);
    projNumCpy.SetName(Form("projNum_%d", iC));
    auto projDen = (TH1F*)hSpGenM->Projection(3);
    //auto projDen = (TH1F*)hGenM->ProjectionY(Form("projDen_%d", iC), centBinsL[iC - 1][0] + 1, centBinsL[iC - 1][1], 1, 1);
    // projNum[iC - 1]->Rebin(2);
    // projDen->Rebin(2);
    projNum[iC - 1]->Divide(projNum[iC - 1], projDen, 1., 1., "B");
    projNum[iC - 1]->GetYaxis()->SetTitle("Efficiency");
    fileOut->cd();
    //projDen->Write();
    //projNumCpy.Write();
    projNum[iC - 1]->Write(Form("effL_%d", iC));
    cL.cd();
    if (iC == 1){
      projNum[iC - 1]->GetYaxis()->SetRangeUser(0, 1.1);
    }
    projNum[iC - 1]->SetLineWidth(2);
    projNum[iC - 1]->SetLineColor(colors[iC - 1]);
    projNum[iC - 1]->Draw(iC == 1 ? "pe" : "samepe");
    lL.AddEntry(projNum[iC - 1], Form("%d-%d%%", centBinsL[iC - 1][0], centBinsL[iC - 1][1]));
  }
  cL.cd();
  lL.Draw("same");
  cL.Print("cL.pdf");

  TLegend lAntiL(0.2, 0.65, 0.5, 0.85);
  lAntiL.SetNColumns(2);
  lAntiL.SetTextFont(44);
  lAntiL.SetTextSize(23);
  TCanvas cAntiL("cAntiL", "cAntiL", 600, 700);
  for (int iC{1}; iC <= 8; ++iC) {
    hSpRec->GetAxis(1)->SetRange(centBinsL[iC - 1][0] + 1, centBinsL[iC - 1][1]);
    hSpGen->GetAxis(1)->SetRange(centBinsL[iC - 1][0] + 1, centBinsL[iC - 1][1]);
    projNum[iC - 1] = (TH1F*)hSpRec->Projection(3);
    //projNum[iC - 1] = (TH1F*)hRec->ProjectionY(Form("eff_%d", iC), centBinsL[iC - 1][0] + 1, centBinsL[iC - 1][1], 1, 1);
    TH1F projNumCpy(*projNum[iC - 1]);
    projNumCpy.SetName(Form("projNum_%d", iC));
    auto projDen = (TH1F*)hSpGen->Projection(3);
    //auto projDen = (TH1F*)hGen->ProjectionY(Form("projDen_%d", iC), centBinsL[iC - 1][0] + 1, centBinsL[iC - 1][1], 1, 1);
    // projNum[iC - 1]->Rebin(2);
    // projDen->Rebin(2);
    projNum[iC - 1]->Divide(projNum[iC - 1], projDen, 1., 1., "B");
    projNum[iC - 1]->GetYaxis()->SetTitle("Efficiency");
    fileOut->cd();
    //projDen->Write();
    //projNumCpy.Write();
    projNum[iC - 1]->Write(Form("effAntiL_%d", iC));
    cAntiL.cd();
    if (iC == 1){
      projNum[iC - 1]->GetYaxis()->SetRangeUser(0, 1.1);
    }
    projNum[iC - 1]->SetLineWidth(2);
    projNum[iC - 1]->SetLineColor(colors[iC - 1]);
    projNum[iC - 1]->Draw(iC == 1 ? "pe" : "samepe");
    lAntiL.AddEntry(projNum[iC - 1], Form("%d-%d%%", centBinsL[iC - 1][0], centBinsL[iC - 1][1]));
  }
  cAntiL.cd();
  lAntiL.Draw("same");
  cAntiL.Print("cAntiL.pdf");

  // deuteron
  TLegend lAntiD(0.5, 0.65, 0.8, 0.85);
  lAntiD.SetNColumns(2);
  lAntiD.SetTextFont(44);
  lAntiD.SetTextSize(23);
  TCanvas cAntiD("cAntiD", "cAntiD", 600, 700);
  TH3F* hRec = (TH3F*)_file0->Get(Form("antid-lambda-ebye%s/recAntid", cutSet));
  TH3F* hRecM = (TH3F*)_file0->Get(Form("antid-lambda-ebye%s/recD", cutSet));
  TH3F* hGen = (TH3F*)_file0->Get(Form("antid-lambda-ebye%s/genAntid", cutSet));
  TH3F* hGenM = (TH3F*)_file0->Get(Form("antid-lambda-ebye%s/genD", cutSet));
  int centBinsD[10][2]{{0, 10}, {10, 20}, {20, 30}, {30, 40}, {40, 50}, {50, 60}, {60, 70}, {70, 80}, {80, 90}, {90, 100}};//{{0, 20}, {0, 20}, {20, 40}, {20, 40}, {40, 60}, {40, 60}, {60, 80}, {60, 80}, {80, 90}};
  for (int iC{1}; iC <= 8; ++iC) {
    projNum[iC - 1] = (TH1F*)hRec->ProjectionY(Form("eff_%d", iC), centBinsD[iC - 1][0] + 1, centBinsD[iC - 1][1], 1, 1);
    TH1F projNumCpy(*projNum[iC - 1]);
    projNumCpy.SetName(Form("projNum_%d", iC));
    auto projDen = (TH1F*)hGen->ProjectionY(Form("projDen_%d", iC), centBinsD[iC - 1][0] + 1, centBinsD[iC - 1][1], 1, 1);
    // projNum[iC - 1]->Rebin(2);
    // projDen->Rebin(2);
    projNum[iC - 1]->Divide(projNum[iC - 1], projDen, 1., 1., "B");
    projNum[iC - 1]->GetYaxis()->SetTitle("Efficiency");
    fileOut->cd();
    //projDen->Write();
    //projNumCpy.Write();
    projNum[iC - 1]->Write(Form("effD_%d", iC));
    cAntiD.cd();
    if (iC == 1){
      projNum[iC - 1]->GetYaxis()->SetRangeUser(0, 1.1);
    }
    projNum[iC - 1]->SetLineWidth(2);
    projNum[iC - 1]->SetLineColor(colors[iC - 1]);
    projNum[iC - 1]->Draw(iC == 1 ? "pe" : "samepe");
    lAntiD.AddEntry(projNum[iC - 1], Form("%d-%d%%", centBinsL[iC - 1][0], centBinsL[iC - 1][1]));
  }
  cAntiD.cd();
  lAntiD.Draw("same");
  cAntiD.Print("cAntiD.pdf");

  // proton
  hRec = (TH3F*)_file0->Get(Form("antid-lambda-ebye%s/recAntip", cutSet));
  hRecM = (TH3F*)_file0->Get(Form("antid-lambda-ebye%s/recP", cutSet));
  hGen = (TH3F*)_file0->Get(Form("antid-lambda-ebye%s/genAntip", cutSet));
  hGenM = (TH3F*)_file0->Get(Form("antid-lambda-ebye%s/genP", cutSet));
  int centBinsP[10][2]{{0, 10}, {10, 20}, {20, 30}, {30, 40}, {40, 50}, {50, 60}, {60, 70}, {70, 80}, {80, 90}, {90, 100}}; //{{0, 20}, {0, 20}, {20, 40}, {20, 40}, {40, 60}, {40, 60}, {60, 80}, {60, 80}, {80, 90}};
  for (int iC{1}; iC <= 8; ++iC) {
    projNum[iC - 1] = (TH1F*)hRec->ProjectionY(Form("eff_%d", iC), centBinsP[iC - 1][0] + 1, centBinsP[iC - 1][1], 1, 1);
    TH1F projNumCpy(*projNum[iC - 1]);
    projNumCpy.SetName(Form("projNum_%d", iC));
    auto projDen = (TH1F*)hGen->ProjectionY(Form("projDen_%d", iC), centBinsP[iC - 1][0] + 1, centBinsP[iC - 1][1], 1, 1);
    // projNum[iC - 1]->Rebin(2);
    // projDen->Rebin(2);
    projNum[iC - 1]->Divide(projNum[iC - 1], projDen, 1., 1., "B");
    projNum[iC - 1]->GetYaxis()->SetTitle("Efficiency");
    fileOut->cd();
    //projDen->Write();
    //projNumCpy.Write();
    projNum[iC - 1]->Write(Form("effP_%d", iC));
  }
  fileOut->Close();
}
