#include "TH1D.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TProfile.h"
#include "TGraphErrors.h"
#include "TF1.h"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TCanvas.h"
#include "THStack.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TStyle.h"
#include "TMath.h"
#include <TSystem.h>
#include "TLorentzVector.h"
#include <vector>
#include <iostream>
#include <string>
#include <sstream>
#include <cstdio>
#include <memory>

using namespace std;


void makeTH2(string inputDir) {

  string nameFile = "dt0_MCP1_MCP2_vs_amp_eff";
  TFile inputFile((inputDir+nameFile+".root").c_str());
  TCanvas* inputCanv = static_cast<TCanvas*>(inputFile.Get((nameFile).c_str()));
  TH2F* inputHist = static_cast<TH2F*>(inputCanv->FindObject((nameFile+"_h").c_str()));

  TFile *file = new TFile((inputDir+"dtMCPs_vs_ampEff.root").c_str(), "recreate"); 
  // inputHist->RebinX(3);
  // TFile *file = new TFile((inputDir+"dtMCPs_vs_ampEff_rebin.root").c_str(), "recreate"); 

  inputHist->Write();  
  file->Close();

}



void fitMCPresolution(){ 
  string inputDir = "~/www/TBJune2018_H4Analysis/hodo/160MHz_18C/";
  // string inputDir = "~/www/TBJune2018_H4Analysis/hodo/160MHz_9C/";
  // string inputDir = "~/www/TBJune2018_H4Analysis/hodo/120MHz_18C/";
  // string inputDir = "~/www/TBJune2018_H4Analysis/hodo/120MHz_9C/";
  cout<<inputDir<<endl;

  makeTH2(inputDir);

  string nameFile = "dtMCPs_vs_ampEff";
  // string nameFile = "dtMCPs_vs_ampEff_rebin";
  TFile inputFile((inputDir+nameFile+".root").c_str());
  TH2F* hist = (TH2F*)inputFile.Get("dt0_MCP1_MCP2_vs_amp_eff_h");

  TCanvas* outputCanv = new TCanvas("outputCanv", "outputCanv", 1100, 900);
  outputCanv->cd();

  hist->FitSlicesY();

  TH2F *hist_2 = (TH2F*)inputFile.Get("dt0_MCP1_MCP2_vs_amp_eff_h_2");
  hist_2->SetStats(0);
  hist_2->SetTitle("");
  hist_2->GetYaxis()->SetTitle("#sigma(t_{MCP1}-t_{MCP2}) (ns)");
  hist_2->GetYaxis()->SetTitleOffset(1.2);
  hist_2->GetYaxis()->SetRangeUser(-0.15, 1.25);
  hist_2->GetXaxis()->SetRangeUser(0., 1600.);
  hist_2->GetXaxis()->SetTitle("amp_eff_MCPs (ADC counts)");
  hist_2->Draw();

  outputCanv->SaveAs((inputDir+"res_"+nameFile+".png").c_str()); 
  outputCanv->SaveAs((inputDir+"res_"+nameFile+".root").c_str()); 


  hist_2->GetXaxis()->SetRangeUser(100., 1600.);
  TF1 *myfit = new TF1("myfit","sqrt( (2*pow([0],2)) + pow(([1]/x) , 2) )");
  myfit->SetParName(0,"a");
  myfit->SetParName(1,"b");
  myfit->SetLineColor(kTeal+4);
  hist_2->Fit("myfit","Q");

  cout << "a = " << myfit->GetParameter(0) << " +/- " << myfit->GetParError(0) << endl;
  cout << "b = " << myfit->GetParameter(1) << " +/- " << myfit->GetParError(1) << endl;

  Double_t chi2 = myfit->GetChisquare();
  Int_t ndf = myfit->GetNDF();
  cout<<"chi2="<<chi2<<", ndf="<<ndf<<", chi2/ndf="<<chi2/ndf<<endl;

  std::ostringstream p0;
  p0 << roundf(myfit->GetParameter(0)*1000 *100)/100;
  std::string p0_str = p0.str();

  std::ostringstream p0e;
  p0e << roundf(myfit->GetParError(0)*1000 *100)/100;
  std::string p0e_str = p0e.str();

  std::ostringstream p1;
  p1 << roundf(myfit->GetParameter(1)*1000 *100)/100;
  std::string p1_str = p1.str();

  std::ostringstream p1e;
  p1e << roundf(myfit->GetParError(1)*1000 *100)/100;
  std::string p1e_str = p1e.str();

  std::ostringstream chi2Ndf;
  chi2Ndf << roundf(chi2/ndf *100)/100;
  std::string chi2Ndf_str = chi2Ndf.str();

  TString p0text = "C = " + p0_str + " #pm " + p0e_str + " ps";
  TString p1text = "N = " + p1_str + " #pm " + p1e_str + " ps";
  TString chi2text = "#chi_{2}/ndf = " + chi2Ndf_str;

  TLatex latex;
  latex.SetNDC();
  latex.SetTextAngle(0);
  latex.SetTextColor(kBlack);

  latex.SetTextFont(42);
  latex.SetTextSize(0.03);
  latex.DrawLatex(0.5,0.8,p0text);
  latex.DrawLatex(0.5,0.74,p1text);
  latex.DrawLatex(0.5,0.68,chi2text);

  float l = outputCanv->GetLeftMargin();
  float t = outputCanv->GetTopMargin(); 
  float r = outputCanv->GetRightMargin();
  float b = outputCanv->GetBottomMargin();
  float x = l + 0.045*(1-l-r);
  float y = 1 - t + 0.01; 

  latex.SetTextFont(61);
  latex.SetTextSize(0.45*t);
  latex.DrawLatex(x,y,"CMS");

  latex.SetTextFont(52);
  latex.SetTextSize(0.4*t);
  latex.DrawLatex(x + 0.09, y, "Preliminary");

  hist_2->GetXaxis()->SetRangeUser(0., 1600.);
  hist_2->GetYaxis()->SetRangeUser(-0.1, 0.3);

  outputCanv->SaveAs((inputDir+"res_"+nameFile+"_fit.png").c_str()); 
  outputCanv->SaveAs((inputDir+"res_"+nameFile+"_fit.root").c_str()); 


  hist_2->GetXaxis()->SetRangeUser(100., 1000.);
  hist_2->GetYaxis()->SetRangeUser(0., 0.12);
  hist_2->Draw("same");

  outputCanv->SaveAs((inputDir+"res_"+nameFile+"_zoom.png").c_str()); 
  outputCanv->SaveAs((inputDir+"res_"+nameFile+"_zoom.root").c_str()); 

}