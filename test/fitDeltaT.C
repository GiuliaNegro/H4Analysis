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


void fitDeltaT(int run){ 

  string inputDir = "~/www/TBJune2018_H4Analysis/hodo/"+to_string(run)+"/";
  cout<<inputDir<<endl;
  
  string nameFile = "dt_fit_VFE_MCP1_mod";
  // string nameFile = "dt_fit_VFE_MCP2_mod";
  // string nameFile = "dt0_MCP1_MCP2";

  TFile inputFile((inputDir+nameFile+".root").c_str());
  TCanvas* inputCanv = static_cast<TCanvas*>(inputFile.Get((nameFile).c_str()));
  TH1F* inputHist = static_cast<TH1F*>(inputCanv->FindObject((nameFile+"_h").c_str()));
  inputHist->SetName((nameFile).c_str());

  TCanvas* outputCanv = new TCanvas("outputCanv", "outputCanv", 1100, 900);
  outputCanv->cd();

  Int_t dtmaxBin = inputHist->GetMaximumBin();
  cout<<"dtmaxBin="<<dtmaxBin<<endl;
  Float_t dtmax_VFE_MCP = inputHist->GetBinCenter(dtmaxBin);
  cout<<"dtmax_VFE_MCP="<<dtmax_VFE_MCP<<endl;

  inputHist->GetXaxis()->SetRangeUser(3.5,5.5);
  inputHist->GetXaxis()->SetTitle("fit_time-t0_MCP1+t0_CLK (ns)");
  inputHist->GetXaxis()->SetTitleOffset(1.1);
  inputHist->GetXaxis()->SetTitleSize(0.04);
  inputHist->GetXaxis()->SetLabelSize(0.04);

  inputHist->GetYaxis()->SetTitle("Events");
  inputHist->GetYaxis()->SetTitleOffset(1.4);
  inputHist->GetYaxis()->SetTitleSize(0.04);
  inputHist->GetYaxis()->SetLabelSize(0.04);

  inputHist->Fit("gaus");
  gStyle->SetOptFit();
  TF1 *fit = inputHist->GetFunction("gaus");
  Double_t chi2 = fit->GetChisquare();
  Int_t ndf = fit->GetNDF();
  Double_t p0 = fit->GetParameter(0);
  Double_t p1 = fit->GetParameter(1);
  Double_t p2 = fit->GetParameter(2);
  Double_t e0 = fit->GetParError(0);
  Double_t e1 = fit->GetParError(1);
  Double_t e2 = fit->GetParError(2);
  cout<<"chi2="<<chi2<<", ndf="<<ndf<<", chi2/ndf="<<chi2/ndf<<", p0="<<p0<<" +/- "<<e0<<", p1="<<p1<<" +/- "<<e1<<", p2="<<p2<<" +/- "<<e2<<endl;
  
  inputHist->Draw();

  TLatex latex;
  latex.SetNDC();
  latex.SetTextAngle(0);
  latex.SetTextColor(kBlack);

  float l = outputCanv->GetLeftMargin();
  float t = outputCanv->GetTopMargin(); 
  float r = outputCanv->GetRightMargin();
  float b = outputCanv->GetBottomMargin();
  float x = l + 0.045*(1-l-r);
  float y = 1 - t + 0.01; 

  latex.SetTextAlign(11); //31 on the right

  latex.SetTextFont(61);
  latex.SetTextSize(0.45*t);
  latex.DrawLatex(x,y,"CMS");

  latex.SetTextFont(52);
  latex.SetTextSize(0.4*t);
  latex.DrawLatex(x + 0.08, y, "Preliminary");

  latex.SetTextFont(42);
  latex.SetTextSize(0.35*t);
  latex.DrawLatex(x,y-0.07,"160 MHz - 18#circ C");
  // latex.DrawLatex(x,y-0.07,"160 MHz - 9#circ C");
  // latex.DrawLatex(x,y-0.07,"120 MHz - 18#circ C");
  // latex.DrawLatex(x,y-0.07,"120 MHz - 9#circ C");

  inputHist->GetXaxis()->SetRangeUser(3.9,4.7);

  outputCanv->SaveAs((inputDir+nameFile+"_fit.png").c_str()); 
  outputCanv->SaveAs((inputDir+nameFile+"_fit.root").c_str()); 

}