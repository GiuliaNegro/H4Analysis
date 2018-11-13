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


void fitAmplPeak(int run){ 

  string inputDir = "~/www/TBJune2018_H4Analysis/hodo/"+to_string(run)+"/";
  cout<<inputDir<<endl;


  string nameNoiseFile = "noise_spectrum_C3";
  TFile inputNoiseFile((inputDir+nameNoiseFile+".root").c_str());

  TCanvas* inputNoiseCanv = static_cast<TCanvas*>(inputNoiseFile.Get((nameNoiseFile).c_str()));
  TH1F* inputNoiseHist = static_cast<TH1F*>(inputNoiseCanv->FindObject((nameNoiseFile+"_h").c_str()));

  Int_t maxBin = inputNoiseHist->GetMaximumBin();
  Float_t noise = inputNoiseHist->GetBinCenter(maxBin);
  cout<<"noise="<<noise<<endl;


  string nameAmpMCP1File = "amp_maximum_MCP1";
  TFile inputAmpMCP1File((inputDir+nameAmpMCP1File+".root").c_str());
  TCanvas* inputAmpMCP1Canv = static_cast<TCanvas*>(inputAmpMCP1File.Get((nameAmpMCP1File).c_str()));
  TH1F* inputAmpMCP1Hist = static_cast<TH1F*>(inputAmpMCP1Canv->FindObject((nameAmpMCP1File+"_h").c_str()));
  Int_t maxAmpMCP1Bin = inputAmpMCP1Hist->GetMaximumBin();
  Float_t ampMCP1 = inputAmpMCP1Hist->GetBinCenter(maxAmpMCP1Bin);
  if (ampMCP1<100) {
    inputAmpMCP1Hist->GetXaxis()->SetRangeUser(100, 8000);
    maxAmpMCP1Bin = inputAmpMCP1Hist->GetMaximumBin();
    ampMCP1 = inputAmpMCP1Hist->GetBinCenter(maxAmpMCP1Bin);
  }
  cout<<"ampMCP1="<<ampMCP1<<endl;


  string nameAmpMCP2File = "amp_maximum_MCP2";
  TFile inputAmpMCP2File((inputDir+nameAmpMCP2File+".root").c_str());
  TCanvas* inputAmpMCP2Canv = static_cast<TCanvas*>(inputAmpMCP2File.Get((nameAmpMCP2File).c_str()));
  TH1F* inputAmpMCP2Hist = static_cast<TH1F*>(inputAmpMCP2Canv->FindObject((nameAmpMCP2File+"_h").c_str()));
  Int_t maxAmpMCP2Bin = inputAmpMCP2Hist->GetMaximumBin();
  Float_t ampMCP2 = inputAmpMCP2Hist->GetBinCenter(maxAmpMCP2Bin);
  if (ampMCP2<100) {
    inputAmpMCP2Hist->GetXaxis()->SetRangeUser(100, 8000);
    maxAmpMCP2Bin = inputAmpMCP2Hist->GetMaximumBin();
    ampMCP2 = inputAmpMCP2Hist->GetBinCenter(maxAmpMCP2Bin);
  }
  cout<<"ampMCP2="<<ampMCP2<<endl;



  string nameFile = "fit_ampl_C3";
  TFile inputFile((inputDir+nameFile+".root").c_str());
  TCanvas* inputCanv = static_cast<TCanvas*>(inputFile.Get((nameFile).c_str()));
  TH1F* inputHist = static_cast<TH1F*>(inputCanv->FindObject((nameFile+"_h").c_str()));
  inputHist->SetName((nameFile).c_str());

  TCanvas* outputCanv = new TCanvas("outputCanv", "outputCanv", 1100, 900);
  outputCanv->cd();

  Int_t qmaxBin = inputHist->GetMaximumBin();
  Float_t qmax_VFE = inputHist->GetBinCenter(qmaxBin);
  cout<<"amp_max_VFE="<<qmax_VFE<<endl;

  Float_t qmax_VFE_first = qmax_VFE*0.3;
  Float_t qmax_VFE_last = qmax_VFE*1.5;

  inputHist->GetXaxis()->SetRangeUser(qmax_VFE_first, qmax_VFE_last);


  TF1 *fitFunc = new TF1("fitFunc","crystalball"); 
  fitFunc->SetParameters(0, qmax_VFE, 20, 1, 3);
  // fitFunc->SetParameters(0, qmax_VFE, 20, 0.5, 5); 
  fitFunc->SetLineColor(kBlue);

  inputHist->Fit("fitFunc", "Q");
  gStyle->SetOptFit();

  inputHist->GetXaxis()->SetRangeUser(qmax_VFE*0.6, qmax_VFE*1.5);
  // inputHist->GetXaxis()->SetRange((int)(inputHist->GetMaximumBin()*0.50), (int)(inputHist->GetMaximumBin()*1.50));

  inputHist->GetXaxis()->SetTitle("fit_amp (ADC counts)");
  inputHist->GetXaxis()->SetTitleOffset(1.1);
  inputHist->GetXaxis()->SetTitleSize(0.04);
  inputHist->GetXaxis()->SetLabelSize(0.04);

  inputHist->GetYaxis()->SetTitle("Events");
  inputHist->GetYaxis()->SetTitleOffset(1.4);
  inputHist->GetYaxis()->SetTitleSize(0.04);
  inputHist->GetYaxis()->SetLabelSize(0.04);

  inputHist->Draw();
  fitFunc->Draw("SAME");

  Float_t mean = fitFunc->GetParameter(1);
  Float_t meanErr = fitFunc->GetParError(1);
  Float_t sigma = fitFunc->GetParameter(2);
  Float_t sigmaErr = fitFunc->GetParError(2);
  cout<<"amp_max_VFE: mean="<<mean<<", mean_err="<<meanErr<<", sigma="<<sigma<<", sigma_err="<<sigmaErr<<endl;


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
  // latex.DrawLatex(x,y-0.07,"160 MHz - 18#circ C");
  // latex.DrawLatex(x,y-0.07,"160 MHz - 9#circ C");
  latex.DrawLatex(x,y-0.07,"120 MHz - 18#circ C");
  // latex.DrawLatex(x,y-0.07,"120 MHz - 9#circ C");

  outputCanv->SaveAs((inputDir+nameFile+"_fitCB.png").c_str()); 
  outputCanv->SaveAs((inputDir+nameFile+"_fitCB.root").c_str()); 
  // inputHist->GetXaxis()->SetRangeUser(3500., 5700.);
  // outputCanv->SaveAs((inputDir+nameFile+"_fitCB_zoom.png").c_str()); 




  // string nameFileAmp = "fit_ampl_C3";
  // TFile inputFileAmp((inputDir+nameFileAmp+".root").c_str());

  // TCanvas* inputCanvAmp = static_cast<TCanvas*>(inputFileAmp.Get((nameFileAmp).c_str()));
  // TH1F* inputHistAmp = static_cast<TH1F*>(inputCanvAmp->FindObject((nameFileAmp+"_h").c_str()));
  // inputHistAmp->SetName((nameFileAmp).c_str());

  // TCanvas* canv = new TCanvas("canv", "canv", 700, 600);
  // canv->cd();
  // gStyle->SetOptStat(0);

  // inputHistAmp->GetXaxis()->SetTitle("fit_amp (ADC counts)");
  // inputHistAmp->GetXaxis()->SetTitleOffset(1.1);
  // inputHistAmp->GetXaxis()->SetTitleSize(0.04);
  // inputHistAmp->GetXaxis()->SetLabelSize(0.04);
  // inputHistAmp->GetYaxis()->SetTitle("Events");
  // inputHistAmp->GetYaxis()->SetTitleOffset(1.4);
  // inputHistAmp->GetYaxis()->SetTitleSize(0.04);
  // inputHistAmp->GetYaxis()->SetLabelSize(0.04);
  // inputHistAmp->Draw();

  // l = canv->GetLeftMargin();
  // t = canv->GetTopMargin(); 
  // r = canv->GetRightMargin();
  // b = canv->GetBottomMargin();
  // x = l + 0.045*(1-l-r);
  // y = 1 - t + 0.01; 

  // latex.SetTextFont(61);
  // latex.SetTextSize(0.45*t);
  // latex.DrawLatex(x,y,"CMS");

  // latex.SetTextFont(52);
  // latex.SetTextSize(0.4*t);
  // latex.DrawLatex(x + 0.09, y, "Preliminary");

  // latex.SetTextAlign(11);
  // latex.SetTextFont(42);
  // latex.SetTextSize(0.35*t);
  // latex.DrawLatex(x,y-0.06,"160 MHz - 18#circ C");

  // canv->SaveAs((inputDir+nameFileAmp+".png").c_str()); 



  // TCanvas* canv1 = new TCanvas("canv1", "canv1", 700, 600);
  // canv1->cd();
  // gStyle->SetOptStat(0);

  // inputNoiseHist->SetName((nameNoiseFile).c_str()); 
  // inputNoiseHist->GetXaxis()->SetTitle("noise (ADC counts)");
  // inputNoiseHist->GetXaxis()->SetTitleOffset(1.1);
  // inputNoiseHist->GetXaxis()->SetTitleSize(0.04);
  // inputNoiseHist->GetXaxis()->SetLabelSize(0.04);
  // inputNoiseHist->GetYaxis()->SetTitle("Events");
  // inputNoiseHist->GetYaxis()->SetTitleOffset(1.4);
  // inputNoiseHist->GetYaxis()->SetTitleSize(0.04);
  // inputNoiseHist->GetYaxis()->SetLabelSize(0.04);
  // inputNoiseHist->Draw();

  // l = canv1->GetLeftMargin();
  // t = canv1->GetTopMargin(); 
  // r = canv1->GetRightMargin();
  // b = canv1->GetBottomMargin();
  // x = l + 0.045*(1-l-r);
  // y = 1 - t + 0.01; 

  // latex.SetTextFont(61);
  // latex.SetTextSize(0.45*t);
  // latex.DrawLatex(x,y,"CMS");

  // latex.SetTextFont(52);
  // latex.SetTextSize(0.4*t);
  // latex.DrawLatex(x + 0.09, y, "Preliminary");

  // latex.SetTextAlign(11);
  // latex.SetTextFont(42);
  // latex.SetTextSize(0.35*t);
  // latex.DrawLatex(1-r-0.25,y-0.06,"160 MHz - 18#circ C");

  // inputNoiseHist->GetXaxis()->SetRangeUser(0., 20.);
  // canv1->SaveAs((inputDir+nameNoiseFile+".png").c_str()); 



  // TCanvas* canv2 = new TCanvas("canv2", "canv2", 700, 600);
  // canv2->cd();
  // gStyle->SetOptStat(0);

  // inputAmpMCP1Hist->SetName((nameAmpMCP1File).c_str()); 
  // inputAmpMCP1Hist->GetXaxis()->SetRangeUser(0., 3000.);
  // inputAmpMCP1Hist->GetXaxis()->SetTitle("amp_maximum (ADC counts)");
  // inputAmpMCP1Hist->GetXaxis()->SetTitleOffset(1.1);
  // inputAmpMCP1Hist->GetXaxis()->SetTitleSize(0.04);
  // inputAmpMCP1Hist->GetXaxis()->SetLabelSize(0.04);
  // inputAmpMCP1Hist->GetYaxis()->SetTitle("Events");
  // inputAmpMCP1Hist->GetYaxis()->SetTitleOffset(1.4);
  // inputAmpMCP1Hist->GetYaxis()->SetTitleSize(0.04);
  // inputAmpMCP1Hist->GetYaxis()->SetLabelSize(0.04);
  // inputAmpMCP1Hist->Draw();

  // l = canv2->GetLeftMargin();
  // t = canv2->GetTopMargin(); 
  // r = canv2->GetRightMargin();
  // b = canv2->GetBottomMargin();
  // x = l + 0.045*(1-l-r);
  // y = 1 - t + 0.01; 

  // latex.SetTextFont(61);
  // latex.SetTextSize(0.45*t);
  // latex.DrawLatex(x,y,"CMS");

  // latex.SetTextFont(52);
  // latex.SetTextSize(0.4*t);
  // latex.DrawLatex(x + 0.09, y, "Preliminary");

  // latex.SetTextAlign(11);
  // latex.SetTextFont(42);
  // latex.SetTextSize(0.35*t);
  // latex.DrawLatex(1-r-0.25,y-0.06,"160 MHz - 18#circ C");

  // canv2->SaveAs((inputDir+nameAmpMCP1File+".png").c_str()); 



  // TCanvas* canv3 = new TCanvas("canv3", "canv3", 700, 600);
  // canv3->cd();
  // gStyle->SetOptStat(0);

  // inputAmpMCP2Hist->SetName((nameAmpMCP2File).c_str()); 
  // inputAmpMCP2Hist->GetXaxis()->SetRangeUser(0., 3000.);
  // inputAmpMCP2Hist->GetXaxis()->SetTitle("amp_maximum (ADC counts)");
  // inputAmpMCP2Hist->GetXaxis()->SetTitleOffset(1.1);
  // inputAmpMCP2Hist->GetXaxis()->SetTitleSize(0.04);
  // inputAmpMCP2Hist->GetXaxis()->SetLabelSize(0.04);
  // inputAmpMCP2Hist->GetYaxis()->SetTitle("Events");
  // inputAmpMCP2Hist->GetYaxis()->SetTitleOffset(1.4);
  // inputAmpMCP2Hist->GetYaxis()->SetTitleSize(0.04);
  // inputAmpMCP2Hist->GetYaxis()->SetLabelSize(0.04);
  // inputAmpMCP2Hist->Draw();

  // l = canv3->GetLeftMargin();
  // t = canv3->GetTopMargin(); 
  // r = canv3->GetRightMargin();
  // b = canv3->GetBottomMargin();
  // x = l + 0.045*(1-l-r);
  // y = 1 - t + 0.01; 

  // latex.SetTextFont(61);
  // latex.SetTextSize(0.45*t);
  // latex.DrawLatex(x,y,"CMS");

  // latex.SetTextFont(52);
  // latex.SetTextSize(0.4*t);
  // latex.DrawLatex(x + 0.09, y, "Preliminary");

  // latex.SetTextAlign(11);
  // latex.SetTextFont(42);
  // latex.SetTextSize(0.35*t);
  // latex.DrawLatex(1-r-0.25,y-0.06,"160 MHz - 18#circ C");

  // canv3->SaveAs((inputDir+nameAmpMCP2File+".png").c_str()); 


}