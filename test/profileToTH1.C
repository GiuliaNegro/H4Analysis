#include "TH1D.h"
#include "TH2D.h"
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


void profileToTH1(int run, string nameFile){

  string inputDir = "~/www/TBJune2018_H4Analysis/hodo/"+to_string(run)+"/";
  cout<<inputDir<<endl;

  TFile inputFile((inputDir+nameFile+".root").c_str());
  TCanvas* inputCanv = static_cast<TCanvas*>(inputFile.Get((nameFile).c_str()));
  TProfile* inputProf = static_cast<TProfile*>(inputCanv->FindObject((nameFile+"_hist").c_str()));

/// for WF not oversampled
  TH1D *prof = new TH1D("prof", "", 1000, -25, 125);

  Double_t half_bin  = inputProf->FindFirstBinAbove(inputProf->GetMaximum()/2.);
  int t_shift = int((inputProf->GetBinCenter(half_bin))*prof->GetBinWidth(1)); 
  // cout<<inputProf->GetBin(inputProf->GetMaximumBin())<<", "<<inputProf->GetBinCenter(inputProf->GetMaximumBin())<<endl;
  // cout<<inputProf->GetMaximum()<<", "<<half_bin<<", "<<inputProf->GetBinCenter(half_bin)<<", "<<prof->GetBinWidth(1)<<", "<<t_shift<<endl;

  for(int ibin=1; ibin<1001; ibin++) {
    int jbin = ibin-t_shift;
    if (jbin<1) continue;
    prof->SetBinContent(jbin, inputProf->GetBinContent(ibin));
    prof->SetBinError(jbin, inputProf->GetBinError(ibin));
    // cout<<"bin "<<ibin<<": "<<inputProf->GetBinCenter(ibin)<<", "<<inputProf->GetBinContent(ibin)<<" -> ";
    // cout<<"bin "<<jbin<<": "<<prof->GetBinCenter(jbin)<<", "<<prof->GetBinContent(jbin)<<endl;
  }
///

  TCanvas* outputCanv = new TCanvas("outputCanv");

  // TH1D *prof = inputProf->ProjectionX();   
  prof->SetDirectory(0);
  prof->Draw();

  // for (int i=0; i<270; i++) 
    // cout<<"bin "<<i<<": time="<<(i*0.15)-25.<<", val="<<prof->GetBinContent(i)<<endl;

  string outputDir = "shapes/";
  outputCanv->SaveAs((outputDir+nameFile+"_TH1.root").c_str()); 

}
