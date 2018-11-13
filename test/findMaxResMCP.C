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


void findMaxResMCP(int run){ 

  string inputDir = "~/www/TBJune2018_H4Analysis/hodo/"+to_string(run)+"/";
  cout<<inputDir<<endl;

  string nameFile = "resMCPs";
  TFile inputFile((inputDir+nameFile+".root").c_str());
  TCanvas* inputCanv = static_cast<TCanvas*>(inputFile.Get((nameFile).c_str()));
  TH1F* inputHist = static_cast<TH1F*>(inputCanv->FindObject((nameFile+"_h").c_str()));

  Int_t maxBin = inputHist->GetMaximumBin();
  Float_t res = inputHist->GetBinCenter(maxBin);
  Float_t res_err = inputHist->GetBinError(maxBin);
  cout<<"resMCPs="<<res<<" +/- "<<res_err<<endl;

}

