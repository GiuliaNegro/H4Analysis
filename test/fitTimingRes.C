#include "TH1F.h"
#include "TH2F.h"
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


void fitTimingRes(){

	TCanvas* mycanvas = new TCanvas( "mycanvas", "", 0, 0, 800, 600 ) ;

	const Int_t n = 6;
	Double_t x[n];
	Double_t ex[n];


/// 11433, 11423, 11416, 11409, 11393, 11357   ---   160 MHz, 18C
	//VFE-MCP1 res
	Double_t y1[n] = {0.0712517, 0.0517905, 0.0437648, 0.0419162, 0.0412569, 0.0405452}; 
	Double_t ey1[n] = {0.0010588, 0.00111595, 0.000703376, 0.000664515, 0.000543421, 0.00206717};
	// //VFE-MCP2 res
	// Double_t y1[n] = {0.0712673, 0.0520744, 0.0438922, 0.0411946, 0.0403581, 0.0395889}; 
	// Double_t ey1[n] = {0.00103976, 0.00107723, 0.000708053, 0.000669626, 0.00055822, 0.0026817};

	//MCP1-MCP2 res (event by event)
	Double_t y2[n] = {0.03225, 0.03165, 0.03135, 0.03135, 0.03165, 0.03255}; 
	Double_t ey2[n] = {0.000824621, 0.000707107, 0.000974679, 0.00116619, 0.00126095, 0.000424264};

	Double_t y[n]; 
	Double_t ey[n];
	for(int i = 0; i < n; i++) {
		y[i] = TMath::Sqrt( pow(y1[i],2) - pow(y2[i],2) );
		ey[i] = TMath::Sqrt( pow(y1[i]/y[i]*ey1[i],2) + pow(y2[i]/y[i]*ey2[i],2) );  
	}

  //amp_fitCB_mean
	Double_t A[n] = {766.785, 1558.51, 3140.44, 4736.72, 6283.76, 7715.79};
	Double_t A_err[n] = {0.432666, 0.640414, 0.807693, 1.09376, .29452, 5.01785};

  //noise from rms samples
	Double_t sig_noise[n] = {4.375, 4.025, 4.375, 4.725, 4.375, 4.025}; 
///


	for(int i = 0; i < n; i++) {
		x[i] = A[i]/sig_noise[i];
		ex[i] = TMath::Sqrt(pow(A_err[i]/A[i], 2) + pow(0/sig_noise[i], 2)) * A[i]/sig_noise[i];
		y[i] = y[i]*1000;
		ey[i] = ey[i]*1000;
		cout << x[i] << "+/-" << ex[i] << ", " << y[i] << "+/-" << ey[i] << endl;
	} 

	TGraphErrors* gr = new TGraphErrors(n,x,y,ex,ey);
	gr->GetXaxis()->SetTitle("A/#sigma_{n}");
	gr->GetYaxis()->SetTitle("#sigma(t_{VFE}-t_{MCP}) (ps)");
	gr->GetYaxis()->SetTitleOffset(1.2);
	gr->SetTitle("Timing Resolution");
	gr->SetMarkerColor(4);
	gr->SetMarkerStyle(21);
	gr->Draw("AP");


	TF1 *myfit = new TF1("myfit","sqrt( (pow([0],2)) + pow(([1]/x) , 2) )"); 
	myfit->SetParName(0,"a");
	myfit->SetParName(1,"b");
	myfit->SetParameter(0, 30);
	myfit->SetParameter(1, 10000);
	myfit->SetLineColor(kTeal+4);
	gr->Fit("myfit","Q");

	cout << "a = " << myfit->GetParameter(0) << " +/- " << myfit->GetParError(0) << endl;
	cout << "b = " << myfit->GetParameter(1) << " +/- " << myfit->GetParError(1) << endl;

    Double_t chi2 = myfit->GetChisquare();
    Int_t ndf = myfit->GetNDF();
    cout<<"chi2="<<chi2<<", ndf="<<ndf<<", chi2/ndf="<<chi2/ndf<<endl;

    std::ostringstream p0;
    p0 << roundf(myfit->GetParameter(0)*100)/100;
    std::string p0_str = p0.str();

    std::ostringstream p0e;
    p0e << roundf(myfit->GetParError(0)*100)/100;
    std::string p0e_str = p0e.str();

    std::ostringstream p1;
    p1 << roundf(myfit->GetParameter(1)*100)/100;
    std::string p1_str = p1.str();

    std::ostringstream p1e;
    p1e << roundf(myfit->GetParError(1)*100)/100;
    std::string p1e_str = p1e.str();

	std::ostringstream chi2Ndf;
	chi2Ndf << roundf(chi2/ndf *100)/100;
	std::string chi2Ndf_str = chi2Ndf.str();

	TString p0text = "C = " + p0_str + " +/- " + p0e_str + " ps";
	TString p1text = "N = " + p1_str + " +/- " + p1e_str + " ps";
	TString chi2text = "#chi_{2}/ndf = " + chi2Ndf_str;


	TLatex latex;
	latex.SetNDC();
	latex.SetTextAngle(0);
	latex.SetTextColor(kBlack);

	latex.SetTextFont(42);
	latex.SetTextSize(0.03);
	latex.DrawLatex(0.5,0.74,p0text);
    latex.DrawLatex(0.5,0.68,p1text);
    latex.DrawLatex(0.5,0.62,chi2text);

	float l = mycanvas->GetLeftMargin();
	float t = mycanvas->GetTopMargin(); 
	float r = mycanvas->GetRightMargin();
	float b = mycanvas->GetBottomMargin();
	float xText = l + 0.045*(1-l-r);
	float yText = 1 - t + 0.01; 

	latex.SetTextFont(61);
	latex.SetTextSize(0.45*t);
	latex.DrawLatex(xText,yText,"CMS");

	latex.SetTextFont(52);
	latex.SetTextSize(0.4*t);
	latex.DrawLatex(xText + 0.09, yText, "Preliminary");

	latex.SetTextAlign(11);
	latex.SetTextFont(42);
	latex.SetTextSize(0.35*t);
	latex.DrawLatex(1-r-0.25,yText-0.06,"160 MHz - 18#circ C");
	// latex.DrawLatex(1-r-0.25,yText-0.06,"160 MHz - 9#circ C");
	// latex.DrawLatex(1-r-0.25,yText-0.06,"120 MHz - 18#circ C");
	// latex.DrawLatex(1-r-0.25,yText-0.06,"120 MHz - 9#circ C");


	TString configuration = "160MHz_18C"; 
	// TString configuration = "160MHz_9C"; 
	// TString configuration = "120MHz_18C"; 
	// TString configuration = "120MHz_9C"; 

	TString outputFile = "~/www/TBJune2018_H4Analysis/hodo/"+configuration+"/timingResolution_C3-MCP1_ampFit_noise_"+configuration+"_minusMCPresFit";
	// TString outputFile = "~/www/TBJune2018_H4Analysis/hodo/"+configuration+"/timingResolution_C3-MCP2_ampFit_noise_"+configuration+"_minusMCPresFit";

	mycanvas->SaveAs(outputFile+".png");
	mycanvas->SaveAs(outputFile+".root");


	TFile *f = new TFile(outputFile+"_TF1.root", "RECREATE");
	gr->Write();
	myfit->Write();
	f->Close();

}

