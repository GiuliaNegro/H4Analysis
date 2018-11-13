#include "TH1.h"
#include "TH2.h"
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
#include "TSystem.h"
#include "TLorentzVector.h"
#include <vector>
#include <iostream>
#include <string>
#include <sstream>
#include <cstdio>
#include <memory>
#include <fstream>
#include <math.h>

using namespace std;


void drawResultsTiming()
{
	TString inputDir = "~/www/TBJune2018_H4Analysis/hodo/";
	TString nameFile = "timingResolution_C3-MCP2_ampFit_noise_";

	TString config1 = "160MHz_18C";
	TString config2 = "160MHz_9C";
	TString config3 = "120MHz_18C";
	TString config4 = "120MHz_9C";

	TFile inputFile1(inputDir+config1+"/"+nameFile+config1+"_minusMCPresFit_TF1.root");
	TFile inputFile2(inputDir+config2+"/"+nameFile+config2+"_minusMCPresFit_TF1.root");
	TFile inputFile3(inputDir+config3+"/"+nameFile+config3+"_minusMCPresFit_TF1.root");
	TFile inputFile4(inputDir+config4+"/"+nameFile+config4+"_minusMCPresFit_TF1.root");

	TGraphErrors* h_tmpl1=(TGraphErrors*)inputFile1.Get("Graph");
	h_tmpl1->SetName("graph1");
	TF1* fit1=(TF1*)inputFile1.Get("myfit");
	fit1->SetName("fit1");

	TGraphErrors* h_tmpl2=(TGraphErrors*)inputFile2.Get("Graph");
	h_tmpl2->SetName("graph2");
	TF1* fit2=(TF1*)inputFile2.Get("myfit");
	fit2->SetName("fit2");

	TGraphErrors* h_tmpl3=(TGraphErrors*)inputFile3.Get("Graph");
	h_tmpl3->SetName("graph3");
	TF1* fit3=(TF1*)inputFile3.Get("myfit");
	fit3->SetName("fit3");

	TGraphErrors* h_tmpl4=(TGraphErrors*)inputFile4.Get("Graph");
	h_tmpl4->SetName("graph4");
	TF1* fit4=(TF1*)inputFile4.Get("myfit");
	fit4->SetName("fit4");

	gStyle->SetOptStat(0);


	TCanvas *canv_output = new TCanvas("canv_output","canv_output", 1100, 900);
	gPad->SetTickx();
	gPad->SetTicky();

	h_tmpl1->GetXaxis()->SetLimits(0., 2400);
	h_tmpl1->GetYaxis()->SetRangeUser(20., 116.);
	h_tmpl1->SetMarkerColor(kTeal+4);
	h_tmpl1->SetMarkerStyle(21);
	h_tmpl1->Draw("AP");
	fit1->SetLineColor(kTeal+4); 
	h_tmpl1->Fit("fit1","Q");

	h_tmpl2->SetMarkerColor(kBlue);
	h_tmpl2->SetMarkerStyle(25);
	h_tmpl2->Draw("PSAME");
	fit2->SetLineColor(kBlue);
	fit2->SetLineStyle(2);
	h_tmpl2->Fit("fit2","Q");

	h_tmpl3->SetMarkerColor(kTeal+4);
	h_tmpl3->SetMarkerStyle(26);
	h_tmpl3->SetMarkerSize(1.5);
	h_tmpl3->Draw("PSAME");
	fit3->SetLineColor(kTeal+4);
	fit3->SetLineStyle(2);
	h_tmpl3->Fit("fit3","Q");

	h_tmpl4->SetMarkerColor(kBlue);
	h_tmpl4->SetMarkerStyle(22);
	h_tmpl4->SetMarkerSize(1.5);
	h_tmpl4->Draw("PSAME");
	fit4->SetLineColor(kBlue);
	h_tmpl4->Fit("fit4","Q");

	TLegend* legend = new TLegend(0.6, 0.65, 0.8, 0.85);
	legend->AddEntry(h_tmpl1, config1, "P");
	legend->AddEntry(h_tmpl2, config2, "P");	
	legend->AddEntry(h_tmpl3, config3, "P");	
	legend->AddEntry(h_tmpl4, config4, "P");
	legend->SetFillColor( kWhite ) ; 
	legend->SetLineColor( kWhite ) ; 
	legend->Draw("SAME");


	TLatex latex;
	latex.SetNDC();
	latex.SetTextAngle(0);
	latex.SetTextColor(kBlack);

	float l = canv_output->GetLeftMargin();
	float t = canv_output->GetTopMargin(); 
	float r = canv_output->GetRightMargin();
	float b = canv_output->GetBottomMargin();
	float xText = l + 0.045*(1-l-r);
	float yText = 1 - t + 0.01; 

	latex.SetTextFont(61);
	latex.SetTextSize(0.45*t);
	latex.DrawLatex(xText,yText,"CMS");

	latex.SetTextFont(52);
	latex.SetTextSize(0.4*t);
	latex.DrawLatex(xText + 0.09, yText, "Preliminary");


	canv_output->cd();
	string outputDir = "~/www/TBJune2018_H4Analysis/hodo/";
	string outputFile = "summaryTimingResolution"; 

	canv_output->SaveAs( (outputDir+outputFile+".png").c_str() );
	canv_output->SaveAs( (outputDir+outputFile+".pdf").c_str() );
	canv_output->SaveAs( (outputDir+outputFile+".root").c_str() );

}
