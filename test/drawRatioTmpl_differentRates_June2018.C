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


void drawRatioTmpl_Rates()
{
	//160 MHz
	TString name = "shape_tmpl_11409_150GeV_C3";
	TString name1 = "shape_tmpl_11433_25GeV_C3";
	TString name2 = "shape_tmpl_11423_50GeV_C3";
	TString name3 = "shape_tmpl_11416_100GeV_C3";
	TString name4 = "shape_tmpl_11393_200GeV_C3";
	TString name5 = "shape_tmpl_11357_250GeV_C3";

	//120 MHz
	TString name_120 = "shape_tmpl_11478_150GeV_C3";
	TString name1_120 = "shape_tmpl_11456_25GeV_C3";
	TString name2_120 = "shape_tmpl_11468_50GeV_C3";
	TString name3_120 = "shape_tmpl_11472_100GeV_C3";
	TString name4_120 = "shape_tmpl_11485_200GeV_C3";
	TString name5_120 = "shape_tmpl_11494_250GeV_C3";


	TFile inputFile("shapes/"+name+"_TH1.root");
	TFile inputFile1("shapes/"+name1+"_TH1.root");
	TFile inputFile2("shapes/"+name2+"_TH1.root");
	TFile inputFile3("shapes/"+name3+"_TH1.root");
	TFile inputFile4("shapes/"+name4+"_TH1.root");
	TFile inputFile5("shapes/"+name5+"_TH1.root");

	TFile inputFile_120("shapes/"+name_120+"_TH1.root");
	TFile inputFile1_120("shapes/"+name1_120+"_TH1.root");
	TFile inputFile2_120("shapes/"+name2_120+"_TH1.root");
	TFile inputFile3_120("shapes/"+name3_120+"_TH1.root");
	TFile inputFile4_120("shapes/"+name4_120+"_TH1.root");
	TFile inputFile5_120("shapes/"+name5_120+"_TH1.root");


	TCanvas* canv = (TCanvas*)inputFile.Get("outputCanv");
	TCanvas* canv1 = (TCanvas*)inputFile1.Get("outputCanv");
	TCanvas* canv2 = (TCanvas*)inputFile2.Get("outputCanv");
	TCanvas* canv3 = (TCanvas*)inputFile3.Get("outputCanv");
	TCanvas* canv4 = (TCanvas*)inputFile4.Get("outputCanv");
	TCanvas* canv5 = (TCanvas*)inputFile5.Get("outputCanv");

	TCanvas* canv_120 = (TCanvas*)inputFile_120.Get("outputCanv");
	TCanvas* canv1_120 = (TCanvas*)inputFile1_120.Get("outputCanv");
	TCanvas* canv2_120 = (TCanvas*)inputFile2_120.Get("outputCanv");
	TCanvas* canv3_120 = (TCanvas*)inputFile3_120.Get("outputCanv");
	TCanvas* canv4_120 = (TCanvas*)inputFile4_120.Get("outputCanv");
	TCanvas* canv5_120 = (TCanvas*)inputFile5_120.Get("outputCanv");


	TH1* h_tmpl=(TH1*)canv->FindObject(name+"_hist_px");
	TH1* h_tmpl1=(TH1*)canv1->FindObject(name1+"_hist_px");
	TH1* h_tmpl2=(TH1*)canv2->FindObject(name2+"_hist_px");
	TH1* h_tmpl3=(TH1*)canv3->FindObject(name3+"_hist_px");
	TH1* h_tmpl4=(TH1*)canv4->FindObject(name4+"_hist_px");
	TH1* h_tmpl5=(TH1*)canv5->FindObject(name5+"_hist_px");

	TH1* h_tmpl_120=(TH1*)canv_120->FindObject(name_120+"_hist_px");
	TH1* h_tmpl1_120=(TH1*)canv1_120->FindObject(name1_120+"_hist_px");
	TH1* h_tmpl2_120=(TH1*)canv2_120->FindObject(name2_120+"_hist_px");
	TH1* h_tmpl3_120=(TH1*)canv3_120->FindObject(name3_120+"_hist_px");
	TH1* h_tmpl4_120=(TH1*)canv4_120->FindObject(name4_120+"_hist_px");
	TH1* h_tmpl5_120=(TH1*)canv5_120->FindObject(name5_120+"_hist_px");


	h_tmpl->SetTitle("");
	h_tmpl1->SetTitle("");
	h_tmpl2->SetTitle("");
	h_tmpl3->SetTitle("");
	h_tmpl4->SetTitle("");
	h_tmpl5->SetTitle("");

	h_tmpl_120->SetTitle("");
	h_tmpl1_120->SetTitle("");
	h_tmpl2_120->SetTitle("");
	h_tmpl3_120->SetTitle("");
	h_tmpl4_120->SetTitle("");
	h_tmpl5_120->SetTitle("");


	gStyle->SetOptStat(0);

	TCanvas *canv_output = new TCanvas("canv_output","canv_output", 900, 900);

	TLegend* legend = new TLegend(0.6, 0.55, 0.8, 0.85);
	legend->AddEntry(h_tmpl1, "25 GeV, 160 MHz");
	legend->AddEntry(h_tmpl2, "50 GeV, 160 MHz");	
	legend->AddEntry(h_tmpl3, "100 GeV, 160 MHz");	
	legend->AddEntry(h_tmpl, "150 GeV, 160 MHz");
	legend->AddEntry(h_tmpl4, "200 GeV, 160 MHz");
	legend->AddEntry(h_tmpl5, "250 GeV, 160 MHz");

	legend->AddEntry(h_tmpl1_120, "25 GeV, 120 MHz");
	legend->AddEntry(h_tmpl2_120, "50 GeV, 120 MHz");	
	legend->AddEntry(h_tmpl3_120, "100 GeV, 120 MHz");	
	legend->AddEntry(h_tmpl_120, "150 GeV, 120 MHz");
	legend->AddEntry(h_tmpl4_120, "200 GeV, 120 MHz");
	legend->AddEntry(h_tmpl5_120, "250 GeV, 120 MHz");

	legend->SetFillColor( kWhite ) ; 
	legend->SetLineColor( kWhite ) ; 


	TPad *cUp  = new TPad("p1","p1",0,0.26,1,1,0); cUp->Draw();
	TPad *cDown = new TPad("p2","p2",0,0,1,0.24,0); cDown->Draw();
	cUp->SetBottomMargin(0); 
	cDown->SetTopMargin(0); 
	cDown->SetBottomMargin(0.5); 
    
    cUp->cd();
	gPad->SetTickx();
	gPad->SetTicky();

	h_tmpl->GetYaxis()->SetRangeUser(-0.1, 1.05);
	h_tmpl->GetXaxis()->SetRangeUser(-15, 25);
	h_tmpl->SetLineColor(kOrange);
	h_tmpl->SetLineWidth(2);
	h_tmpl->Draw("histl");
	h_tmpl->GetYaxis()->SetTitle("Normalized amplitude (A. U.)");
	h_tmpl->GetYaxis()->SetTitleOffset(0.8); 
	h_tmpl->GetXaxis()->SetLabelSize(0);	

	h_tmpl1->GetYaxis()->SetRangeUser(-0.1, 1.05);
	h_tmpl1->GetXaxis()->SetRangeUser(-15, 25);
	h_tmpl1->SetLineColor(kBlack);
	h_tmpl1->SetLineWidth(2);
	h_tmpl1->Draw("histlSAME");

	h_tmpl2->GetYaxis()->SetRangeUser(-0.1, 1.05);
	h_tmpl2->GetXaxis()->SetRangeUser(-15, 25);
	h_tmpl2->SetLineColor(kBlue);
	h_tmpl2->SetLineWidth(2);
	h_tmpl2->Draw("histlSAME");

	h_tmpl3->GetYaxis()->SetRangeUser(-0.1, 1.05);
	h_tmpl3->GetXaxis()->SetRangeUser(-15, 25);
	h_tmpl3->SetLineColor(kGreen+1);
	h_tmpl3->SetLineWidth(2);
	h_tmpl3->Draw("histlSAME");

	h_tmpl4->GetYaxis()->SetRangeUser(-0.1, 1.05);
	h_tmpl4->GetXaxis()->SetRangeUser(-15, 25);
	h_tmpl4->SetLineColor(kRed);
	h_tmpl4->SetLineWidth(2);
	h_tmpl4->Draw("histlSAME");

	h_tmpl5->GetYaxis()->SetRangeUser(-0.1, 1.05);
	h_tmpl5->GetXaxis()->SetRangeUser(-15, 25);
	h_tmpl5->SetLineColor(kViolet);
	h_tmpl5->SetLineWidth(2);
	h_tmpl5->Draw("histlSAME");



	h_tmpl_120->SetLineColor(kOrange);
	h_tmpl_120->SetLineWidth(2);
	h_tmpl_120->SetLineStyle(2);
	h_tmpl_120->Draw("histlSAME");

	h_tmpl1_120->SetLineColor(kBlack);
	h_tmpl1_120->SetLineWidth(2);
	h_tmpl1_120->SetLineStyle(2);
	h_tmpl1_120->Draw("histlSAME");

	h_tmpl2_120->SetLineColor(kBlue);
	h_tmpl2_120->SetLineWidth(2);
	h_tmpl2_120->SetLineStyle(2);
	h_tmpl2_120->Draw("histlSAME");

	h_tmpl3_120->SetLineColor(kGreen+1);
	h_tmpl3_120->SetLineWidth(2);
	h_tmpl3_120->SetLineStyle(2);
	h_tmpl3_120->Draw("histlSAME");

	h_tmpl4_120->SetLineColor(kRed);
	h_tmpl4_120->SetLineWidth(2);
	h_tmpl4_120->SetLineStyle(2);
	h_tmpl4_120->Draw("histlSAME");

	h_tmpl5_120->SetLineColor(kViolet);
	h_tmpl5_120->SetLineWidth(2);
	h_tmpl5_120->SetLineStyle(2);
	h_tmpl5_120->Draw("histlSAME");



	legend->Draw("SAME");

	canv_output->cd();
	cDown->cd();
	cDown->SetGridy(); 

	TH1 *ratio1 = (TH1*)h_tmpl1->Clone();
	TH1 *ratio2 = (TH1*)h_tmpl2->Clone();
	TH1 *ratio3 = (TH1*)h_tmpl3->Clone();
	TH1 *ratio4 = (TH1*)h_tmpl4->Clone();
	TH1 *ratio5 = (TH1*)h_tmpl5->Clone();
	TH1 *ratio6 = (TH1*)h_tmpl->Clone();

	ratio1->Divide(h_tmpl1_120);
	ratio2->Divide(h_tmpl2_120);
	ratio3->Divide(h_tmpl3_120);
	ratio4->Divide(h_tmpl4_120);
	ratio5->Divide(h_tmpl5_120);
	ratio6->Divide(h_tmpl_120);

	ratio1->GetYaxis()->SetRangeUser(0.7, 1.3);
	ratio1->SetLineColor(kBlack);
	ratio1->SetLineWidth(2);
	ratio1->Draw("histl");
	ratio1->GetYaxis()->SetTitle("ratio");
	ratio1->GetYaxis()->SetTitleSize(0.15);
	ratio1->GetYaxis()->SetTitleOffset(0.3); 
	ratio1->SetLabelSize(0.12,"y");	
	ratio1->GetXaxis()->SetTitle("Time (ns)");
	ratio1->GetXaxis()->SetTickSize(0.2);
	ratio1->GetXaxis()->SetTitleSize(0.18);
	ratio1->SetLabelSize(0.15,"x");	
	ratio1->GetYaxis()->SetNdivisions(505);

	ratio2->SetLineColor(kBlue);
	ratio2->SetLineWidth(2);
	ratio2->Draw("histlSAME");

	ratio3->SetLineColor(kGreen+1);
	ratio3->SetLineWidth(2);
	ratio3->Draw("histlSAME");

	ratio6->SetLineColor(kOrange);
	ratio6->SetLineWidth(2);
	ratio6->Draw("histlSAME");

	ratio4->SetLineColor(kRed);
	ratio4->SetLineWidth(2);
	ratio4->Draw("histlSAME");

	ratio5->SetLineColor(kViolet);
	ratio5->SetLineWidth(2);
	ratio5->Draw("histlSAME");


	canv_output->cd();
	string outputDir = "~/www/TBJune2018_H4Analysis/";
	string outputFile = "ratio_shapes_C3_June2018_18C_160MHz_wrt_120MHz"; 


	canv_output->SaveAs( (outputDir+outputFile+".png").c_str() );
	canv_output->SaveAs( (outputDir+outputFile+".pdf").c_str() );
	canv_output->SaveAs( (outputDir+outputFile+".root").c_str() );

}
