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


void drawTmpl_En()
{
	// TString name = "shape_tmpl_11409_150GeV_C3_T";
	// TString name1 = "shape_tmpl_11433_25GeV_C3_T";
	// TString name2 = "shape_tmpl_11423_50GeV_C3_T";
	// TString name3 = "shape_tmpl_11416_100GeV_C3_T";
	// TString name4 = "shape_tmpl_11393_200GeV_C3_T";
	// TString name5 = "shape_tmpl_11357_250GeV_C3_T";

	TString name = "shape_tmpl_11478_150GeV_C3_T";
	TString name1 = "shape_tmpl_11456_25GeV_C3_T";
	TString name2 = "shape_tmpl_11468_50GeV_C3_T";
	TString name3 = "shape_tmpl_11472_100GeV_C3_T";
	TString name4 = "shape_tmpl_11485_200GeV_C3_T";
	TString name5 = "shape_tmpl_11494_250GeV_C3_T";

	TFile inputFile("shapes/"+name+"_TH1.root");
	TFile inputFile1("shapes/"+name1+"_TH1.root");
	TFile inputFile2("shapes/"+name2+"_TH1.root");
	TFile inputFile3("shapes/"+name3+"_TH1.root");
	TFile inputFile4("shapes/"+name4+"_TH1.root");
	TFile inputFile5("shapes/"+name5+"_TH1.root");

	TCanvas* canv = (TCanvas*)inputFile.Get("outputCanv");
	TCanvas* canv1 = (TCanvas*)inputFile1.Get("outputCanv");
	TCanvas* canv2 = (TCanvas*)inputFile2.Get("outputCanv");
	TCanvas* canv3 = (TCanvas*)inputFile3.Get("outputCanv");
	TCanvas* canv4 = (TCanvas*)inputFile4.Get("outputCanv");
	TCanvas* canv5 = (TCanvas*)inputFile5.Get("outputCanv");

	TH1* h_tmpl=(TH1*)canv->FindObject("prof");
	TH1* h_tmpl1=(TH1*)canv1->FindObject("prof");
	TH1* h_tmpl2=(TH1*)canv2->FindObject("prof");
	TH1* h_tmpl3=(TH1*)canv3->FindObject("prof");
	TH1* h_tmpl4=(TH1*)canv4->FindObject("prof");
	TH1* h_tmpl5=(TH1*)canv5->FindObject("prof");

	h_tmpl->SetTitle("");
	h_tmpl1->SetTitle("");
	h_tmpl2->SetTitle("");
	h_tmpl3->SetTitle("");
	h_tmpl4->SetTitle("");
	h_tmpl5->SetTitle("");


	gStyle->SetOptStat(0);

	TCanvas *canv_output = new TCanvas("canv_output","canv_output", 900, 900);

	// TLegend* legend = new TLegend(0.6, 0.55, 0.8, 0.85);
	TLegend* legend = new TLegend(0.6, 0.5, 0.8, 0.8);
	legend->AddEntry(h_tmpl1, "25 GeV");
	legend->AddEntry(h_tmpl2, "50 GeV");	
	legend->AddEntry(h_tmpl3, "100 GeV");	
	legend->AddEntry(h_tmpl, "150 GeV");
	legend->AddEntry(h_tmpl4, "200 GeV");
	legend->AddEntry(h_tmpl5, "250 GeV");
	legend->SetFillColor( kWhite ) ; 
	legend->SetLineColor( kWhite ) ; 

	gPad->SetTickx();
	gPad->SetTicky();

	h_tmpl->GetYaxis()->SetRangeUser(-0.1, 1.05);
	h_tmpl->SetLineColor(kOrange);
	h_tmpl->SetLineWidth(2);
	h_tmpl->Draw("histl");
	h_tmpl->GetYaxis()->SetTitle("Normalized amplitude (A. U.)");
	h_tmpl->GetYaxis()->SetTitleSize(0.04);
	h_tmpl->GetYaxis()->SetLabelSize(0.03);	
	h_tmpl->GetYaxis()->SetTitleOffset(0.999); 
	h_tmpl->GetXaxis()->SetTitle("time (ns)");
	h_tmpl->GetXaxis()->SetTitleSize(0.04);
	h_tmpl->SetLabelSize(0.03,"x");	

	h_tmpl1->GetYaxis()->SetRangeUser(-0.1, 1.05);
	h_tmpl1->SetLineColor(kBlack);
	h_tmpl1->SetLineWidth(2);
	h_tmpl1->Draw("histlSAME");

	h_tmpl2->GetYaxis()->SetRangeUser(-0.1, 1.05);
	h_tmpl2->SetLineColor(kBlue);
	h_tmpl2->SetLineWidth(2);
	h_tmpl2->Draw("histlSAME");

	h_tmpl3->GetYaxis()->SetRangeUser(-0.1, 1.05);
	h_tmpl3->SetLineColor(kGreen+1);
	h_tmpl3->SetLineWidth(2);
	h_tmpl3->Draw("histlSAME");

	h_tmpl4->GetYaxis()->SetRangeUser(-0.1, 1.05);
	h_tmpl4->SetLineColor(kRed);
	h_tmpl4->SetLineWidth(2);
	h_tmpl4->Draw("histlSAME");

	h_tmpl5->GetYaxis()->SetRangeUser(-0.1, 1.05);
	h_tmpl5->SetLineColor(kViolet);
	h_tmpl5->SetLineWidth(2);
	h_tmpl5->Draw("histlSAME");

	legend->Draw("SAME");
	

	TLatex latex;
    latex.SetNDC();
    latex.SetTextAngle(0);
    latex.SetTextColor(kBlack);

    float l = canv_output->GetLeftMargin();
    float t = canv_output->GetTopMargin(); 
    float r = canv_output->GetRightMargin();
    float b = canv_output->GetBottomMargin();
    float x = l + 0.045*(1-l-r);
    float y = 1 - t + 0.01; 

    latex.SetTextAlign(11); //31 on the right
 
    latex.SetTextFont(61);
    latex.SetTextSize(0.45*t);
    latex.DrawLatex(x,y,"CMS");

    latex.SetTextFont(52);
    latex.SetTextSize(0.4*t);  
    latex.DrawLatex(x + 0.1, y, "Preliminary");

    latex.SetTextFont(42);
    latex.SetTextSize(0.35*t);
    // latex.DrawLatex(1-r-0.3,y-0.08,"160 MHz - 18#circ C");
    latex.DrawLatex(1-r-0.3,y-0.08,"120 MHz - 18#circ C");

	canv_output->cd();
	string outputDir = "~/www/TBJune2018_H4Analysis/hodo/templates/";
	// string outputFile = "shapes_tmpl_C3_June2018_160MHz_18C"; 
	string outputFile = "shapes_tmpl_C3_June2018_120MHz_18C"; 

	canv_output->SaveAs( (outputDir+outputFile+".png").c_str() );
	canv_output->SaveAs( (outputDir+outputFile+".pdf").c_str() );
	canv_output->SaveAs( (outputDir+outputFile+".root").c_str() );

}
