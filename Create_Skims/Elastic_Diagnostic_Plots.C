/***********************************************************************************
 *
 * Author: Derek Holmberg
 *
 * Date Created: 10/30/2026
 *
 * Last Modified: 11/1/2026
 *
 * Purpose:
 * This program parses the ROOT file skims created from the gmn trains from Noemie's
 * elastic analysis code. This program loops over the output ROOT files with the 
 * fiducial cuts already applied and makes plots of several key kinematic quantities
 * to check for any differences in behavior between run periods. 
 *
***********************************************************************************/

#include "../Header_Files/Binning_Classes.h"

using namespace std;

// Object used to hold the plots for different data sets
struct DiagnosticPlots{

    DiagnosticPlots() = default;
    DiagnosticPlots( string nameExt ){	
	Ee_calc->SetName( string("Ee_calc"+nameExt).c_str() );
 	mm2_total->SetName( string("mm2_total"+nameExt).c_str() );
	dPhi->SetName( string("dPhi"+nameExt).c_str() );
	delta_pp->SetName( string("delta_pp"+nameExt).c_str() );
	delta_pe->SetName( string("delta_pe"+nameExt).c_str() );
	Q2->SetName( string("Q2"+nameExt).c_str() );
	W2->SetName( string("W2"+nameExt).c_str() );
	dPhi_vs_dpp->SetName( string("dPhi_vs_dpp"+nameExt).c_str() );
	dPhi_vs_dpe->SetName( string("dPhi_vs_dpe"+nameExt).c_str() );

    }

    // Various kinematic plots used for diagnostics
    TH1D* Ee_calc = new TH1D("Ee_calc", "E_{Beam} Distribution; E_{Beam} [GeV];",2000,0,20);
    TH1D* mm2_total = new TH1D("mm2_total", "MM_{ep}^{2} Distribution; MM_{ep}^{2} [GeV^{2}];",1000,-0.5,0.5);
    TH1D* dPhi = new TH1D("dPhi", "#Delta #phi Distribution; #Delta #phi [deg];",400,170,190);
    TH1D* delta_pp = new TH1D("delta_pp", "#delta P_{p} Distribution; #delta P_{p} [GeV/c];",2000,-1,1);
    TH1D* delta_pe = new TH1D("delta_pe", "#delta P_{e-} Distribution; #delta P_{e-} [Gev/c];",2000,-1,1);
    TH1D* Q2 = new TH1D("Q2", "Q^{2} Distribution; Q^{2} [GeV^{2}];",500,0,5);
    TH1D* W2 = new TH1D("W2", "W^{2} Distribution; W^{2} [GeV^{2}];",400,0,4);

    TH2D* dPhi_vs_dpp = new TH2D("dPhi_vs_dpp", "#Delta#phi vs. #delta P_{p} Distribution; #Delta P_{p} [GeV/c]; #Delta #phi [deg]",200,-1,1,40,170,190 );
    TH2D* dPhi_vs_dpe = new TH2D("dPhi_vs_dpe", "#Delta#phi vs. #delta P_{e-} Distribution; #Delta P_{e-} [GeV/c]; #Delta #phi [deg]",200,-1,1,40,170,190 );

};

// Takes the DiagnosticPlots objects and plots them together.
// First entry is Su22, then Fa22Neg, Fa22Pos, Sp23Inb.
TCanvas* PlotElasticDiagnosticPlots( vector<DiagnosticPlots>& Plots ){

    TCanvas* c = new TCanvas("c","c",1200,600); c->Divide(3,2);

    // Used for colors
    vector<int> colors = {1, 2, 3, 4};

    int RunPeriods = Plots.size();

    for(int i=0; i<RunPeriods; i++){
	c->cd(1);
	Plots[i].Ee_calc->SetLineColor( colors[i] ); Plots[i].Ee_calc->SetStats(false);
	Plots[i].Ee_calc->Scale( 1.0 / Plots[i].Ee_calc->Integral() );
	Plots[i].Ee_calc->GetXaxis()->SetRangeUser(8,12);
	Plots[i].Ee_calc->GetYaxis()->SetRangeUser(0,0.003);
	if( i==0 ) Plots[i].Ee_calc->Draw("hist"); else Plots[i].Ee_calc->Draw("hist same");

	c->cd(2);
	Plots[i].mm2_total->SetLineColor( colors[i] ); Plots[i].mm2_total->SetStats(false);
	Plots[i].mm2_total->Scale( 1.0 / Plots[i].mm2_total->Integral() );
	Plots[i].mm2_total->GetXaxis()->SetRangeUser(-0.05,0.05);
	Plots[i].mm2_total->GetYaxis()->SetRangeUser(0,0.01);
	if( i==0 ) Plots[i].mm2_total->Draw("hist"); else Plots[i].mm2_total->Draw("hist same");

	c->cd(3);
	Plots[i].dPhi->SetLineColor( colors[i] ); Plots[i].dPhi->SetStats(false);
	Plots[i].dPhi->Scale( 1.0 / Plots[i].dPhi->Integral() );
	Plots[i].dPhi->GetXaxis()->SetRangeUser(175,185);
	Plots[i].dPhi->GetYaxis()->SetRangeUser(0.002,0.006);
	if( i==0 ) Plots[i].dPhi->Draw("hist"); else Plots[i].dPhi->Draw("hist same");

	c->cd(4);
	Plots[i].delta_pp->SetLineColor( colors[i] ); Plots[i].delta_pp->SetStats(false);
	Plots[i].delta_pp->Scale( 1.0 / Plots[i].delta_pp->Integral() );
	if( i==0 ) Plots[i].delta_pp->Draw("hist"); else Plots[i].delta_pp->Draw("hist same");

	c->cd(5);
	Plots[i].delta_pe->SetLineColor( colors[i] ); Plots[i].delta_pe->SetStats(false);
	Plots[i].delta_pe->Scale( 1.0 / Plots[i].delta_pe->Integral() );
	if( i==0 ) Plots[i].delta_pe->Draw("hist"); else Plots[i].delta_pe->Draw("hist same");

	c->cd(6);
	Plots[i].Q2->SetLineColor( colors[i] ); Plots[i].Q2->SetStats(false);
	Plots[i].Q2->Scale( 1.0 / Plots[i].Q2->Integral() );
	if( i==0 ) Plots[i].Q2->Draw("hist"); else Plots[i].Q2->Draw("hist same");

    }

    TLegend* leg = new TLegend(0.65,0.6,0.9,0.9);
    leg->AddEntry( Plots[0].Ee_calc, "Su22",    "lp" );
    leg->AddEntry( Plots[1].Ee_calc, "Fa22Neg", "lp" );
    leg->AddEntry( Plots[2].Ee_calc, "Fa22Pos", "lp" );
    leg->AddEntry( Plots[3].Ee_calc, "Sp23Inb", "lp" );
    c->cd(1);
    leg->Draw("same");

    return c;

}

TCanvas* Plot2DHistos( vector<DiagnosticPlots>& Plots, string plotType ){

    string canName = "c"+plotType;
    TCanvas* c = new TCanvas(canName.c_str(),canName.c_str(),1600,1200); c->Divide(2,2);

    // Used for colors
    vector<int> colors = {1, 2, 3, 4};
    // Used for setting the title
    vector<string> title = {"Summer","Fall (Neg. Sol.)","Fall (Pos. Sol.)","Spring (Inb.)"};

    int RunPeriods = Plots.size();

    if( plotType == "Proton" ){
        for(int i=0; i<RunPeriods; i++){
	    c->cd(i+1);
	    //Plots[i].dPhi_vs_dpp->SetStats(false);
	    Plots[i].dPhi_vs_dpp->SetTitle( string(title[i]+"#Delta#phi vs. #delta P_{p} Distribution; #Delta P_{p} [GeV/c]; #Delta #phi [deg]").c_str() );
	    Plots[i].dPhi_vs_dpp->SetStats(false);
	    Plots[i].dPhi_vs_dpp->Draw("colz");
	    gPad->SetLogz();
        }
    }
    else if( plotType == "Electron" ){
        for(int i=0; i<RunPeriods; i++){
	    c->cd(i+1);
	    //Plots[i].dPhi_vs_dpe->SetStats(false);
	    Plots[i].dPhi_vs_dpe->SetTitle( string(title[i]+"#Delta#phi vs. #delta P_{e-} Distribution; #Delta P_{e-} [GeV/c]; #Delta #phi [deg]").c_str() );
	    Plots[i].dPhi_vs_dpe->SetStats(false);
	    Plots[i].dPhi_vs_dpe->Draw("colz");
	    gPad->SetLogz();
        }
    }

    return c;

}

void Elastic_Diagnostic_Plots(){

    DiagnosticPlots Su22_Plots("Su22");
    DiagnosticPlots Fa22Neg_Plots("Fa22Neg");
    DiagnosticPlots Fa22Pos_Plots("Fa22Pos");
    DiagnosticPlots Sp23Inb_Plots("Sp23");

    RunPeriod RP;

    //string fileExt = "_";
    string fileExt = "_correctedP_";

    for(int i=1; i<24; i++){ //Loop over NH3 epochs
	
      string epoch = "P";
      if( i<10 ) epoch += "0"+to_string(i);
      else epoch += to_string(i);

      vector<int> Epoch_Runs = RP.getRunEpoch( epoch );

      for(int Run : Epoch_Runs){
    //for(int j=16100; j<18000; j++){

        string filePath = "../Latest_Skims/Diagnostic_Elastic/Elastic_Diagnostic"+ fileExt + to_string(Run) +".root";
        ifstream fin( filePath.c_str() );
        if( !fin.fail() ){
	    TFile* file = new TFile( filePath.c_str(), "read" );
	    auto Ee_calc     = (TH1D*)file->Get("Ee_calc");
	    auto mm2_total   = (TH1D*)file->Get("mm2_total");
	    auto dPhi        = (TH1D*)file->Get("dPhi");
	    auto delta_pp    = (TH1D*)file->Get("delta_pp");
	    auto delta_pe    = (TH1D*)file->Get("delta_pe");
	    auto Q2          = (TH1D*)file->Get("Q2");
	    auto W2          = (TH1D*)file->Get("W2");
	    auto dPhi_vs_dpp = (TH2D*)file->Get("dPhi_vs_dpp");
	    auto dPhi_vs_dpe = (TH2D*)file->Get("dPhi_vs_dpe");

	    if( Run < 16850 ){
		Su22_Plots.Ee_calc->Add( Ee_calc ); Su22_Plots.mm2_total->Add( mm2_total ); Su22_Plots.dPhi->Add( dPhi );
		Su22_Plots.delta_pp->Add( delta_pp ); Su22_Plots.delta_pe->Add( delta_pe );
		Su22_Plots.Q2->Add( Q2 ); Su22_Plots.W2->Add( W2 );
		Su22_Plots.dPhi_vs_dpp->Add( dPhi_vs_dpp ); Su22_Plots.dPhi_vs_dpe->Add( dPhi_vs_dpe );
	    }
	    else if( Run > 16850 && Run < 17185 ){
		Fa22Neg_Plots.Ee_calc->Add( Ee_calc ); Fa22Neg_Plots.mm2_total->Add( mm2_total ); Fa22Neg_Plots.dPhi->Add( dPhi );
		Fa22Neg_Plots.delta_pp->Add( delta_pp ); Fa22Neg_Plots.delta_pe->Add( delta_pe );
		Fa22Neg_Plots.Q2->Add( Q2 ); Fa22Neg_Plots.W2->Add( W2 );
		Fa22Neg_Plots.dPhi_vs_dpp->Add( dPhi_vs_dpp ); Fa22Neg_Plots.dPhi_vs_dpe->Add( dPhi_vs_dpe );
	    }
	    else if( Run > 17185 && Run < 17450 ){
		Fa22Pos_Plots.Ee_calc->Add( Ee_calc ); Fa22Pos_Plots.mm2_total->Add( mm2_total ); Fa22Pos_Plots.dPhi->Add( dPhi );
		Fa22Pos_Plots.delta_pp->Add( delta_pp ); Fa22Pos_Plots.delta_pe->Add( delta_pe );
		Fa22Pos_Plots.Q2->Add( Q2 ); Fa22Pos_Plots.W2->Add( W2 );
		Fa22Pos_Plots.dPhi_vs_dpp->Add( dPhi_vs_dpp ); Fa22Pos_Plots.dPhi_vs_dpe->Add( dPhi_vs_dpe );
	    }
	    else if( Run> 17450 && Run <= 17768 ){
		Sp23Inb_Plots.Ee_calc->Add( Ee_calc ); Sp23Inb_Plots.mm2_total->Add( mm2_total ); Sp23Inb_Plots.dPhi->Add( dPhi );
		Sp23Inb_Plots.delta_pp->Add( delta_pp ); Sp23Inb_Plots.delta_pe->Add( delta_pe );
		Sp23Inb_Plots.Q2->Add( Q2 ); Sp23Inb_Plots.W2->Add( W2 );
		Sp23Inb_Plots.dPhi_vs_dpp->Add( dPhi_vs_dpp ); Sp23Inb_Plots.dPhi_vs_dpe->Add( dPhi_vs_dpe );
	    }

        }
        fin.close();
      } // End loop over this individual epoch
    } // End loop over ALL NH3 epochs

    vector<DiagnosticPlots> Plots = { Su22_Plots, Fa22Neg_Plots, Fa22Pos_Plots, Sp23Inb_Plots };

    TCanvas* All_Plots = PlotElasticDiagnosticPlots( Plots ); All_Plots->Print( string("PDF_Plots/Elastic_Diagnostic"+ fileExt +"Plots.pdf(").c_str() );
    TCanvas* Dpp = Plot2DHistos( Plots, "Proton" );   Dpp->Print( string("PDF_Plots/Elastic_Diagnostic"+ fileExt +"Plots.pdf").c_str() );
    TCanvas* Dpe = Plot2DHistos( Plots, "Electron" ); Dpe->Print( string("PDF_Plots/Elastic_Diagnostic"+ fileExt +"Plots.pdf)").c_str() );

}

