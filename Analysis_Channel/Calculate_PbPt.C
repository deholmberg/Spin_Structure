/***********************************************************************************
 *
 * Author: Derek Holmberg
 *
 * Date Created: 9/25/2026
 *
 * Last Modified: 9/25/2026
 *
 * Purpose:
 * This program calculates the DIS values of PbPt using the DFs from "Test.C".
 * The program reads in files from "Output_Data/All_DF_Data_Epoch_X.txt" and loads
 * them into a DataSet class. For each run in that epoch, a unique DataSet class 
 * with the appropriate DF values is created, and PbPt is calculated for that run.
 * These PbPt are plotted as a function of run number.
 *
***********************************************************************************/

#include "../Header_Files/Binning_Classes.h"

using namespace std;

// This function returns a TCanvas object holding the normalized and unnormalized PbPt vs. run number 
// for a specified range of epochs
TCanvas* Plot_PbPt_Epoch( RunPeriod& Period, int FirstEpoch, int LastEpoch, string Period_Title ){

    // Hold all the relevant PbPt info
    // A_FC is the FC charge asymmetry for that run (with sign correction to TPol), and
    // Delta_PbPt is the difference between the PbPt with and without FC-normalization
    vector<double> PbPt, PbPtErr, UnnormPbPt, UnnormPbPtErr, Runs, A_FC, Delta_PbPt, Delta_PbPtErr;
    // The positive targets
    vector<double> Pos_PbPt, Pos_PbPtErr, Pos_UnnormPbPt, Pos_UnnormPbPtErr, Pos_Runs, Pos_A_FC, Pos_Delta_PbPt, Pos_Delta_PbPtErr;
    // The negative targets
    vector<double> Neg_PbPt, Neg_PbPtErr, Neg_UnnormPbPt, Neg_UnnormPbPtErr, Neg_Runs, Neg_A_FC, Neg_Delta_PbPt, Neg_Delta_PbPtErr;

    vector<double> NMR_PbPt, NMR_PbPtErr; // Offline NMR values
    vector<double> Pos_NMR_PbPt, Pos_NMR_PbPtErr, Neg_NMR_PbPt, Neg_NMR_PbPtErr;
    
    // Track any missed runs
    vector<int> MissedRuns;

    // Start the loop over all epochs in the specified range
    for( int i=FirstEpoch; i<=LastEpoch; i++ ){
	string epoch;
	if( i < 10 ) epoch = "P0"+to_string(i);
	else         epoch = "P" +to_string(i);

	vector<int> EpochRuns = Period.getRunEpoch( epoch );
	for( int Run : EpochRuns ){
	    DataSet thisRun;
	    // Read in the DF data
	    string DFdata = "Output_Data/All_DF_Data_Epoch_"+to_string(i)+".txt";
	    thisRun.ReadDFfromTXT( DFdata );
	    if( !thisRun.AddToBins( string("../Latest_Skims/Text_Files/NH3_"+to_string(Run)+"_DF_Data.txt"), "NH3" ) ) MissedRuns.push_back( Run );
	    
	    // Get the unnormalized PbPt
	    thisRun.MaxLikelihoodPbPt( Run, Period, false );
	    double unnormPbPt    = thisRun.getBT_Pol();
	    double unnormPbPtErr = thisRun.getBT_Pol_Err();
	    double afc  	 = thisRun.getAFC("NH3");

	    // Get the normalized PbPt
	    thisRun.MaxLikelihoodPbPt( Run, Period, true );
	    double normPbPt    = thisRun.getBT_Pol();
	    double normPbPtErr = thisRun.getBT_Pol_Err();

	    // Add to the vectors for plotting later
	    Runs.push_back( Run );
	    PbPt.push_back( normPbPt ); PbPtErr.push_back( normPbPtErr );
	    UnnormPbPt.push_back( unnormPbPt ); UnnormPbPtErr.push_back( unnormPbPtErr );
	    A_FC.push_back( afc );
	    Delta_PbPt.push_back( (unnormPbPt - normPbPt)/unnormPbPt );

	    // Add the offline NMR values
	    NMR_PbPt.push_back( Period.getOfflineTPol( Run ) );
	    NMR_PbPtErr.push_back( Period.getOfflineTPolErr( Run ) );

	    // Positive runs
	    if( NMR_PbPt > 0 ){
	        // Add to the vectors for plotting later
	        Pos_Runs.push_back( Run );
	        Pos_PbPt.push_back( normPbPt ); Pos_PbPtErr.push_back( normPbPtErr );
	        Pos_UnnormPbPt.push_back( unnormPbPt ); Pos_UnnormPbPtErr.push_back( unnormPbPtErr );
	        Pos_A_FC.push_back( afc );
	        Pos_Delta_PbPt.push_back( (unnormPbPt - normPbPt)/unnormPbPt );

	        // Add the offline NMR values
	        Pos_NMR_PbPt.push_back( Period.getOfflineTPol( Run ) );
	        Pos_NMR_PbPtErr.push_back( Period.getOfflineTPolErr( Run ) );
	    }
	    else if( NMR_PbPt < 0 ){
	        // Add to the vectors for plotting later
	        Neg_Runs.push_back( Run );
	        Neg_PbPt.push_back( normPbPt ); Neg_PbPtErr.push_back( normPbPtErr );
	        Neg_UnnormPbPt.push_back( unnormPbPt ); Neg_UnnormPbPtErr.push_back( unnormPbPtErr );
	        Neg_A_FC.push_back( afc );
	        Neg_Delta_PbPt.push_back( (unnormPbPt - normPbPt)/unnormPbPt );

	        // Add the offline NMR values
	        Neg_NMR_PbPt.push_back( Period.getOfflineTPol( Run ) );
	        Neg_NMR_PbPtErr.push_back( Period.getOfflineTPolErr( Run ) );
	    }

	}
	
    }

    // Kill the attempt if runs are missing
    if( MissedRuns.size() > 0 ){
	cout <<"ERROR: Missed runs were found! Aborting calculation...\n";
	cout <<"Runs Missing:\n";
	for(int r : MissedRuns) cout <<"---> "<< r << endl;
	TCanvas* blank;
	return blank;
    }

    // Make plots if everything is OK
/*
    TLegend* legUnnorm    = new TLegend(0.18,0.8,0.48,0.9);
    TLegend* legNorm      = new TLegend(0.18,0.8,0.48,0.9);
    TMultiGraph* mgUnnorm = new TMultiGraph();
    TMultiGraph* mgNorm   = new TMultiGraph();

    // Plot the offline NMR values
    TGraphErrors* grNMR = new TGraphErrors( Runs.size(), Runs.data(), NMR_PbPt.data(), nullptr, NMR_PbPtErr.data() );
    grNMR->SetMarkerStyle( kFullSquare );
    grNMR->SetMarkerColor( kBlack );
    mgNorm->Add( grNMR, "p" );
    mgUnnorm->Add( grNMR, "p" );
    legNorm->AddEntry( grNMR, "Offline NMR P_{B}P_{T}", "p" );
    legUnnorm->AddEntry( grNMR, "Offline NMR P_{B}P_{T}", "p" );

    // Plot Normalized PbPt
    TGraphErrors* grNormPbPt = new TGraphErrors( Runs.size(), Runs.data(), PbPt.data(), nullptr, PbPtErr.data() );
    grNormPbPt->SetMarkerStyle( kFullCircle );
    grNormPbPt->SetMarkerColor( kBlue );
    mgNorm->Add( grNormPbPt, "p" );
    legNorm->AddEntry( grNormPbPt, "FC-Norm. P_{B}P_{T}", "p" );

    mgNorm->SetTitle( string( Period_Title+" P_{B}P_{T} (Norm. to FC Charge); Run Number; P_{B}P_{T}").c_str() );

    // Plot Unnormalized PbPt
    TGraphErrors* grUnnormPbPt = new TGraphErrors( Runs.size(), Runs.data(), UnnormPbPt.data(), nullptr, UnnormPbPtErr.data() );
    grUnnormPbPt->SetMarkerStyle( kFullCircle );
    grUnnormPbPt->SetMarkerColor( kViolet );
    mgUnnorm->Add( grUnnormPbPt, "p" );
    legUnnorm->AddEntry( grUnnormPbPt, "Unnorm. P_{B}P_{T}", "p" );

    mgUnnorm->SetTitle( string( Period_Title+" P_{B}P_{T} (Unnorm. to FC Charge); Run Number; P_{B}P_{T}").c_str() );
*/
    TLegend* Pos_legUnnorm    = new TLegend(0.18,0.8,0.48,0.9);
    TLegend* Pos_legNorm      = new TLegend(0.18,0.8,0.48,0.9);
    TMultiGraph* Pos_mgUnnorm = new TMultiGraph();
    TMultiGraph* Pos_mgNorm   = new TMultiGraph();

    // Plot the offline NMR values
    TGraphErrors* Pos_grNMR = new TGraphErrors( Pos_Runs.size(), Pos_Runs.data(), Pos_NMR_PbPt.data(), nullptr, Pos_NMR_PbPtErr.data() );
    Pos_grNMR->SetMarkerStyle( kFullSquare );
    Pos_grNMR->SetMarkerColor( kBlack );
    Pos_mgNorm->Add( Pos_grNMR, "p" );
    Pos_mgUnnorm->Add( Pos_grNMR, "p" );
    Pos_legNorm->AddEntry( Pos_grNMR, "Offline NMR P_{B}P_{T}", "p" );
    Pos_legUnnorm->AddEntry( Pos_grNMR, "Offline NMR P_{B}P_{T}", "p" );

    // Plot Normalized PbPt
    TGraphErrors* Pos_grNormPbPt = new TGraphErrors( Pos_Runs.size(), Pos_Runs.data(), Pos_PbPt.data(), nullptr, Pos_PbPtErr.data() );
    Pos_grNormPbPt->SetMarkerStyle( kFullCircle );
    Pos_grNormPbPt->SetMarkerColor( kBlue );
    Pos_mgNorm->Add( Pos_grNormPbPt, "p" );
    Pos_legNorm->AddEntry( Pos_grNormPbPt, "FC-Norm. P_{B}P_{T}", "p" );

    Pos_mgNorm->SetTitle( string( Period_Title+" P_{B}P_{T} (Norm. to FC Charge); Run Number; P_{B}P_{T}").c_str() );

    // Plot Unnormalized PbPt
    TGraphErrors* grUnnormPbPt = new TGraphErrors( Runs.size(), Runs.data(), UnnormPbPt.data(), nullptr, UnnormPbPtErr.data() );
    grUnnormPbPt->SetMarkerStyle( kFullCircle );
    grUnnormPbPt->SetMarkerColor( kViolet );
    mgUnnorm->Add( grUnnormPbPt, "p" );
    legUnnorm->AddEntry( grUnnormPbPt, "Unnorm. P_{B}P_{T}", "p" );

    mgUnnorm->SetTitle( string( Period_Title+" P_{B}P_{T} (Unnorm. to FC Charge); Run Number; P_{B}P_{T}").c_str() );


    // Plot the %-difference between normalized PbPt and unnormalized PbPt versus FC asymmetry
    TGraphErrors* deltaFC = new TGraphErrors( A_FC.size(), A_FC.data(), Delta_PbPt.data(), nullptr, nullptr );
    deltaFC->SetMarkerStyle( kFullCircle );
    deltaFC->SetMarkerColor( kRed );
    deltaFC->SetTitle( string( Period_Title+" Percent Diff. in P_{B}P_{T} versus A_{FC}; A_{FC}; Percent Difference").c_str() );

    TCanvas* c = new TCanvas("c","c",2200,600);
    c->Divide(3,1);

    c->cd(1); 
    gPad->SetLeftMargin(0.18); gPad->SetRightMargin(0.02);
    mgNorm->Draw("ap"); legNorm->Draw("same");

    c->cd(2); 
    gPad->SetLeftMargin(0.18); gPad->SetRightMargin(0.02);
    mgUnnorm->Draw("ap"); legUnnorm->Draw("same");

    c->cd(3); 
    gPad->SetLeftMargin(0.18); gPad->SetRightMargin(0.02);
    deltaFC->Draw("ap");

    return c;

}

void Calculate_PbPt(){

    RunPeriod Period; // Holds the run info

    auto Plot_Su22 = Plot_PbPt_Epoch( Period, 1, 10, "Summer" );
/*
    // Hold all the relevant PbPt info
    // A_FC is the FC charge asymmetry for that run (with sign correction to TPol), and
    // Delta_PbPt is the difference between the PbPt with and without FC-normalization
    vector<double> PbPt, PbPtErr, UnnormPbPt, UnnormPbPtErr, Runs, A_FC, Delta_PbPt, Delta_PbPtErr;
    
    // Track any missed runs
    vector<int> MissedRuns;

    // Start the loop over all epochs in Summer 2022
    for( int i=1; i<=2; i++ ){
	string epoch;
	if( i < 10 ) epoch = "P0"+to_string(i);
	else         epoch = "P" +to_string(i);

	vector<int> EpochRuns = Period.getRunEpoch( epoch );
	for( int Run : EpochRuns ){
	    DataSet thisRun;
	    // Read in the DF data
	    string DFdata = "Output_Data/All_DF_Data_Epoch_"+to_string(i)+".txt";
	    thisRun.ReadDFfromTXT( DFdata );
	    if( !thisRun.AddToBins( string("../Latest_Skims/Text_Files/NH3_"+to_string(Run)+"_DF_Data.txt"), "NH3" ) ) MissedRuns.push_back( Run );
	    
	    // Get the unnormalized PbPt
	    thisRun.MaxLikelihoodPbPt( Run, Period, false );
	    double unnormPbPt    = thisRun.getBT_Pol();
	    double unnormPbPtErr = thisRun.getBT_Pol_Err();
	    double afc  	 = thisRun.getAFC("NH3");

	    // Get the normalized PbPt
	    thisRun.MaxLikelihoodPbPt( Run, Period, true );
	    double normPbPt    = thisRun.getBT_Pol();
	    double normPbPtErr = thisRun.getBT_Pol_Err();

	    // Add to the vectors for plotting later
	    Runs.push_back( Run );
	    PbPt.push_back( normPbPt ); PbPtErr.push_back( normPbPtErr );
	    UnnormPbPt.push_back( unnormPbPt ); UnnormPbPtErr.push_back( unnormPbPtErr );
	    A_FC.push_back( afc );
	    Delta_PbPt.push_back( (unnormPbPt - normPbPt)/unnormPbPt );
	}
	
    }

    // Kill the attempt if runs are missing
    if( MissedRuns.size() > 0 ){
	cout <<"ERROR: Missed runs were found! Aborting calculation...\n";
	cout <<"Runs Missing:\n";
	for(int r : MissedRuns) cout <<"---> "<< r << endl;
	return;
    }

    // Make plots if everything is OK

    // Plot Normalized PbPt
    TGraphErrors* grNormPbPt = new TGraphErrors( Runs.size(), Runs.data(), PbPt.data(), nullptr, PbPtErr.data() );
    grNormPbPt->SetMarkerStyle( kFullCircle );
    grNormPbPt->SetMarkerColor( kBlue );
    grNormPbPt->SetTitle("Summer P_{B}P_{T} (Norm. to FC Charge); Run Number; P_{B}P_{T}");

    // Plot Unnormalized PbPt
    TGraphErrors* grUnnormPbPt = new TGraphErrors( Runs.size(), Runs.data(), UnnormPbPt.data(), nullptr, UnnormPbPtErr.data() );
    grUnnormPbPt->SetMarkerStyle( kFullCircle );
    grUnnormPbPt->SetMarkerColor( kViolet );
    grUnnormPbPt->SetTitle("Summer P_{B}P_{T} (Unnorm. to FC Charge); Run Number; P_{B}P_{T}");

    // Plot the %-difference between normalized PbPt and unnormalized PbPt versus FC asymmetry
    TGraphErrors* deltaFC = new TGraphErrors( A_FC.size(), A_FC.data(), Delta_PbPt.data(), nullptr, nullptr );
    deltaFC->SetMarkerStyle( kFullCircle );
    deltaFC->SetMarkerColor( kRed );
    deltaFC->SetTitle("Summer Percent Diff. in P_{B}P_{T} versus A_{FC}; A_{FC}; Percent Difference");

    TCanvas* c = new TCanvas("c","c",2200,600);
    c->Divide(3,1);

    c->cd(1); grNormPbPt->Draw("ap");
    c->cd(2); grUnnormPbPt->Draw("ap");
    c->cd(3); deltaFC->Draw("ap");
*/
/*
    DataSet Test;
    Test.ReadDFfromTXT("Output_Data/All_DF_Data_Epoch_10.txt");
    //Test.Print(true);
    if( Test.AddToBins("../Latest_Skims/Text_Files/NH3_16772_DF_Data.txt", "NH3") ) cout << "Opened run 16772!\n";

    // Calculate it twice to make sure nothing is overwritten upon each call
    Test.MaxLikelihoodPbPt( 16772, Period, true );
    Test.MaxLikelihoodPbPt( 16772, Period, true );

    // Calculate it a couple times without FC charge normalization
    Test.MaxLikelihoodPbPt( 16772, Period, false );
    Test.MaxLikelihoodPbPt( 16772, Period, false );
*/

}
