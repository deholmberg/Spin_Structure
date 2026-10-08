/***********************************************************************************
 *
 * Author: Derek Holmberg
 *
 * Date Created: 9/25/2026
 *
 * Last Modified: 9/29/2026
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

// This function appends a string object to a file of a given path
void AppendToFile( string fileName, string data ){

    ofstream fout;
    fout.open( fileName, ios_base::app ); // Set the file to append
    fout << data;

}

// This function takes two vectors, one of values and one of errors and calculates a weighted average.
// This value is returned in a vector, where the first entry is the average and the second is the error.
vector<double> WeightedAverage( vector<double>& Vals, vector<double>& Errs ){

    vector<double> Avg = {0,0};

    // Kill attempt if vectors aren't same size
    int sizeVal = Vals.size();
    int sizeErr = Errs.size();

    if( sizeVal != sizeErr ){
	cout <<"ERROR: Value and error vectors aren't same size. No weighted average computed.\n";
	return Avg;
    }

    // Used for weighted averages
    double num = 0; double den = 0;
    for( int i=0; i<sizeVal; i++ ){
	
	if( Errs[i] != 0 ){
	    num += Vals[i] / pow( Errs[i], 2 );
	    den += 1.0 / pow( Errs[i], 2 );
	}
	else{
	    cout <<"ERROR: Invalid (zero) value for the error. No weighted average computed.\n";
	    return Avg;
	}

    }

    if( den > 0 ){
	Avg[0] = num / den;
	Avg[1] = sqrt( 1.0 / den );
    }
    else{
        cout <<"ERROR: No weighted average computed.\n";
    }

    return Avg;

}

// This function returns a TCanvas object holding the normalized and unnormalized PbPt vs. run number 
// for a specified range of epochs
TCanvas* Plot_PbPt_Epoch( RunPeriod& Period, int FirstEpoch, int LastEpoch, string Period_Title ){

    // Set the Xaxis run ranges based on the run period
    double runMin = 0; double runMax = 0;
    if( Period_Title == "Summer" ){ runMin = 16000; runMax = 16800; }
    if( Period_Title == "Fall (Neg. Sol.)" ){ runMin = 16800; runMax = 17250; }
    if( Period_Title == "Fall (Pos. Sol.)" ){ runMin = 17150; runMax = 17500; }
    if( Period_Title == "Spring" ){ runMin = 17450; runMax = 17800; }

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
	    double normPbPt     = thisRun.getBT_Pol();
	    double normPbPtErr  = thisRun.getBT_Pol_Err();

	    double deltaPbPt    = (unnormPbPt - normPbPt)/unnormPbPt;
	    double deltaPbPtErr = 0.5*sqrt( pow(unnormPbPtErr,2) + pow(normPbPtErr,2) );

	    // Add to the vectors for plotting later
	    Runs.push_back( Run );
	    PbPt.push_back( normPbPt ); PbPtErr.push_back( normPbPtErr );
	    UnnormPbPt.push_back( unnormPbPt ); UnnormPbPtErr.push_back( unnormPbPtErr );
	    A_FC.push_back( afc );
	    Delta_PbPt.push_back( deltaPbPt );
	    Delta_PbPtErr.push_back( deltaPbPtErr );

	    // Add the offline NMR values
	    double nmrPbPt = Period.getOfflineTPol( Run );
	    double nmrPbPtErr = Period.getOfflineTPolErr( Run );
	    NMR_PbPt.push_back( nmrPbPt );
	    NMR_PbPtErr.push_back( nmrPbPtErr );

	    // Positive runs
	    if( nmrPbPt > 0 ){
	        // Add to the vectors for plotting later
	        Pos_Runs.push_back( Run );
	        Pos_PbPt.push_back( normPbPt ); Pos_PbPtErr.push_back( normPbPtErr );
	        Pos_UnnormPbPt.push_back( unnormPbPt ); Pos_UnnormPbPtErr.push_back( unnormPbPtErr );
	        Pos_A_FC.push_back( afc );
	        Pos_Delta_PbPt.push_back( deltaPbPt );
		Pos_Delta_PbPtErr.push_back( deltaPbPtErr );

	        // Add the offline NMR values
	        Pos_NMR_PbPt.push_back( nmrPbPt );
	        Pos_NMR_PbPtErr.push_back( nmrPbPtErr );

	    }
	    else if( nmrPbPt < 0 ){
	        // Add to the vectors for plotting later
	        Neg_Runs.push_back( Run );
	        Neg_PbPt.push_back( normPbPt ); Neg_PbPtErr.push_back( normPbPtErr );
	        Neg_UnnormPbPt.push_back( unnormPbPt ); Neg_UnnormPbPtErr.push_back( unnormPbPtErr );
	        Neg_A_FC.push_back( afc );
	        Neg_Delta_PbPt.push_back( deltaPbPt );
		Neg_Delta_PbPtErr.push_back( deltaPbPtErr );

	        // Add the offline NMR values
	        Neg_NMR_PbPt.push_back( nmrPbPt );
	        Neg_NMR_PbPtErr.push_back( nmrPbPtErr );

	    }

	    // Write the PbPt data to a file
	    string dataStream = to_string( Run ) +"   "+ to_string( unnormPbPt ) +"   "+ to_string( unnormPbPtErr ) +"   ";
	    dataStream += to_string( normPbPt ) +"   "+ to_string( normPbPtErr ) + "\n";
	    AppendToFile( "Output_Data/All_NH3_PbPt.txt", dataStream );

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

    // Plot the %-difference between normalized PbPt and unnormalized PbPt versus FC asymmetry
    TGraphErrors* deltaFC = new TGraphErrors( A_FC.size(), A_FC.data(), Delta_PbPt.data(), nullptr, nullptr );
    deltaFC->SetMarkerStyle( kFullCircle );
    deltaFC->SetMarkerColor( kRed );
    deltaFC->SetTitle( string( Period_Title+" Percent Diff. in P_{B}P_{T} versus A_{FC}; A_{FC}; Percent Difference").c_str() );

*/
 
    string canName = Period_Title+"_c";
    TCanvas* c = new TCanvas(canName.c_str(),canName.c_str(),2200,1200);
    c->Divide(3,2);

    // Used for the fits
    TFitResultPtr r;

    // Positive plots
    TLegend* Pos_legUnnorm    = new TLegend(0.18,0.7,0.7,0.9); 
    Pos_legUnnorm->SetTextFont(43);
    Pos_legUnnorm->SetTextSize(15);
    Pos_legUnnorm->SetBorderSize(0);
    Pos_legUnnorm->SetFillStyle(0);
    TLegend* Pos_legNorm      = new TLegend(0.18,0.7,0.7,0.9); 
    Pos_legNorm->SetTextFont(43);
    Pos_legNorm->SetTextSize(15);
    Pos_legNorm->SetBorderSize(0);
    Pos_legNorm->SetFillStyle(0);
    TLegend* Pos_legDelta     = new TLegend(0.18,0.7,0.7,0.9);
    Pos_legDelta->SetTextFont(43);
    Pos_legDelta->SetTextSize(15);
    Pos_legDelta->SetBorderSize(0);
    Pos_legDelta->SetFillStyle(0);
    Pos_legDelta->SetHeader("Fit = [0] + [1]*A_{FC}","C");

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

    Pos_mgNorm->SetTitle( string( "Positive Target "+Period_Title+" P_{B}P_{T} (Norm. to FC Charge); Run Number; P_{B}P_{T}").c_str() );

    // Plot Unnormalized PbPt
    TGraphErrors* Pos_grUnnormPbPt = new TGraphErrors( Pos_Runs.size(), Pos_Runs.data(), Pos_UnnormPbPt.data(), nullptr, Pos_UnnormPbPtErr.data() );
    Pos_grUnnormPbPt->SetMarkerStyle( kFullCircle );
    Pos_grUnnormPbPt->SetMarkerColor( kViolet );
    Pos_mgUnnorm->Add( Pos_grUnnormPbPt, "p" );
    Pos_legUnnorm->AddEntry( Pos_grUnnormPbPt, "Unnorm. P_{B}P_{T}", "p" );

    Pos_mgUnnorm->SetTitle( string( "Positive Target "+Period_Title+" P_{B}P_{T} (Unnorm. to FC Charge); Run Number; P_{B}P_{T}").c_str() );

    // Plot the %-difference between normalized PbPt and unnormalized PbPt versus FC asymmetry
    TGraphErrors* Pos_deltaFC = new TGraphErrors( Pos_A_FC.size(), Pos_A_FC.data(), Pos_Delta_PbPt.data(), nullptr, Pos_Delta_PbPtErr.data() );
    Pos_deltaFC->SetMarkerStyle( kFullCircle );
    Pos_deltaFC->SetMarkerColor( kRed );
    r = Pos_deltaFC->Fit("pol1","S","Q");
    stringstream Pos_en1;  Pos_en1 << "[0] = "<< setprecision(3) << r->Value(0) <<" #pm "<< r->Error(0);
    stringstream Pos_en2;  Pos_en2 << "[1] = "<< setprecision(3) << r->Value(1) <<" #pm "<< r->Error(1);
    stringstream Pos_chi2; Pos_chi2<< "#chi^{2}_{red} = "<< setprecision(3) << r->Chi2() / r->Ndf();
    Pos_legDelta->AddEntry( (TObject*)0, Pos_en1.str().c_str(),  "p" );
    Pos_legDelta->AddEntry( (TObject*)0, Pos_en2.str().c_str(),  "p" );
    Pos_legDelta->AddEntry( (TObject*)0, Pos_chi2.str().c_str(), "p" );
    Pos_deltaFC->SetTitle( string( "Positive Target "+Period_Title+" Percent Diff. in P_{B}P_{T} versus A_{FC}; A_{FC}; Percent Difference").c_str() );

    // Create the lines of averages
    vector<double> Pos_AvgUnnorm = WeightedAverage( Pos_UnnormPbPt, Pos_UnnormPbPtErr );
    TF1* Pos_UnnormLine = new TF1("Pos_UnnormLine",to_string( Pos_AvgUnnorm[0] ).c_str(), Runs.front(), Runs.back() );
    Pos_UnnormLine->SetLineColor( kViolet );
    vector<double> Pos_AvgNorm = WeightedAverage( Pos_PbPt, Pos_PbPtErr );
    TF1* Pos_NormLine = new TF1("Pos_NormLine",to_string( Pos_AvgNorm[0] ).c_str(), Runs.front(), Runs.back() );
    Pos_NormLine->SetLineColor( kBlue );

    stringstream Pos_AvgUnVal; Pos_AvgUnVal <<"Average Unnorm. P_{B}P_{T} = "<< setprecision(3) << Pos_AvgUnnorm[0] <<" #pm "<< Pos_AvgUnnorm[1];
    Pos_legUnnorm->SetHeader( Pos_AvgUnVal.str().c_str(), "C" );
    stringstream Pos_AvgNVal; Pos_AvgNVal <<"Average Norm. P_{B}P_{T} = "<< setprecision(3) << Pos_AvgNorm[0] <<" #pm "<< Pos_AvgNorm[1];
    Pos_legNorm->SetHeader( Pos_AvgNVal.str().c_str(), "C" );
  

    // Negative plots
    TLegend* Neg_legUnnorm    = new TLegend(0.18,0.7,0.7,0.9); 
    Neg_legUnnorm->SetTextFont(43);
    Neg_legUnnorm->SetTextSize(15);
    Neg_legUnnorm->SetBorderSize(0);
    Neg_legUnnorm->SetFillStyle(0);
    TLegend* Neg_legNorm      = new TLegend(0.18,0.7,0.7,0.9); 
    Neg_legNorm->SetTextFont(43);
    Neg_legNorm->SetTextSize(15);
    Neg_legNorm->SetBorderSize(0);
    Neg_legNorm->SetFillStyle(0);
    TLegend* Neg_legDelta     = new TLegend(0.18,0.7,0.7,0.9);
    Neg_legDelta->SetTextFont(43);
    Neg_legDelta->SetTextSize(15);
    Neg_legDelta->SetBorderSize(0);
    Neg_legDelta->SetFillStyle(0);
    Neg_legDelta->SetHeader("Fit = [0] + [1]*A_{FC}","C");

    TMultiGraph* Neg_mgUnnorm = new TMultiGraph();
    TMultiGraph* Neg_mgNorm   = new TMultiGraph();

    // Plot the offline NMR values
    TGraphErrors* Neg_grNMR = new TGraphErrors( Neg_Runs.size(), Neg_Runs.data(), Neg_NMR_PbPt.data(), nullptr, Neg_NMR_PbPtErr.data() );
    Neg_grNMR->SetMarkerStyle( kFullSquare );
    Neg_grNMR->SetMarkerColor( kBlack );
    Neg_mgNorm->Add( Neg_grNMR, "p" );
    Neg_mgUnnorm->Add( Neg_grNMR, "p" );
    Neg_legNorm->AddEntry( Neg_grNMR, "Offline NMR P_{B}P_{T}", "p" );
    Neg_legUnnorm->AddEntry( Neg_grNMR, "Offline NMR P_{B}P_{T}", "p" );

    // Plot Normalized PbPt
    TGraphErrors* Neg_grNormPbPt = new TGraphErrors( Neg_Runs.size(), Neg_Runs.data(), Neg_PbPt.data(), nullptr, Neg_PbPtErr.data() );
    Neg_grNormPbPt->SetMarkerStyle( kFullCircle );
    Neg_grNormPbPt->SetMarkerColor( kBlue );
    Neg_mgNorm->Add( Neg_grNormPbPt, "p" );
    Neg_legNorm->AddEntry( Neg_grNormPbPt, "FC-Norm. P_{B}P_{T}", "p" );

    Neg_mgNorm->SetTitle( string( "Negative Target "+Period_Title+" P_{B}P_{T} (Norm. to FC Charge); Run Number; P_{B}P_{T}").c_str() );

    // Plot Unnormalized PbPt
    TGraphErrors* Neg_grUnnormPbPt = new TGraphErrors( Neg_Runs.size(), Neg_Runs.data(), Neg_UnnormPbPt.data(), nullptr, Neg_UnnormPbPtErr.data() );
    Neg_grUnnormPbPt->SetMarkerStyle( kFullCircle );
    Neg_grUnnormPbPt->SetMarkerColor( kViolet );
    Neg_mgUnnorm->Add( Neg_grUnnormPbPt, "p" );
    Neg_legUnnorm->AddEntry( Neg_grUnnormPbPt, "Unnorm. P_{B}P_{T}", "p" );
    Neg_mgUnnorm->SetTitle( string( "Negative Target "+Period_Title+" P_{B}P_{T} (Unnorm. to FC Charge); Run Number; P_{B}P_{T}").c_str() );

    // Plot the %-difference between normalized PbPt and unnormalized PbPt versus FC asymmetry
    TGraphErrors* Neg_deltaFC = new TGraphErrors( Neg_A_FC.size(), Neg_A_FC.data(), Neg_Delta_PbPt.data(), nullptr, Neg_Delta_PbPtErr.data() );
    Neg_deltaFC->SetMarkerStyle( kFullCircle );
    Neg_deltaFC->SetMarkerColor( kRed );
    r = Neg_deltaFC->Fit("pol1","S","Q");
    stringstream Neg_en1;  Neg_en1 << "[0] = "<< setprecision(3) << r->Value(0) <<" #pm "<< r->Error(0);
    stringstream Neg_en2;  Neg_en2 << "[1] = "<< setprecision(3) << r->Value(1) <<" #pm "<< r->Error(1);
    stringstream Neg_chi2; Neg_chi2<< "#chi^{2}_{red} = "<< setprecision(3) << r->Chi2() / r->Ndf();
    Neg_legDelta->AddEntry( "", Neg_en1.str().c_str(),  "p" );
    Neg_legDelta->AddEntry( "", Neg_en2.str().c_str(),  "p" );
    Neg_legDelta->AddEntry( "", Neg_chi2.str().c_str(), "p" );

    Neg_deltaFC->SetTitle( string( "Negative Target "+Period_Title+" Percent Diff. in P_{B}P_{T} versus A_{FC}; A_{FC}; Percent Difference").c_str() );

    // Create the lines of averages
    vector<double> Neg_AvgUnnorm = WeightedAverage( Neg_UnnormPbPt, Neg_UnnormPbPtErr );
    TF1* Neg_UnnormLine = new TF1("Neg_UnnormLine",to_string( Neg_AvgUnnorm[0] ).c_str(), Runs.front(), Runs.back() );
    Neg_UnnormLine->SetLineColor( kViolet );
    vector<double> Neg_AvgNorm = WeightedAverage( Neg_PbPt, Neg_PbPtErr );
    TF1* Neg_NormLine = new TF1("Neg_NormLine",to_string( Neg_AvgNorm[0] ).c_str(), Runs.front(), Runs.back() );
    Neg_NormLine->SetLineColor( kBlue );

    stringstream Neg_AvgUnVal; Neg_AvgUnVal <<"Average Unnorm. P_{B}P_{T} = "<< setprecision(3) << Neg_AvgUnnorm[0] <<" #pm "<< Neg_AvgUnnorm[1];
    Neg_legUnnorm->SetHeader( Neg_AvgUnVal.str().c_str(), "C" );
    stringstream Neg_AvgNVal; Neg_AvgNVal <<"Average Norm. P_{B}P_{T} = "<< setprecision(3) << Neg_AvgNorm[0] <<" #pm "<< Neg_AvgNorm[1];
    Neg_legNorm->SetHeader( Neg_AvgNVal.str().c_str(), "C" );


    c->cd(1); 
    gPad->SetLeftMargin(0.18); gPad->SetRightMargin(0.02);
    // Center the titles to match the margin
    double marginLeft   = gPad->GetLeftMargin();
    double marginRight  = gPad->GetRightMargin();
    double marginCenter = marginLeft + (1.0 - marginLeft - marginRight ) / 2.0;
    gStyle->SetTitleX( marginCenter ); gStyle->SetTitleAlign(23);
    Pos_mgNorm->GetXaxis()->SetNdivisions(505,true);
    Pos_mgNorm->GetXaxis()->SetRangeUser(runMin, runMax);
    Pos_mgNorm->GetYaxis()->SetRangeUser(0,1.2);
    Pos_mgNorm->Draw("ap"); Pos_legNorm->Draw("same"); Pos_NormLine->Draw("same");

    c->cd(2); 
    gPad->SetLeftMargin(0.18); gPad->SetRightMargin(0.02);
    gStyle->SetTitleX( marginCenter ); gStyle->SetTitleAlign(23);
    Pos_mgUnnorm->GetXaxis()->SetNdivisions(505,true);
    Pos_mgUnnorm->GetXaxis()->SetRangeUser(runMin, runMax);
    Pos_mgUnnorm->GetYaxis()->SetRangeUser(0,1.2);
    Pos_mgUnnorm->Draw("ap"); Pos_legUnnorm->Draw("same"); Pos_UnnormLine->Draw("same");

    c->cd(3); 
    gPad->SetLeftMargin(0.18); gPad->SetRightMargin(0.02);
    gStyle->SetTitleX( marginCenter ); gStyle->SetTitleAlign(23);
    Pos_deltaFC->GetXaxis()->SetLimits(-0.017,0.017);
    Pos_deltaFC->GetYaxis()->SetRangeUser(-0.5,0.5);
    Pos_deltaFC->Draw("ap"); Pos_legDelta->Draw("same");

    c->cd(4); 
    gPad->SetLeftMargin(0.18); gPad->SetRightMargin(0.02);
    gStyle->SetTitleX( marginCenter ); gStyle->SetTitleAlign(23);
    Neg_mgNorm->GetXaxis()->SetNdivisions(505,true);
    Neg_mgNorm->GetXaxis()->SetRangeUser(runMin, runMax);
    Neg_mgNorm->GetYaxis()->SetRangeUser(-1.2,0);
    Neg_mgNorm->Draw("ap"); Neg_legNorm->Draw("same"); Neg_NormLine->Draw("same");

    c->cd(5); 
    gPad->SetLeftMargin(0.18); gPad->SetRightMargin(0.02);
    gStyle->SetTitleX( marginCenter ); gStyle->SetTitleAlign(23);
    Neg_mgUnnorm->GetXaxis()->SetNdivisions(505,true);
    Neg_mgUnnorm->GetXaxis()->SetRangeUser(runMin, runMax);
    Neg_mgUnnorm->GetYaxis()->SetRangeUser(-1.2,0);
    Neg_mgUnnorm->Draw("ap"); Neg_legUnnorm->Draw("same"); Neg_UnnormLine->Draw("same");

    c->cd(6); 
    gPad->SetLeftMargin(0.18); gPad->SetRightMargin(0.02);
    gStyle->SetTitleX( marginCenter ); gStyle->SetTitleAlign(23);
    Neg_deltaFC->GetXaxis()->SetLimits(-0.017,0.017);
    Neg_deltaFC->GetYaxis()->SetRangeUser(-0.5,0.5);
    Neg_deltaFC->Draw("ap"); Neg_legDelta->Draw("same");

    return c;

}

// Holds boundaries for elastic Q2 bin values
vector<double> El_Q2_Bin_Bounds = {1.8, 1.923765, 2.056040, 2.197410, 2.348501, 2.509980, 2.682562, 2.867011, 3.064142, 3.274828, 3.5};


// This dumb struct holds the elastic scattering information for a given run
struct ElasticInfo{

    double A_el      = 0; // Theoretical elastic asymmetry from Noemie's code
    double N_P       = 0; // Number of helicity+1 scattering events
    double N_M       = 0; // Number of helicity-1 scattering events
    double FC_P      = 0; // Faraday Cup charge for helicity+1 data
    double FC_M      = 0; // Faraday Cup charge for helicity-1 data
    double Q2_Mean   = 0; // Mean value of Q2 for this bin
    double Q2_Low    = 0; // Lower boundary of Q2 bin (Q2 must be >= this value)
    double Q2_High   = 0; // Upper boundary of Q2 bin (Q2 must be < this value)
    double DF_el     = 0; // Elastic dilution factor
    double DF_el_err = 0; // Error in the elastic dilution factor

    void SetValues( string data ){
	stringstream sin(data); double theta_mean;
	sin >> Q2_Mean >> theta_mean >> A_el >> N_P >> N_M >> FC_P >> FC_M;
	int size = El_Q2_Bin_Bounds.size();
	for(int i=0; i<size-1; i++){
	    double qlow = El_Q2_Bin_Bounds[i]; double qhigh = El_Q2_Bin_Bounds[i+1];
	    if( Q2_Mean >= qlow && Q2_Mean < qhigh ){
		Q2_Low  = qlow;
		Q2_High = qhigh;
		return;
	    }
	}
        cout <<"ERROR: Couldn't find bin bounds for Q2_Mean = "<< Q2_Mean <<". Bin bounds not set!!\n";
        cout << data << endl;
    }

    // Checks if a given Q2 midpoint is within this boundary
    bool isInBin( double Q2Mid ) const{ return Q2Mid >= Q2_Low && Q2Mid < Q2_High; }

    // Sets DF info
    void SetDF( double df, double dferr ){
	DF_el = df;
	DF_el_err = dferr;
    }
    
    void Print() const{
	cout <<"****** For "<< Q2_Low <<" <= Q2 < "<< Q2_High <<" ******\n";
	cout <<"----> A_el    = "<< A_el << endl;
	cout <<"----> N_P     = "<< N_P << endl;
	cout <<"----> N_M     = "<< N_M << endl;
	cout <<"----> FC_P    = "<< FC_P << endl;
	cout <<"----> FC_M    = "<< FC_M << endl;
	cout <<"----> Q2_Mean = "<< Q2_Mean << endl;
	cout <<"----> DF_el   = "<< DF_el <<" +- "<< DF_el_err << endl;
	cout <<"***********************************\n";
    }

};

// This dumb struct just holds the PbPt information
struct PbPt{

    int RunNumber        = 0; // Run number
    double PbPt_Raw      = 0; // PbPt unnormalized to FC charge
    double PbPt_Raw_Err  = 0;
    double PbPt_Norm     = 0; // PbPt normalized to the FC charge
    double PbPt_Norm_Err = 0;
    // Holds the elastic info for convenience; vector for each elastic Q2 bin
    vector<ElasticInfo> Elastic_Bins;

    void SetValues( string data, string target="NH3" ){
	stringstream sin(data);
	sin >> RunNumber >> PbPt_Raw >> PbPt_Raw_Err >> PbPt_Norm >> PbPt_Norm_Err;
	// Once the run number is obtained, get all the elastic run info
	ifstream fin( string("../Latest_Skims/Elastic_Text_Files/"+ target +"_"+to_string(RunNumber) +"_Elastic.txt") );
	if( fin.fail() ){ cout <<"ERROR: Couldn't set elastic info in PbPt struct for run "<< RunNumber << endl; return; }

	string line;
	getline( fin, line ); // Throw away header row
	while( getline( fin, line ) ){
	    ElasticInfo EI;
	    EI.SetValues( line );
	    Elastic_Bins.push_back( EI );
	}
	fin.close();

    }

    // Based on the string query, return the given elastic bin info for the Q2 bin
    double GetElasticBinValue( string item, double q2Mid ) const{
	for(auto EI : Elastic_Bins){
	   if( EI.isInBin( q2Mid ) ){
		if( item == "A_el" ) return EI.A_el;
		else if( item == "N_P" ) return EI.N_P;
		else if( item == "N_M" ) return EI.N_M;
		else if( item == "FC_P" ) return EI.FC_P;
		else if( item == "FC_M" ) return EI.FC_M;
		else if( item == "Q2_Mean" ) return EI.Q2_Mean;
	   } 
	}
	cout << "ERROR: Unable to find query for "<< item <<" at Q2 midpoint "<< q2Mid <<" in run "<< RunNumber <<". Returning zero...\n";
	return 0;
    }

    // Sets the elastic DF based on the run number
    void SetElasticDF(){
	string Period, Tpol; // run period and target polarization
	// Set run period
	if( RunNumber < 16800 ) Period = "Su22";
	else if( RunNumber > 16800 && RunNumber < 17185 ) Period = "Fa22Neg";
	else if( RunNumber > 17185 && RunNumber < 17450 ) Period = "Fa22Pos";
	else if( RunNumber > 17450 && RunNumber <=17768 ) Period = "Sp23Inb";
	else Period = "Sp23Outb";
	// Set target polarization
	if( PbPt_Raw > 0 ) Tpol = "Pos";
	else if( PbPt_Raw < 0 ) Tpol = "Neg"; // Allows for Tpol to be unititialized for debugging...

	ifstream fin( string("Output_Data/Elastic_DF_Data_"+ Tpol +"_"+ Period +".txt") );
	if( fin.fail() ){ cout <<"ERROR: Unable to find elastic DF data for "<< Period <<" "<< Tpol <<" at run "<< RunNumber <<". Elastic DF not set!\n"; return; }

	// Bin_Q2 Bin_X DF_NH3 Err_DF_NH3	
	string line; 
	getline( fin, line ); // Throw away header row
	while( getline( fin, line ) ){
	    stringstream sin(line);
	    double q2mean, dummyX, df, dferr;
	    sin >> q2mean >> dummyX >> df >> dferr;
	    int size = Elastic_Bins.size();
	    for( int i=0; i<size; i++ ){
		if( Elastic_Bins[i].isInBin( q2mean ) ) Elastic_Bins[i].SetDF( df, dferr );
	    }
	}
	fin.close();

    }

    void Print() const{
	cout <<"==== Run "<< RunNumber <<" PbPt Info ====\n";
	cout <<"--> PbPt_Raw  = "<< PbPt_Raw <<" +- "<< PbPt_Raw_Err << endl;
	cout <<"--> PbPt_Norm = "<< PbPt_Norm <<" +- "<< PbPt_Norm_Err << endl;
	for(auto EI : Elastic_Bins ) EI.Print(); // Print elastic bin contents too
	cout <<"=============================\n\n";
    }

    // First entry is the numerator term for PNF, the second is the denominator term. Returns zeros if denominator is invalid
    vector<double> GetTermsPNF( string useNorm="Raw" ) const{
	double PNF_num = 0; double PNF_den = 0;
	for(auto EI : Elastic_Bins){
	    if( EI.Q2_Mean != 0 ){
		double PbPt = PbPt_Raw;
		if( useNorm != "Raw" ) PbPt = PbPt_Norm;
		PNF_num += EI.DF_el * (EI.N_P - EI.N_M) * EI.A_el * PbPt;
		PNF_den += pow( EI.DF_el, 2 ) * (EI.N_P + EI.N_M) * pow( (EI.A_el * PbPt ), 2 );
	    }
	}
	vector<double> PNF_terms = { PNF_num, PNF_den };
	return PNF_terms;
    }

};
// Returns a vector of all PbPt info
vector<PbPt> All_PbPt_Info(){
    
    vector<PbPt> Values;

    ifstream fin("Output_Data/All_NH3_PbPt.txt");
    if( fin.fail() ){ cout <<"ERROR: Couldn't open PbPt file. No PbPt info set! Check inputs.\n"; return Values; }

    string line;
    while( getline( fin, line ) ){
	PbPt thisRun;
	thisRun.SetValues( line ); // Sets the values for PbPt for this run specifically
	thisRun.SetElasticDF(); // Sets the elastic DF data
	Values.push_back( thisRun );
    }
    fin.close();

    return Values;
}

// This function calculates the PNF for a given run period and target polarization. Returns a vector with the PNF info.
// The first entry is the actual PNF value, the second is the statistical error on it
vector<double> Calculate_PNF( vector<PbPt>& AllPbPt, RunPeriod& RP, string Period, string Target, int targetPol ){

    vector<double> PNF = {0,0};

    // Get the run info
    auto Runs = RP.getElasticEpoch( Period, Target, targetPol );
    //for(int r : Runs) cout << r << endl;

    // PNF terms
    double PNF_num = 0; double PNF_den = 0;

    // Loop to calculate the PNF
    // Get a vector of all PbPt to use...
    for( int r : Runs ){
	// Find the appropriate value of PbPt
	for( PbPt p : AllPbPt ){
	    if( p.RunNumber == r ){
		auto PNF_terms = p.GetTermsPNF(); // Defalts to using raw counts
		PNF_num += PNF_terms[0];
		PNF_den += PNF_terms[1];
	    }
	}
    }

    if( PNF_den > 0 ){
	PNF[0] = PNF_num / PNF_den;
	PNF[1] = 1.0 / sqrt( PNF_den );
	cout <<"PNF for "<< Period <<", "<< Target <<", Pt ~ "<< targetPol <<": PNF = "<< PNF[0] <<" +- "<< PNF[1] << endl;
    }
    else cout <<"ERROR: Invalid calculation of PNF for "<< Period <<". Check inputs.\n";

    return PNF;

}

void Calculate_PbPt(){

    // If there's already a file holding the PbPt info, make sure it's deleted
    if( remove("Output_Data/All_NH3_PbPt.txt") == 0 ) cout <<"Removed NH3 PbPt file...\n";

    RunPeriod Period; // Holds the run info

    auto Plot_Su22    = Plot_PbPt_Epoch( Period, 1, 10, "Summer" );            //Plot_Su22->Print("PDF_Plots/PbPt_Vs_Runs.pdf(");
    auto Plot_Fa22Neg = Plot_PbPt_Epoch( Period, 11, 15, "Fall (Neg. Sol.)" ); //Plot_Fa22Neg->Print("PDF_Plots/PbPt_Vs_Runs.pdf");
    //auto Plot_Fa22Pos = Plot_PbPt_Epoch( Period, 16, 19, "Fall (Pos. Sol.)" ); Plot_Fa22Pos->Print("PDF_Plots/PbPt_Vs_Runs.pdf");
    auto Plot_Sp23Inb = Plot_PbPt_Epoch( Period, 20, 23, "Spring" );           //Plot_Sp23Inb->Print("PDF_Plots/PbPt_Vs_Runs.pdf)");

    // Now calculate the normalization factors

    // Read in the just-calculated DIS PbPt values
    auto DIS_PbPt = All_PbPt_Info();
    //for( auto r : DIS_PbPt ) r.Print();

    auto PNFs = Calculate_PNF( DIS_PbPt, Period, "Su22", "NH3", 1 );
    PNFs = Calculate_PNF( DIS_PbPt, Period, "Su22", "NH3",-1 );

    PNFs = Calculate_PNF( DIS_PbPt, Period, "Fa22Neg", "NH3", 1 );
    PNFs = Calculate_PNF( DIS_PbPt, Period, "Fa22Neg", "NH3",-1 );

    PNFs = Calculate_PNF( DIS_PbPt, Period, "Sp23Inb", "NH3", 1 );
    PNFs = Calculate_PNF( DIS_PbPt, Period, "Sp23Inb", "NH3",-1 );


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
