/***********************************************************************************
 *
 * Author: Derek Holmberg
 *
 * Date Created: 9/25/2026
 *
 * Last Modified: 10/8/2026
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

int globIt = 0; // Used in TObject names to prevent duplicate names

// This function appends a string object to a file of a given path
void AppendToFile( string fileName, string data ){

    ofstream fout;
    fout.open( fileName, ios_base::app ); // Set the file to append
    fout << data;

}

// This function takes a mean and a sigma and returns the value by the given sigma
double GaussianSmear( double val, double val_err ){

    TRandom3 rand;
    rand.SetSeed(0);

    return rand.Gaus( val, val_err );
    //double variance = rand.Gaus( val, val_err );
    //return variance;

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
TCanvas* Plot_PbPt_Epoch( RunPeriod& Period, int FirstEpoch, int LastEpoch, string Period_Title, vector<double> PNF_Pos={1.0,0}, vector<double> PNF_Neg={1.0,0} ){

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

		// Apply the correction for the PNF if applicable
		double PNF     = PNF_Pos[0];
	 	double PNF_Err = PNF_Pos[1];
		unnormPbPt = unnormPbPt * PNF;
	        unnormPbPtErr = sqrt( pow((PNF*unnormPbPtErr),2) + pow((unnormPbPt*PNF_Err),2) );

		normPbPt = normPbPt * PNF;
	        normPbPtErr = sqrt( pow((PNF*normPbPtErr),2) + pow((normPbPt*PNF_Err),2) );

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

		// Apply the correction for the PNF
		double PNF     = PNF_Neg[0];
	 	double PNF_Err = PNF_Neg[1];
		unnormPbPt = unnormPbPt * PNF;
	        unnormPbPtErr = sqrt( pow((PNF*unnormPbPtErr),2) + pow((unnormPbPt*PNF_Err),2) );

		normPbPt = normPbPt * PNF;
	        normPbPtErr = sqrt( pow((PNF*normPbPtErr),2) + pow((normPbPt*PNF_Err),2) );

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
	    if( PNF_Pos[1] != 0 || PNF_Neg[1] != 0 ) AppendToFile( "Output_Data/All_NH3_PbPt_PNF.txt", dataStream );
	    else AppendToFile( "Output_Data/All_NH3_PbPt.txt", dataStream );

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

    string canName = Period_Title+"_c"+to_string( globIt );
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
    TF1* Neg_UnnormLine = new TF1( string("Neg_UnnormLine"+to_string(globIt)).c_str(),to_string( Neg_AvgUnnorm[0] ).c_str(), Runs.front(), Runs.back() );
    Neg_UnnormLine->SetLineColor( kViolet );
    vector<double> Neg_AvgNorm = WeightedAverage( Neg_PbPt, Neg_PbPtErr );
    TF1* Neg_NormLine = new TF1( string("Neg_NormLine"+to_string(globIt)).c_str(),to_string( Neg_AvgNorm[0] ).c_str(), Runs.front(), Runs.back() );
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

    // Reset the gPad margins
    gPad->SetLeftMargin(0); gPad->SetRightMargin(0);

    globIt++;

    return c;

}

// This function is used in the next script. It checks whether or not a given run is within the
// specified run period.
bool IsInPeriod(int Run, string Period){

    if( Period == "Summer" )                return Run < 16800;
    else if( Period == "Fall (Neg. Sol.)" ) return Run > 16800 && Run < 17185;
    else if( Period == "Fall (Pos. Sol.)" ) return Run > 17185 && Run < 17450;
    else if( Period == "Spring" )           return Run > 17450 && Run <=17768;
    return false;

}

// This function plots the PbPt for a given epoch both with and without the PNF correction factor
// applied. Plots include both counts normalized to FC charge and unnormalized to FC charge. The PNF
// were calculated without normalization to FC charge.
TCanvas* Plot_PNF_PbPt( RunPeriod& RP, string Period, string Target ){

    cout << "Get the files\n";
    // Check that the files exist
    ifstream fin( string("Output_Data/All_"+Target+"_PbPt.txt") );
    ifstream finPNF( string("Output_Data/All_"+Target+"_PbPt_PNF.txt") );

    if( fin.fail() || finPNF.fail() ){
	cout <<"ERROR: At least one input file is missing. Returning empty canvas for "<< Period << endl;
	TCanvas* err = new TCanvas();
	return err;
    }

    // Read in the files and put them in vectors
    string line;
    vector<double> Runs, P_Unnorm, P_Unnorm_Err, P_Norm, P_Norm_Err;
    vector<double> RunsPNF, P_UnnormPNF, P_Unnorm_ErrPNF, P_NormPNF, P_Norm_ErrPNF;

    cout << "Read in from the files\n";
    // Read in the PbPt without PNF applied
    while( getline( fin, line ) ){ // No header row to remove for now...
	stringstream sin(line);
	double run, p_unnorm, p_unnormerr, p_norm, p_normerr;
	sin >> run >> p_unnorm >> p_unnormerr >> p_norm >> p_normerr;
	// Only add the runs if it's in the appropriate run period
	if( IsInPeriod( run, Period ) ){
	    Runs.push_back( run ); P_Unnorm.push_back( p_unnorm ); P_Unnorm_Err.push_back( p_unnormerr );
	    P_Norm.push_back( p_norm ); P_Norm_Err.push_back( p_normerr );
	}
    }
    fin.close();

    // Read in the PbPt with PNF applied
    while( getline( finPNF, line ) ){ // No header row to remove for now...
	stringstream sin(line);
	double run, p_unnorm, p_unnormerr, p_norm, p_normerr;
	sin >> run >> p_unnorm >> p_unnormerr >> p_norm >> p_normerr;
	// Only add the runs if it's in the appropriate run period
	if( IsInPeriod( run, Period ) ){
	    RunsPNF.push_back( run ); P_UnnormPNF.push_back( p_unnorm ); P_Unnorm_ErrPNF.push_back( p_unnormerr );
	    P_NormPNF.push_back( p_norm ); P_Norm_ErrPNF.push_back( p_normerr );
	}
    }
    finPNF.close();

    // Don't plot anything if there's a discrepancy in the number of runs between each variation
    if( (Runs.size() != RunsPNF.size()) || Runs.size() == 0 ){
	cout <<"ERROR: Difference in runs between PNF normalized and unnormalized runs. Returning empty canvas for "<< Period << endl;
	TCanvas* err = new TCanvas();
	return err;
    }

    // Now move on to create the plots
    cout <<"Start the plots\n";
    TMultiGraph* mg_Unnorm = new TMultiGraph();
    mg_Unnorm->SetTitle( string( Period+" Unnormalized P_{B}P_{T}; Run Number; PbPt" ).c_str() );
    TLegend* leg_Unnorm = new TLegend(0.1, 0.4, 0.7, 0.6);
    leg_Unnorm->SetTextFont(43);
    leg_Unnorm->SetTextSize(15);
    leg_Unnorm->SetBorderSize(0);
    leg_Unnorm->SetFillStyle(0);

    TGraphErrors* grUnnorm = new TGraphErrors( Runs.size(), Runs.data(), P_Unnorm.data(), nullptr, P_Unnorm_Err.data() );
    grUnnorm->SetMarkerStyle( kFullCircle );
    grUnnorm->SetMarkerColor( kRed );
    mg_Unnorm->Add( grUnnorm, "p" );
    leg_Unnorm->AddEntry( grUnnorm, "Unnorm. P_{B}P_{T}", "p" );

    TGraphErrors* grUnnormPNF = new TGraphErrors( RunsPNF.size(), RunsPNF.data(), P_UnnormPNF.data(), nullptr, P_Unnorm_ErrPNF.data() );
    grUnnormPNF->SetMarkerStyle( kFullCircle );
    grUnnormPNF->SetMarkerColor( kBlue );
    mg_Unnorm->Add( grUnnormPNF, "p" );
    leg_Unnorm->AddEntry( grUnnormPNF, "PNF Unnorm. P_{B}P_{T}", "p" );

    TMultiGraph* mg_Norm = new TMultiGraph();
    mg_Norm->SetTitle( string( Period+" P_{B}P_{T} Normalized to FC Charge; Run Number; PbPt" ).c_str() );
    TLegend* leg_Norm = new TLegend(0.1, 0.4, 0.7, 0.6);
    leg_Norm->SetTextFont(43);
    leg_Norm->SetTextSize(15);
    leg_Norm->SetBorderSize(0);
    leg_Norm->SetFillStyle(0);

    TGraphErrors* grNorm = new TGraphErrors( Runs.size(), Runs.data(), P_Norm.data(), nullptr, P_Norm_Err.data() );
    grNorm->SetMarkerStyle( kFullCircle );
    grNorm->SetMarkerColor( kRed );
    mg_Norm->Add( grNorm, "p" );
    leg_Norm->AddEntry( grNorm, "Norm. P_{B}P_{T}", "p" );

    TGraphErrors* grNormPNF = new TGraphErrors( RunsPNF.size(), RunsPNF.data(), P_NormPNF.data(), nullptr, P_Norm_ErrPNF.data() );
    grNormPNF->SetMarkerStyle( kFullCircle );
    grNormPNF->SetMarkerColor( kBlue );
    mg_Norm->Add( grNormPNF, "p" );
    leg_Norm->AddEntry( grNormPNF, "PNF Norm. P_{B}P_{T}", "p" );

    mg_Unnorm->GetYaxis()->SetRangeUser(-1,1);
    mg_Norm->GetYaxis()->SetRangeUser(-1,1);

    TCanvas* c = new TCanvas("c","c",1600,600); c->Divide(2,1);
    c->cd(1);

    mg_Unnorm->Draw("ap");
    leg_Unnorm->Draw("same");
    
    c->cd(2);

    mg_Norm->Draw("ap");
    leg_Norm->Draw("same");

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
    vector<double> GetTermsPNF( string useNorm="Raw", bool useGaussSmear=false ) const{
	double PNF_num = 0; double PNF_den = 0;
	for(auto EI : Elastic_Bins){
	    if( EI.Q2_Mean != 0 ){
		double PbPt = PbPt_Raw;
		if( useNorm != "Raw" ) PbPt = PbPt_Norm;
		// Apply a sign correction if it's in the positive solenoid period
		double corr_factor = 1.0;
		double thisDF_el = EI.DF_el;

		// Vary the dilution factor by a gaussian smear if a systematic study is being done
		if( useGaussSmear ) thisDF_el = GaussianSmear( EI.DF_el, EI.DF_el_err );

		if( RunNumber > 17185 && RunNumber < 17450 ) corr_factor = -1.0;
		PNF_num += thisDF_el * (EI.N_P - EI.N_M) * EI.A_el * PbPt * corr_factor;
		PNF_den += pow( thisDF_el, 2 ) * (EI.N_P + EI.N_M) * pow( (EI.A_el * PbPt ), 2 );
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
vector<double> Calculate_PNF( vector<PbPt>& AllPbPt, RunPeriod& RP, string Period, string Target, int targetPol, bool useGaussSmear=false ){

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
		auto PNF_terms = p.GetTermsPNF("Raw",useGaussSmear); // Defalts to using raw counts
		PNF_num += PNF_terms[0];
		PNF_den += PNF_terms[1];
	    }
	}
    }

    if( PNF_den != 0 ){
	PNF[0] = PNF_num / PNF_den;
	PNF[1] = 1.0 / sqrt( PNF_den );
	if( !useGaussSmear ) cout <<"PNF for "<< Period <<", "<< Target <<", Pt ~ "<< targetPol <<": PNF = "<< PNF[0] <<" +- "<< PNF[1] << endl;
    }

    else cout <<"ERROR: Invalid calculation of PNF for "<< Period <<". Check inputs.\n";

    return PNF;

}

// This function calculates the variance in the PNF for the elastic DFs for all run periods and target polarization states.
// The vector returns the TCanvas with all the histograms and writes the final values to the "Variations" vector.
TCanvas* PNF_ElasticDF_Stat_Err( vector<PbPt>& AllPbPt, RunPeriod& RP, string Period, string Target, int numVariations, vector<double>& Variations ){

    // Create histograms for each run period. There are three histograms for each run period, one for pos Tpol, neg Tpol, and the weighted average
    TH1D* Pos_PNF = new TH1D( string(Period+"_"+Target+"_Pos").c_str(), "", 2000, 0.8, 1.2);
    Pos_PNF->SetStats(false);
    Pos_PNF->SetLineColor( kRed );
    Pos_PNF->SetTitle( string( Period+" PNF Distributions for Variance in d_{el}; PNF;" ).c_str() );
    TH1D* Neg_PNF = new TH1D( string(Period+"_"+Target+"_Neg").c_str(), "", 2000, 0.8, 1.2);
    Neg_PNF->SetStats(false);
    Neg_PNF->SetLineColor( kBlue );
    Neg_PNF->SetTitle( string( Period+" PNF Distributions for Variance in d_{el}; PNF;" ).c_str() );
    TH1D* Avg_PNF = new TH1D( string(Period+"_"+Target+"_Avg").c_str(), "", 2000, 0.8, 1.2);
    Avg_PNF->SetStats(false);
    Avg_PNF->SetLineColor( kViolet );
    Avg_PNF->SetTitle( string( Period+" PNF Distributions for Variance in d_{el}; PNF;" ).c_str() );

    // Start the data loop
    for(int i=0; i<numVariations; i++){

	if( i%1000 == 0 ) cout <<"Processed "<< i <<" variations for "<< Period <<", "<< Target <<"...\n";
	
	// Positive target data
	vector<double> thisVarPos = Calculate_PNF( AllPbPt, RP, Period, Target, 1, true );
	Pos_PNF->Fill( thisVarPos[0] );

	// Negitive target data
	vector<double> thisVarNeg = Calculate_PNF( AllPbPt, RP, Period, Target,-1, true );
	Neg_PNF->Fill( thisVarNeg[0] );

	// Average target data
	vector<double> allPNFs    = { thisVarPos[0], thisVarNeg[0] };
	vector<double> allPNFErrs = { thisVarPos[1], thisVarNeg[1] };
	vector<double> thisVarAvg = WeightedAverage( allPNFs, allPNFErrs );
	Avg_PNF->Fill( thisVarAvg[0] );

    }

    // Plot everything and set the values
    TLegend* leg = new TLegend(0.65,0.7,0.9,0.9);
    
    // Make the strings for the legend
    stringstream legAvgPNF; legAvgPNF <<"Avg. PNF = "<< setprecision(3) << Avg_PNF->GetMean() <<" #pm "<< Avg_PNF->GetStdDev();
    leg->AddEntry( Avg_PNF, legAvgPNF.str().c_str(), "lp" );

    stringstream legPosPNF; legPosPNF <<"P_{t} > 0 PNF = "<< setprecision(3) << Pos_PNF->GetMean() <<" #pm "<< Pos_PNF->GetStdDev();
    leg->AddEntry( Pos_PNF, legPosPNF.str().c_str(), "lp" );

    stringstream legNegPNF; legNegPNF <<"P_{t} < 0 PNF = "<< setprecision(3) << Neg_PNF->GetMean() <<" #pm "<< Neg_PNF->GetStdDev();
    leg->AddEntry( Neg_PNF, legNegPNF.str().c_str(), "lp" );

    // Set the errors
    Variations = { Pos_PNF->GetStdDev(), Neg_PNF->GetStdDev(), Avg_PNF->GetStdDev() };

    // Adjust the ranges based on the positions of the peaks
    double xLow  = Avg_PNF->GetMean() - 100*Avg_PNF->GetStdDev();
    double xHigh = Avg_PNF->GetMean() + 100*Avg_PNF->GetStdDev();
    // Find the largest bin size to adjust the y-axis
    int yMax = 0;
    for(int i=1; i<=2000; i++) if( Avg_PNF->GetBinContent(i) > yMax ) yMax = Avg_PNF->GetBinContent(i);

    Pos_PNF->GetXaxis()->SetRangeUser( xLow, xHigh ); Pos_PNF->GetYaxis()->SetRangeUser( 0, yMax*1.1 );
    Neg_PNF->GetXaxis()->SetRangeUser( xLow, xHigh ); Neg_PNF->GetYaxis()->SetRangeUser( 0, yMax*1.1 );
    Avg_PNF->GetXaxis()->SetRangeUser( xLow, xHigh ); Avg_PNF->GetYaxis()->SetRangeUser( 0, yMax*1.1 );

    gStyle->SetTitleX(0.5);
    gStyle->SetTitleAlign(23);

    string canName = "PNF"+Period+Target;
    TCanvas* c = new TCanvas( canName.c_str(), canName.c_str(), 800, 600 );
    c->cd();
    Pos_PNF->Draw("colz"); Neg_PNF->Draw("same"); Avg_PNF->Draw("same"); leg->Draw("same");
    
    return c;

}

// This function takes vectors of all the positive and negative PNFs for each run period and makes a plot out of them,
// as well as calculating the weighted average for each run period
TCanvas* Plot_All_PNFs( vector<double>& All_Pos_PNF, vector<double>& All_Pos_PNF_Err, vector<double>& All_Neg_PNF, vector<double>& All_Neg_PNF_Err ){

    // Kill the program if the sizes aren't correct
    if( All_Pos_PNF.size() != 4 or All_Pos_PNF_Err.size() != 4 or All_Neg_PNF.size() != 4 or All_Neg_PNF_Err.size() != 4 ){
	cout <<"ERROR: Not enough PNFs were passed in. Returning empty canvas...\n";
	TCanvas* err = new TCanvas(); return err;
    }

    // Make the vectors for each run period
    vector<double> Pos_Periods = {1, 2, 3, 4};
    vector<double> Neg_Periods = {1.1, 2.1, 3.1, 4.1};
    vector<double> Avg_Periods = {1.2, 2.2, 3.2, 4.2};

    // Calculate the weighted averages for each of the run periods
    vector<double> Su22_Vals    = {All_Pos_PNF[0], All_Neg_PNF[0]}; vector<double> Su22_Errs    = {All_Pos_PNF_Err[0], All_Neg_PNF_Err[0]};
    vector<double> Fa22Neg_Vals = {All_Pos_PNF[1], All_Neg_PNF[1]}; vector<double> Fa22Neg_Errs = {All_Pos_PNF_Err[1], All_Neg_PNF_Err[1]};
    vector<double> Fa22Pos_Vals = {All_Pos_PNF[2], All_Neg_PNF[2]}; vector<double> Fa22Pos_Errs = {All_Pos_PNF_Err[2], All_Neg_PNF_Err[2]};
    vector<double> Sp23Inb_Vals = {All_Pos_PNF[3], All_Neg_PNF[3]}; vector<double> Sp23Inb_Errs = {All_Pos_PNF_Err[3], All_Neg_PNF_Err[3]};

    vector<double> Su22_Avg_PNF    = WeightedAverage( Su22_Vals, Su22_Errs );
    vector<double> Fa22Neg_Avg_PNF = WeightedAverage( Fa22Neg_Vals, Fa22Neg_Errs );
    vector<double> Fa22Pos_Avg_PNF = WeightedAverage( Fa22Pos_Vals, Fa22Pos_Errs );
    vector<double> Sp23Inb_Avg_PNF = WeightedAverage( Sp23Inb_Vals, Sp23Inb_Errs );

//    vector<double> Fa22Neg_Avg_PNF = WeightedAverage( {All_Pos_PNF[1], All_Neg_PNF[1]}, {All_Pos_PNF_Err[1], All_Neg_PNF_Err[1]} );
//    vector<double> Fa22Pos_Avg_PNF = WeightedAverage( {All_Pos_PNF[2], All_Neg_PNF[2]}, {All_Pos_PNF_Err[2], All_Neg_PNF_Err[2]} );
//    vector<double> Sp23Inb_Avg_PNF = WeightedAverage( {All_Pos_PNF[3], All_Neg_PNF[3]}, {All_Pos_PNF_Err[3], All_Neg_PNF_Err[3]} );

    // Put the averages in a vector so they can be plotted
    vector<double> All_Avg_PNF     = { Su22_Avg_PNF[0], Fa22Neg_Avg_PNF[0], Fa22Pos_Avg_PNF[0], Sp23Inb_Avg_PNF[0] };    
    vector<double> All_Avg_PNF_Err = { Su22_Avg_PNF[1], Fa22Neg_Avg_PNF[1], Fa22Pos_Avg_PNF[1], Sp23Inb_Avg_PNF[1] };    

    // Make a graph and legend to hold the information
    TMultiGraph* mg = new TMultiGraph();
    TLegend* leg = new TLegend(0.65, 0.7, 0.9, 0.9);

    TGraphErrors* grPos = new TGraphErrors( All_Pos_PNF.size(), Pos_Periods.data(), All_Pos_PNF.data(), nullptr, All_Pos_PNF_Err.data() );
    grPos->SetMarkerStyle( kFullCircle );
    grPos->SetMarkerColor( kRed );
    mg->Add( grPos, "p" );
    leg->AddEntry( grPos, "PNF for P_{t} > 0", "p" );

    TGraphErrors* grNeg = new TGraphErrors( All_Neg_PNF.size(), Neg_Periods.data(), All_Neg_PNF.data(), nullptr, All_Neg_PNF_Err.data() );
    grNeg->SetMarkerStyle( kFullCircle );
    grNeg->SetMarkerColor( kBlue );
    mg->Add( grNeg, "p" );
    leg->AddEntry( grNeg, "PNF for P_{t} < 0", "p" );

    TGraphErrors* grAvg = new TGraphErrors( All_Avg_PNF.size(), Avg_Periods.data(), All_Avg_PNF.data(), nullptr, All_Avg_PNF_Err.data() );
    grAvg->SetMarkerStyle( kFullCircle );
    grAvg->SetMarkerColor( kViolet );
    mg->Add( grAvg, "p" );
    leg->AddEntry( grAvg, "Avg. PNF All P_{t}", "p" );

    mg->SetTitle("All PNF vs. Run Period; Run Period; PNF");

    TCanvas* Plt = new TCanvas("PNF","PNF",800,600);
    Plt->cd();
    mg->Draw("ap");
    leg->Draw("same");

    return Plt;
}

void Calculate_PbPt(){


    // If there's already a file holding the PbPt info, make sure it's deleted
    if( remove("Output_Data/All_NH3_PbPt.txt")     == 0 ) cout <<"Removed NH3 PbPt file...\n";
    if( remove("Output_Data/All_NH3_PbPt_PNF.txt") == 0 ) cout <<"Removed NH3 PNF-corrected PbPt file...\n";

    RunPeriod Period; // Holds the run info

    // Plots the DIS PbPt for each of the epochs
    auto Plot_Su22    = Plot_PbPt_Epoch( Period, 1, 10, "Summer" );            Plot_Su22->Print("PDF_Plots/PbPt_Vs_Runs.pdf(");
    auto Plot_Fa22Neg = Plot_PbPt_Epoch( Period, 11, 15, "Fall (Neg. Sol.)" ); Plot_Fa22Neg->Print("PDF_Plots/PbPt_Vs_Runs.pdf");
    auto Plot_Fa22Pos = Plot_PbPt_Epoch( Period, 16, 19, "Fall (Pos. Sol.)" ); Plot_Fa22Pos->Print("PDF_Plots/PbPt_Vs_Runs.pdf");
    auto Plot_Sp23Inb = Plot_PbPt_Epoch( Period, 20, 23, "Spring" );           Plot_Sp23Inb->Print("PDF_Plots/PbPt_Vs_Runs.pdf)");

    // Delete canvases to free memory
    Plot_Su22->Close();    delete Plot_Su22;
    Plot_Fa22Neg->Close(); delete Plot_Fa22Neg;
    Plot_Fa22Pos->Close(); delete Plot_Fa22Pos;
    Plot_Sp23Inb->Close(); delete Plot_Sp23Inb;

    // Used for visualizing the elastic DF error propagation
    TCanvas* Elastic_DF_Errors = new TCanvas("ElDF","ElDF",1600,1200); Elastic_DF_Errors->Divide(2,2);
    vector<double> Variations; // This vector holds the elastic DF errors
			       // The first entry 

    // Now calculate the normalization factors

    // Read in the just-calculated DIS PbPt values
    auto DIS_PbPt = All_PbPt_Info();
    //for( auto r : DIS_PbPt ) r.Print();

    // Su22 data
    auto Su22_Pos_PNFs = Calculate_PNF( DIS_PbPt, Period, "Su22", "NH3", 1 );
    auto Su22_Neg_PNFs = Calculate_PNF( DIS_PbPt, Period, "Su22", "NH3",-1 );
        // Plot the elastic DF errors
        auto PNF_Su22_Plot = PNF_ElasticDF_Stat_Err( DIS_PbPt, Period, "Su22", "NH3", 1000, Variations );
        Elastic_DF_Errors->cd(1);
        PNF_Su22_Plot->DrawClonePad();
    // Add the elastic DF errors in quadrature with the other statistical errors on PNF
    Su22_Pos_PNFs[1] = sqrt( pow( Su22_Pos_PNFs[1], 2) + pow( Variations[0], 2 ) );
    Su22_Neg_PNFs[1] = sqrt( pow( Su22_Neg_PNFs[1], 2) + pow( Variations[1], 2 ) );

    auto Plot_Su22_PNF = Plot_PbPt_Epoch( Period, 1, 10, "Summer", Su22_Pos_PNFs, Su22_Neg_PNFs );
    Plot_Su22_PNF->Print("PDF_Plots/PNF_PbPt_Vs_Runs.pdf(");


    // Fa22Neg data
    auto Fa22Neg_Pos_PNFs = Calculate_PNF( DIS_PbPt, Period, "Fa22Neg", "NH3", 1 );
    auto Fa22Neg_Neg_PNFs = Calculate_PNF( DIS_PbPt, Period, "Fa22Neg", "NH3",-1 );
        // Plot the elastic DF errors
        auto PNF_Fa22Neg_Plot = PNF_ElasticDF_Stat_Err( DIS_PbPt, Period, "Fa22Neg", "NH3", 1000, Variations );
        Elastic_DF_Errors->cd(2);
        PNF_Fa22Neg_Plot->DrawClonePad();
    // Add the elastic DF errors in quadrature with the other statistical errors on PNF
    Fa22Neg_Pos_PNFs[1] = sqrt( pow( Fa22Neg_Pos_PNFs[1], 2) + pow( Variations[0], 2 ) );
    Fa22Neg_Neg_PNFs[1] = sqrt( pow( Fa22Neg_Neg_PNFs[1], 2) + pow( Variations[1], 2 ) );

    auto Plot_Fa22Neg_PNF = Plot_PbPt_Epoch( Period, 11, 15, "Fall (Neg. Sol.)", Fa22Neg_Pos_PNFs, Fa22Neg_Neg_PNFs );
    Plot_Fa22Neg_PNF->Print("PDF_Plots/PNF_PbPt_Vs_Runs.pdf");


    // Fa22Pos data
    auto Fa22Pos_Pos_PNFs = Calculate_PNF( DIS_PbPt, Period, "Fa22Pos", "NH3", 1 );
    auto Fa22Pos_Neg_PNFs = Calculate_PNF( DIS_PbPt, Period, "Fa22Pos", "NH3",-1 );
        // Plot the elastic DF errors
        auto PNF_Fa22Pos_Plot = PNF_ElasticDF_Stat_Err( DIS_PbPt, Period, "Fa22Pos", "NH3", 1000, Variations );
        Elastic_DF_Errors->cd(3);
        PNF_Fa22Pos_Plot->DrawClonePad();
    // Add the elastic DF errors in quadrature with the other statistical errors on PNF
    Fa22Pos_Pos_PNFs[1] = sqrt( pow( Fa22Pos_Pos_PNFs[1], 2) + pow( Variations[0], 2 ) );
    Fa22Pos_Neg_PNFs[1] = sqrt( pow( Fa22Pos_Neg_PNFs[1], 2) + pow( Variations[1], 2 ) );

    auto Plot_Fa22Pos_PNF = Plot_PbPt_Epoch( Period, 16, 19, "Fall (Pos. Sol.)", Fa22Pos_Pos_PNFs, Fa22Pos_Neg_PNFs );
    Plot_Fa22Pos_PNF->Print("PDF_Plots/PNF_PbPt_Vs_Runs.pdf");


    // Sp23Inb data
    auto Sp23Inb_Pos_PNFs = Calculate_PNF( DIS_PbPt, Period, "Sp23Inb", "NH3", 1 );
    auto Sp23Inb_Neg_PNFs = Calculate_PNF( DIS_PbPt, Period, "Sp23Inb", "NH3",-1 );
        // Plot the elastic DF errors
        auto PNF_Sp23Inb_Plot = PNF_ElasticDF_Stat_Err( DIS_PbPt, Period, "Sp23Inb", "NH3", 1000, Variations );
        Elastic_DF_Errors->cd(4);
        PNF_Sp23Inb_Plot->DrawClonePad();
    // Add the elastic DF errors in quadrature with the other statistical errors on PNF
    Sp23Inb_Pos_PNFs[1] = sqrt( pow( Sp23Inb_Pos_PNFs[1], 2) + pow( Variations[0], 2 ) );
    Sp23Inb_Neg_PNFs[1] = sqrt( pow( Sp23Inb_Neg_PNFs[1], 2) + pow( Variations[1], 2 ) );

    auto Plot_Sp23Inb_PNF = Plot_PbPt_Epoch( Period, 20, 23, "Spring", Sp23Inb_Pos_PNFs, Sp23Inb_Neg_PNFs );
    Plot_Sp23Inb_PNF->Print("PDF_Plots/PNF_PbPt_Vs_Runs.pdf)");

    // Plot and draw the comparisons with and without the PNF applied
    auto Su22_PNF_Plt    = Plot_PNF_PbPt( Period, "Summer", "NH3" );           Su22_PNF_Plt->Print("PDF_Plots/PNF_Comparison.pdf(");
    auto Fa22Neg_PNF_Plt = Plot_PNF_PbPt( Period, "Fall (Neg. Sol.)", "NH3" ); Fa22Neg_PNF_Plt->Print("PDF_Plots/PNF_Comparison.pdf");
    auto Fa22Pos_PNF_Plt = Plot_PNF_PbPt( Period, "Fall (Pos. Sol.)", "NH3" ); Fa22Pos_PNF_Plt->Print("PDF_Plots/PNF_Comparison.pdf");
    auto Sp23_PNF_Plt    = Plot_PNF_PbPt( Period, "Spring", "NH3" );           Sp23_PNF_Plt->Print("PDF_Plots/PNF_Comparison.pdf)");


    // This plots the PNFs for each run period
    vector<double> All_Pos_PNF     = { Su22_Pos_PNFs[0], Fa22Neg_Pos_PNFs[0], Fa22Pos_Pos_PNFs[0], Sp23Inb_Pos_PNFs[0] };
    vector<double> All_Pos_PNF_Err = { Su22_Pos_PNFs[1], Fa22Neg_Pos_PNFs[1], Fa22Pos_Pos_PNFs[1], Sp23Inb_Pos_PNFs[1] };
    vector<double> All_Neg_PNF     = { Su22_Neg_PNFs[0], Fa22Neg_Neg_PNFs[0], Fa22Pos_Neg_PNFs[0], Sp23Inb_Neg_PNFs[0] };
    vector<double> All_Neg_PNF_Err = { Su22_Neg_PNFs[1], Fa22Neg_Neg_PNFs[1], Fa22Pos_Neg_PNFs[1], Sp23Inb_Neg_PNFs[1] };

    auto All_PNF_Plt = Plot_All_PNFs( All_Pos_PNF, All_Pos_PNF_Err, All_Neg_PNF, All_Neg_PNF_Err );
    All_PNF_Plt->Print("PDF_Plots/Weighted_Avg_PNF.pdf");

   
    /*------------------------------------------------- Print all the PNF values with the full errrors ----------------------------------------------------*/

    cout <<"Summer PNFs:\n";
    cout <<"---> Su22 Pos PNF = "<< Su22_Pos_PNFs[0] <<" +- "<< Su22_Pos_PNFs[1] << endl;
    cout <<"---> Su22 Neg PNF = "<< Su22_Neg_PNFs[0] <<" +- "<< Su22_Neg_PNFs[1] << endl;
    cout << endl;

    cout <<"Fall (Neg. Sol.) PNFs:\n";
    cout <<"---> Fa22Neg Pos PNF = "<< Fa22Neg_Pos_PNFs[0] <<" +- "<< Fa22Neg_Pos_PNFs[1] << endl;
    cout <<"---> Fa22Neg Neg PNF = "<< Fa22Neg_Neg_PNFs[0] <<" +- "<< Fa22Neg_Neg_PNFs[1] << endl;
    cout << endl;

    cout <<"Fall (Pos. Sol.) PNFs:\n";
    cout <<"---> Fa22Pos Pos PNF = "<< Fa22Pos_Pos_PNFs[0] <<" +- "<< Fa22Pos_Pos_PNFs[1] << endl;
    cout <<"---> Fa22Pos Neg PNF = "<< Fa22Pos_Neg_PNFs[0] <<" +- "<< Fa22Pos_Neg_PNFs[1] << endl;
    cout << endl;

    cout <<"Spring (Inb.) PNFs:\n";
    cout <<"---> Sp23Inb Pos PNF = "<< Sp23Inb_Pos_PNFs[0] <<" +- "<< Sp23Inb_Pos_PNFs[1] << endl;
    cout <<"---> Sp23Inb Neg PNF = "<< Sp23Inb_Neg_PNFs[0] <<" +- "<< Sp23Inb_Neg_PNFs[1] << endl;
    cout << endl;


    // Calculate the error propagation for the elastic dilution factors
//TCanvas* PNF_ElasticDF_Stat_Err( vector<PbPt>& AllPbPt, RunPeriod& RP, string Period, string Target, int numVariations, vector<double>& Variations ){
/*
    TCanvas* Elastic_DF_Errors = new TCanvas("ElDF","ElDF",1600,1200); Elastic_DF_Errors->Divide(2,2);

    vector<double> Variations;

    auto PNF_Su22_Plot = PNF_ElasticDF_Stat_Err( DIS_PbPt, Period, "Su22", "NH3", 1000, Variations );
    Elastic_DF_Errors->cd(1);
    PNF_Su22_Plot->DrawClonePad();
    cout <<"Su22 errors:\n";
    for(auto v : Variations ) cout << v <<" "<< endl;

    auto PNF_Fa22Neg_Plot = PNF_ElasticDF_Stat_Err( DIS_PbPt, Period, "Fa22Neg", "NH3", 1000, Variations );
    Elastic_DF_Errors->cd(2);
    PNF_Fa22Neg_Plot->DrawClonePad();
    cout <<"Fa22Neg errors:\n";
    for(auto v : Variations ) cout << v <<" "<< endl;

    auto PNF_Fa22Pos_Plot = PNF_ElasticDF_Stat_Err( DIS_PbPt, Period, "Fa22Pos", "NH3", 1000, Variations );
    Elastic_DF_Errors->cd(3);
    PNF_Fa22Pos_Plot->DrawClonePad();
    cout <<"Fa22Pos errors:\n";
    for(auto v : Variations ) cout << v <<" "<< endl;

    auto PNF_Sp23Inb_Plot = PNF_ElasticDF_Stat_Err( DIS_PbPt, Period, "Sp23Inb", "NH3", 1000, Variations );
    Elastic_DF_Errors->cd(4);
    PNF_Sp23Inb_Plot->DrawClonePad();
    cout <<"Sp23Inb errors:\n";
    for(auto v : Variations ) cout << v <<" "<< endl;

    Elastic_DF_Errors->Print("PDF_Plots/PNF_Elastic_DF_Error_Propagation.pdf");
*/
//TCanvas* Plot_PNF_PbPt( RunPeriod& RP, string Period, string Target ){

}
