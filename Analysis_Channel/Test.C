/***********************************************************************************
 *
 * Author: Derek Holmberg
 *
 * Date Created: 9/18/2026
 *
 * Last Modified: 9/25/2026
 *
 * Purpose:
 * This program is used to calculate the dilution factors (DF), packing fractions
 * (PF), and beam-target polarizations (PbPt) for the NH3 data. How I want to 
 * structure things is still up in the air, so I'm keeping everything in "Test" for
 * now...
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

// This function returns a TCanvas object holding all the dilution factors/packing fractions
// for a given set of epochs.
// "RP" is a reference to the RunPeriod object that holds the run info
// "FirstEpoch" is the start of the epoch to look at
// "LastEpoch" is the last epoch number
// "Target" is the target type (NH3 or ND3)
// "Period" specifies the summer, fall, spring, etc.
// "DForPF" specifies whether to plot DF or PF values
// "BathOrCell" tells the program whether to plot values for the Bath or Cell PFs; defaults to
// bath if unspecified by the user
// "writeDFPFtoFile" controls whether to write both the bath and cell PF data to a default file

TCanvas* RunPeriodPlots( RunPeriod& RP, int FirstEpoch, int LastEpoch, string Target, string Period, string DForPF, string BathOrCell="", bool writeDFPFtoFile=false ){

    // Determine the most appropriate division of the canvas
    int Rows = 1; int Cols = 1;
    int nPlots = LastEpoch - FirstEpoch + 1;
    while( Rows*Cols < nPlots ){
	Cols++;
	if( Rows*Cols < nPlots ) Rows++;
    }

    // Used for the canvas size
    int width = 400*Cols; int height = 300*Rows;

    // Create the canvas
    string canName = Period +"_"+ DForPF + BathOrCell +"_Ep_"+ to_string(FirstEpoch) +"-"+ to_string(LastEpoch);
    TCanvas* c = new TCanvas( canName.c_str(), canName.c_str(), width, height );
    c->Divide( Cols, Rows );

    // Start the loop over the epochs
    int it = 1; // Iterator used for plotting
    for(int ep=FirstEpoch; ep<=LastEpoch; ep++){

	string epID; // Used to grab the run numbers for a given epoch
	if( ep < 10 ) epID = "P0"+to_string(ep);
	else epID = "P"+to_string(ep);

	vector<int> Epoch = RP.getRunEpoch( epID );
	vector<int> MissedRuns;
	
	// Now start the loop for this particular epoch:
	DataSet thisEpoch;

	thisEpoch.SetEpochValues( Period, Target, Epoch, MissedRuns ); // Read in the runs
	thisEpoch.CalculateAvgXQ2( Period, Epoch, Target ); // Set the kinematic values for every bin
	thisEpoch.CalculateDF( Period, true, true ); // Calculate DF and PF
	string outTXT = "Output_Data/All_DF_Data_Epoch_"+to_string(ep);
	thisEpoch.WriteToCSV( outTXT, Period, " " ); // Writes data to a text or CSV file

	//thisEpoch.Print();

	if( MissedRuns.size() > 0 ){ // If runs were missed, abort the plot and return what's there
	    cout <<"ERROR: Missing runs for Epoch "<< ep <<"! Runs missing are:\n";
	    for(int r : MissedRuns) cout <<"--> "<< r << endl;
	    cout << endl;
	    return c;
	}

	// If the user wants PF information, this vector holds it
	vector<double> PF_Info = {0,0,0,0,0,0};
	// If the user wants to make DF data against which PF is plotted for
	// each epoch, this vector holds it
	vector<double> DF_Info = {0,0};

	string epName = "Epoch "+to_string(ep);

	if( DForPF == "DF" ){
	    c->cd(it);
	    DF_Legend_Plot( thisEpoch, epName, Target, DF_Info, gPad );
	    c->Modified();
	    c->Update();
	    // Write the DF info
	    stringstream DFout; DFout << ep <<"	";
	    for(auto val : DF_Info) DFout << val <<"	";
	    DFout << endl; 
	    if( writeDFPFtoFile ) AppendToFile( string("Output_Data/All_"+ Target +"_DF_Data.txt"), DFout.str() );

	}
	else if( DForPF == "PF" ){
	   
	    c->cd(it); 
	    PF_Legend_Plot( thisEpoch, epName, Target, PF_Info, gPad, BathOrCell );

	    c->Modified();
	    c->Update();
	    // Write the PF info
	    stringstream PFout; PFout << ep <<"	";
	    for(auto val : PF_Info) PFout << val <<"	";
	    PFout << endl;	    
	    if( writeDFPFtoFile ) AppendToFile( string("Output_Data/All_"+ Target +"_PF_Data.txt"), PFout.str() );

	}

	it++;

    }

    return c;

}

// This function makes a plot of the Q2-averaged DFs for the specified epochs in the RGC data set
TCanvas* PlotAverage( RunPeriod& RP, string Target, string DForPF, int FirstEpoch, int LastEpoch, string BathOrCell="" ){

    // Create the canvas
    string canName = "Average_"+ BathOrCell + DForPF +"s_All_Epochs_"+ to_string(FirstEpoch) +"-"+ to_string(LastEpoch);
    TCanvas* c = new TCanvas( canName.c_str(), canName.c_str(), 2200, 600 );
    c->Divide( 3, 1 );

    // Start the loop over the epochs, starting with the summer data
    vector<DataSet> Summer, Fall, Spring; // Holds the data for each epoch
    vector<int> MissedRuns; // Tracks any runs that are missing

    // Summer data
    for(int ep=FirstEpoch; ep<=LastEpoch; ep++){

	string epID; // Used to grab the run numbers for a given epoch
	if( ep < 10 ) epID = "P0"+to_string(ep);
	else epID = "P"+to_string(ep);

	string Period;
	if( ep <= 10 ) Period = "Su22";
	else if( ep > 10 && ep <= 15 ) Period = "Fa22Neg";
	else if( ep > 15 && ep <= 19 ) Period = "Fa22Pos";
	else if( ep > 19 && ep <= 23 ) Period = "Sp23Inb";

	vector<int> Epoch = RP.getRunEpoch( epID );
	
	// Now start the loop for this particular epoch:
	DataSet thisEpoch;

	thisEpoch.SetEpochValues( Period, Target, Epoch, MissedRuns ); // Read in the runs
	thisEpoch.CalculateAvgXQ2( Period, Epoch, Target ); // Set the kinematic values for every bin
	thisEpoch.CalculateDF( Period, true, true ); // Calculate DF and PF
	if( ep <= 10 ) Summer.push_back( thisEpoch );
	else if( ep > 10 && ep <= 19 ) Fall.push_back( thisEpoch );
	else if( ep > 19 && ep <= 23 ) Spring.push_back( thisEpoch );

    }

    //TCanvas* DF_Q2_Averaged_Plot( vector<DataSet>& Epochs, string epoch, string targetType, TVirtualPad* pad, double ymin=0.1, double ymax=0.3 ){
    // Plot summer data
    // Plot the DF data
    if( DForPF == "DF" ){
        c->cd(1);
        DF_Q2_Averaged_Plot( Summer, "Summer", "NH3", gPad );
        c->cd(2);
        DF_Q2_Averaged_Plot( Fall, "Fall", "NH3", gPad );
        c->cd(3);
        DF_Q2_Averaged_Plot( Spring, "Spring", "NH3", gPad );
    }
    // Plot the PF data
    else if( DForPF == "PF" ){
        c->cd(1);
        PF_Q2_Averaged_Plot( Summer, "Summer", "NH3", BathOrCell, gPad );
        c->cd(2);
        PF_Q2_Averaged_Plot( Fall, "Fall", "NH3", BathOrCell, gPad );
        c->cd(3);
        PF_Q2_Averaged_Plot( Spring, "Spring", "NH3", BathOrCell, gPad );
    }

    if( MissedRuns.size() > 0 ){
	cout <<"ERROR: There are missing runs! No plots were made.\n";
	cout <<"       Missing runs:\n";
	for(int r : MissedRuns) cout <<"--> "<< r << endl;
	cout << endl;
	TCanvas* fail;
	return fail;
    }

    return c;

}

// This function creates a linearity test to check the behavior of the PF vs. the DF.
// The errors on the DF and PF are statistical only and are added in quadrature and 
// plotted as an error on the PF.
TCanvas* Linearity_DF_PF( string Target="NH3" ){

    string line;

    // Read in the DF data
    vector<double> DFepoch, DFs, DF_Errs, ScaleDFs, ScaleDF_Errs;

    ifstream DFin("Output_Data/All_"+ Target +"_DF_Data.txt");
    while( getline( DFin, line ) ){
	stringstream sin(line);
	double epoch, df, dferr;
	sin >> epoch >> df >> dferr;
	DFepoch.push_back( epoch );
	DFs.push_back( df );
	DF_Errs.push_back( dferr );
	ScaleDFs.push_back( 3*df );
	ScaleDF_Errs.push_back( 3*dferr );
    }

    // Read in the PF data
    vector<double> PF_epSu22, PF_Su22, PF_ErrsSu22;
    vector<double> PF_epFa22Neg, PF_Fa22Neg, PF_ErrsFa22Neg;
    vector<double> PF_epFa22Pos, PF_Fa22Pos, PF_ErrsFa22Pos;
    vector<double> PF_epSp23Inb, PF_Sp23Inb, PF_ErrsSp23Inb;
    vector<double> PFepoch, PFs, PF_Errs;

    ifstream PFin("Output_Data/All_"+ Target +"_PF_Data.txt");
    while( getline( PFin, line ) ){
	stringstream sin(line);
	double epoch, pf, pferr;
	sin >> epoch >> pf >> pferr;
	PFepoch.push_back( epoch );
	PFs.push_back( pf );
	PF_Errs.push_back( pferr );

	if( epoch <= 10 ){
	    PF_epSu22.push_back( epoch );
	    PF_Su22.push_back( pf );
	    PF_ErrsSu22.push_back( pferr );
	}
	else if( epoch > 10 && epoch <= 15 ){
	    PF_epFa22Neg.push_back( epoch );
	    PF_Fa22Neg.push_back( pf );
	    PF_ErrsFa22Neg.push_back( pferr );
	}
	else if( epoch > 15 && epoch <= 19 ){
	    PF_epFa22Pos.push_back( epoch );
	    PF_Fa22Pos.push_back( pf );
	    PF_ErrsFa22Pos.push_back( pferr );
	}
	else if( epoch > 19 && epoch <= 23 ){
	    PF_epSp23Inb.push_back( epoch );
	    PF_Sp23Inb.push_back( pf );
	    PF_ErrsSp23Inb.push_back( pferr );
	}

    }

    if( PFepoch.size() != DFepoch.size() ){
	cout <<"ERROR: DF and PF files don't cover the same epoch range. Check inputs in 'Output_Data/'\n";
	TCanvas* fail;
	return fail;
    }

    // Plot the PF and scaled DF together as a function of epoch
    TMultiGraph* mgEpoch = new TMultiGraph();
    TLegend* legEpoch = new TLegend(0.1,0.7,0.3,0.9);
    // Scaled DFs
    TGraphErrors* grSDF = new TGraphErrors( DFepoch.size(), DFepoch.data(), ScaleDFs.data(), nullptr, ScaleDF_Errs.data() );
    grSDF->SetMarkerStyle( kFullCircle );
    grSDF->SetMarkerColor( kBlack );
    mgEpoch->Add( grSDF, "p" );
    legEpoch->AddEntry( grSDF, "3*D_{F}", "p" );
    // PF info
    // Su22
    TGraphErrors* grPFSu22 = new TGraphErrors( PF_epSu22.size(), PF_epSu22.data(), PF_Su22.data(), nullptr, PF_ErrsSu22.data() );
    grPFSu22->SetMarkerStyle( kFullSquare );
    grPFSu22->SetMarkerColor( kBlue );
    mgEpoch->Add( grPFSu22, "p" );
    legEpoch->AddEntry( grPFSu22, "P_{F} Su22", "p" );

    // Fa22Neg
    TGraphErrors* grPFFa22Neg = new TGraphErrors( PF_epFa22Neg.size(), PF_epFa22Neg.data(), PF_Fa22Neg.data(), nullptr, PF_ErrsFa22Neg.data() );
    grPFFa22Neg->SetMarkerStyle( kFullSquare );
    grPFFa22Neg->SetMarkerColor( kGreen );
    mgEpoch->Add( grPFFa22Neg, "p" );
    legEpoch->AddEntry( grPFFa22Neg, "P_{F} Fa22Neg", "p" );

    // Fa22Pos
    TGraphErrors* grPFFa22Pos = new TGraphErrors( PF_epFa22Pos.size(), PF_epFa22Pos.data(), PF_Fa22Pos.data(), nullptr, PF_ErrsFa22Pos.data() );
    grPFFa22Pos->SetMarkerStyle( kFullSquare );
    grPFFa22Pos->SetMarkerColor( kViolet );
    mgEpoch->Add( grPFFa22Pos, "p" );
    legEpoch->AddEntry( grPFFa22Pos, "P_{F} Fa22Pos", "p" );

    // Sp23Inb
    TGraphErrors* grPFSp23Inb = new TGraphErrors( PF_epSp23Inb.size(), PF_epSp23Inb.data(), PF_Sp23Inb.data(), nullptr, PF_ErrsSp23Inb.data() );
    grPFSp23Inb->SetMarkerStyle( kFullSquare );
    grPFSp23Inb->SetMarkerColor( kRed );
    mgEpoch->Add( grPFSp23Inb, "p" );
    legEpoch->AddEntry( grPFSp23Inb, "P_{F} Sp23Inb", "p" );

    // Now make the plots to show the PF vs. unscaled DF for each inbending run period.
    // Create DF vectors for each run period, treating the error as the quadrature sum
    // of the DF and PF statistical errors.
    vector<double> DF_Su22, DF_Fa22Neg, DF_Fa22Pos, DF_Sp23Inb;
    vector<double> DF_Err_Su22, DF_Err_Fa22Neg, DF_Err_Fa22Pos, DF_Err_Sp23Inb;
    for( size_t i=0; i<PFepoch.size(); i++ ){
	double ep        = PFepoch[i];
	double thisPF    = PFs[i];
	double thisPFErr = PF_Errs[i];
	double totalErr  = 0;
	if( ep <= 10 ){
	    DF_Su22.push_back( DFs[i] );
	    totalErr = sqrt( pow(thisPFErr,2) + pow( (thisPF*DF_Errs[i]/DFs[i]),2) );
	    DF_Err_Su22.push_back( totalErr );
	}
	else if( ep > 10 && ep <= 15 ){
	    DF_Fa22Neg.push_back( DFs[i] );
	    totalErr = sqrt( pow(thisPFErr,2) + pow( (thisPF*DF_Errs[i]/DFs[i]),2) );
	    DF_Err_Fa22Neg.push_back( totalErr );
	}
	else if( ep > 15 && ep <= 19 ){
	    DF_Fa22Pos.push_back( DFs[i] );
	    totalErr = sqrt( pow(thisPFErr,2) + pow( (thisPF*DF_Errs[i]/DFs[i]),2) );
	    DF_Err_Fa22Pos.push_back( totalErr );
	}
	else if( ep > 19 && ep <= 23 ){
	    DF_Sp23Inb.push_back( DFs[i] );
	    totalErr = sqrt( pow(thisPFErr,2) + pow( (thisPF*DF_Errs[i]/DFs[i]),2) );
	    DF_Err_Sp23Inb.push_back( totalErr );
	}
    }
   
    // Now make the PF vs. DF plot
    TMultiGraph* mgPF = new TMultiGraph();
    TLegend* legPF = new TLegend(0.1,0.6,0.7,0.9);
    // Plot and fit each to a line
    TFitResultPtr r;
    // PF vs. DF Su22 plot
    TGraphErrors* grSu22 = new TGraphErrors( DF_Su22.size(), DF_Su22.data(), PF_Su22.data(), nullptr, DF_Err_Su22.data() );
    grSu22->SetMarkerStyle( kFullCircle );
    grSu22->SetMarkerColor( kBlue );
    r = grSu22->Fit("pol1","S","Q");
    grSu22->GetFunction("pol1")->SetLineColor( kBlue );
    stringstream sSu22; sSu22 << "Su22, A = "<< setprecision(2) << r->Value(1) <<" #pm "<< r->Error(1) <<
				  ", B = "<< setprecision(2) << r->Value(0) <<" #pm "<< r->Error(0) << endl;
    mgPF->Add( grSu22, "p" );
    legPF->AddEntry( grSu22, sSu22.str().c_str(), "p" );

    // PF vs. DF Fa22Neg plot
    TGraphErrors* grFa22Neg = new TGraphErrors( DF_Fa22Neg.size(), DF_Fa22Neg.data(), PF_Fa22Neg.data(), nullptr, DF_Err_Fa22Neg.data() );
    grFa22Neg->SetMarkerStyle( kFullCircle );
    grFa22Neg->SetMarkerColor( kGreen );
    r = grFa22Neg->Fit("pol1","S","Q");
    grFa22Neg->GetFunction("pol1")->SetLineColor( kGreen );
    stringstream sFa22Neg; sFa22Neg << "Fa22Neg, A = "<< setprecision(2) << r->Value(1) <<" #pm "<< r->Error(1) <<
				  ", B = "<< setprecision(2) << r->Value(0) <<" #pm "<< r->Error(0) << endl;
    mgPF->Add( grFa22Neg, "p" );
    legPF->AddEntry( grFa22Neg, sFa22Neg.str().c_str(), "p" );

    // PF vs. DF Fa22Pos plot
    TGraphErrors* grFa22Pos = new TGraphErrors( DF_Fa22Pos.size(), DF_Fa22Pos.data(), PF_Fa22Pos.data(), nullptr, DF_Err_Fa22Pos.data() );
    grFa22Pos->SetMarkerStyle( kFullCircle );
    grFa22Pos->SetMarkerColor( kViolet );
    r = grFa22Pos->Fit("pol1","S","Q");
    grFa22Pos->GetFunction("pol1")->SetLineColor( kViolet );
    stringstream sFa22Pos; sFa22Pos << "Fa22Pos, A = "<< setprecision(2) << r->Value(1) <<" #pm "<< r->Error(1) <<
				  ", B = "<< setprecision(2) << r->Value(0) <<" #pm "<< r->Error(0) << endl;
    mgPF->Add( grFa22Pos, "p" );
    legPF->AddEntry( grFa22Pos, sFa22Pos.str().c_str(), "p" );

    // PF vs. DF Sp23Inb plot
    TGraphErrors* grSp23Inb = new TGraphErrors( DF_Sp23Inb.size(), DF_Sp23Inb.data(), PF_Sp23Inb.data(), nullptr, DF_Err_Sp23Inb.data() );
    grSp23Inb->SetMarkerStyle( kFullCircle );
    grSp23Inb->SetMarkerColor( kRed );
    r = grSp23Inb->Fit("pol1","S","Q");
    grSp23Inb->GetFunction("pol1")->SetLineColor( kRed );
    stringstream sSp23Inb; sSp23Inb << "Sp23Inb, A = "<< setprecision(2) << r->Value(1) <<" #pm "<< r->Error(1) <<
				  ", B = "<< setprecision(2) << r->Value(0) <<" #pm "<< r->Error(0) << endl;
    mgPF->Add( grSp23Inb, "p" );
    legPF->AddEntry( grSp23Inb, sSp23Inb.str().c_str(), "p" );


    // Plot everything together
    TCanvas* c = new TCanvas("DF_linearity","DF_linearity",1600,600);
    c->Divide(2,1);

    c->cd(1);
    mgEpoch->SetTitle("D_{F}, P_{F} for (x #approx 0.21, Q^{2} #approx 2.89 GeV^{2}) vs. Epoch; Epoch; 3*D_{F}, P_{F}");
    mgEpoch->GetYaxis()->SetRangeUser(0.4,0.65);
    mgEpoch->Draw("ap");
    legEpoch->Draw("same");

    c->cd(2);
    mgPF->SetTitle( string("P_{F} vs. D_{F} for All "+ Target +" Epochs; D_{F}; P_{F}").c_str() );
    legPF->SetHeader("P_{F} = A*D_{F} + B","C");
    mgPF->GetYaxis()->SetRangeUser(0.4,0.65);
    mgPF->Draw("ap");
    legPF->Draw("same");

    return c;
}


void Test(){

    RunPeriod Period; // Holds the information on all RGC runs

    // Draw the summer DF data
    auto DF_Su22  = RunPeriodPlots( Period, 1, 10, "NH3", "Su22", "DF", "", true );
    auto PFc_Su22 = RunPeriodPlots( Period, 1, 10, "NH3", "Su22", "PF", "Cell", true );
    auto PFb_Su22 = RunPeriodPlots( Period, 1, 10, "NH3", "Su22", "PF", "Bath" );

    DF_Su22->Print("PDF_Plots/DF.pdf(");
    PFc_Su22->Print("PDF_Plots/PFc.pdf(");
    PFb_Su22->Print("PDF_Plots/PFb.pdf(");

    // Draw the summer DF data
    auto DF_Fa22Neg  = RunPeriodPlots( Period, 11, 15, "NH3", "Fa22Neg", "DF", "", true );
    auto PFc_Fa22Neg = RunPeriodPlots( Period, 11, 15, "NH3", "Fa22Neg", "PF", "Cell", true );
    auto PFb_Fa22Neg = RunPeriodPlots( Period, 11, 15, "NH3", "Fa22Neg", "PF", "Bath" );

    DF_Fa22Neg->Print("PDF_Plots/DF.pdf");
    PFc_Fa22Neg->Print("PDF_Plots/PFc.pdf");
    PFb_Fa22Neg->Print("PDF_Plots/PFb.pdf");


    // Draw the summer DF data
    auto DF_Fa22Pos  = RunPeriodPlots( Period, 16, 19, "NH3", "Fa22Pos", "DF", "", true );
    auto PFc_Fa22Pos = RunPeriodPlots( Period, 16, 19, "NH3", "Fa22Pos", "PF", "Cell", true );
    auto PFb_Fa22Pos = RunPeriodPlots( Period, 16, 19, "NH3", "Fa22Pos", "PF", "Bath" );

    DF_Fa22Pos->Print("PDF_Plots/DF.pdf");
    PFc_Fa22Pos->Print("PDF_Plots/PFc.pdf");
    PFb_Fa22Pos->Print("PDF_Plots/PFb.pdf");


    // Draw the summer DF data
    auto DF_Sp23Inb  = RunPeriodPlots( Period, 20, 23, "NH3", "Sp23Inb", "DF", "", true );
    auto PFc_Sp23Inb = RunPeriodPlots( Period, 20, 23, "NH3", "Sp23Inb", "PF", "Cell", true );
    auto PFb_Sp23Inb = RunPeriodPlots( Period, 20, 23, "NH3", "Sp23Inb", "PF", "Bath" );

    DF_Sp23Inb->Print("PDF_Plots/DF.pdf)");
    PFc_Sp23Inb->Print("PDF_Plots/PFc.pdf)");
    PFb_Sp23Inb->Print("PDF_Plots/PFb.pdf)");



    auto DF_Q2_Plt  = PlotAverage( Period, "NH3", "DF", 1, 23 );
    auto PFb_Q2_Plt = PlotAverage( Period, "NH3", "PF", 1, 23, "Bath" );
    auto PFc_Q2_Plt = PlotAverage( Period, "NH3", "PF", 1, 23, "Cell" );
    DF_Q2_Plt->Print( "PDF_Plots/Q2_Averaged_DF_and_PF_NH3.pdf(");
    PFb_Q2_Plt->Print("PDF_Plots/Q2_Averaged_DF_and_PF_NH3.pdf" );
    PFc_Q2_Plt->Print("PDF_Plots/Q2_Averaged_DF_and_PF_NH3.pdf)");

    auto DF_Linearity = Linearity_DF_PF();
    DF_Linearity->Print("PDF_Plots/DF_vs_PF_Plots.pdf");

//TCanvas* RunPeriodPlots( RunPeriod& RP, int FirstEpoch, int LastEpoch, string Target, string Period, string DForPF, string BathOrCell="Bath" ){

/*
    DataSet Test;

    cout << "Header file compiles.\n";

    RunPeriod Period;

    vector<int> Epoch2 = Period.getRunEpoch("P02");
    vector<int> MissedRuns;
    // Read in the data from epoch 2
    Test.SetEpochValues( "Su22", "NH3", Epoch2, MissedRuns );

    // Calculate the average x, Q2 bins for Epoch2..........
    Test.CalculateAvgXQ2( "Su22", Epoch2, "NH3" );
    Test.CalculateDF("Su22", true, true);

    // Now make the DF and PF plots for the epoch:
// TCanvas* PF_Legend_Plot( DataSet& AllData, string sector, string targetType, vector<double>& PFs, string Period="Default", bool useScaling =true, bool usePseudoData=true,
// double ymin=0.48, double ymax=0.6, bool dummy=false ){

// TCanvas* DF_Legend_Plot( DataSet& AllData, string sector, string targetType, string Period, bool useScaling=true, bool usePseudoData=true ){

    vector<double> PFs = {0,0,0,0,0,0};

    auto DF_Plt = DF_Legend_Plot( Test, "Epoch 2", "NH3", "Su22" );
    auto PF_Plt = PF_Legend_Plot( Test, "Epoch 2", "NH3", PFs );

    //Test.Print();

    if( MissedRuns.size() > 0 ){
	cout << "S*** F*****G D*RN IT!! Runs were missed! Missing runs:\n";
	for(int run : MissedRuns ) cout <<"--> "<< run << endl;
	cout << endl;
    }
*/

}
