/***********************************************************************************
 *
 * Author: Derek Holmberg
 *
 * Date Created: 9/18/2026
 *
 * Last Modified: 9/18/2026
 *
 * Purpose:
 * This program loops over the ROOT files from the output of the "Data_Loop.C" 
 * script and generates text files that serve as the main input loop to the spin
 * structure analysis. These text files are stored in "Latest_Skims/Text_Files/".
 *
***********************************************************************************/

#include "../Header_Files/Binning_Classes.h"

using namespace std;

// Object used to hold the plots for different data sets
struct DiagnosticPlots{

    DiagnosticPlots() = default;
    DiagnosticPlots( string nameExt ) {	
	Theta->SetName( string("Theta"+nameExt).c_str() );
 	Phi->SetName( string("Phi"+nameExt).c_str() );
	Nphe->SetName( string("Nphe"+nameExt).c_str() );
	Vx->SetName( string("Vx"+nameExt).c_str() );
	Vy->SetName( string("Vy"+nameExt).c_str() );
	Vz->SetName( string("Vz"+nameExt).c_str() );
	Q2->SetName( string("Q2"+nameExt).c_str() );
	W->SetName( string("W"+nameExt).c_str() );
	X->SetName( string("X"+nameExt).c_str() );
	M2_Dist->SetName( string("M2_Dist"+nameExt).c_str() );
    }

    // Various kinematic plots used for diagnostics
    TH1D* Theta = new TH1D("Theta", "Theta Distribution; #theta [deg];",450,0,45);
    TH1D* Phi = new TH1D("Phi", "Phi Distribution; #phi [deg];",3600,-180,180);
    TH1D* Nphe = new TH1D("nphe", "Photoelectrons; Nphe;;",50,0,50);
    TH1D* Vx = new TH1D("Vx", "X-Vertex Distribution; V_{X} [cm];",2000,-15,5);
    TH1D* Vy = new TH1D("Vy", "Y-Vertex Distribution; V_{Y} [cm];",2000,-15,5);
    TH1D* Vz = new TH1D("Vz", "Z-Vertex Distribution; V_{Z} [cm];",2000,-15,5);
    TH1D* Q2 = new TH1D("Q2", "Q^{2} Distribution; Q^{2} [GeV^{2}];",1100,0,11);
    TH1D* W = new TH1D("W", "W Distribution; W [GeV];",500,0,5);
    TH1D* X = new TH1D("X", "X Distribution; X;",100,0,1);
    TH1D* M2_Dist = new TH1D("M2_Dist", "ECAL Second Moments Distribution; M2 [cm^{2}];",4500,0,450);

};

// This function takes a set of DiagnosticPlots and draws them together
TCanvas* DrawDiagnosticPlots( vector<DiagnosticPlots>& Plots, string nameExt ){

    // Entry 1 is the summer data, followed by fallneg, fallpos, and then spring inbending...
    vector<int> palette = {1, 2, 3, 4}; // Sets the colors
    TCanvas* AllPlots = new TCanvas( string("AllPlots"+nameExt).c_str(), string("AllPlots"+nameExt).c_str(), 2400, 1200 );
    AllPlots->Divide(5,2);

    // Loop over each entry
    int index = 1;
    AllPlots->cd(index);
    for(size_t i=0; i<Plots.size(); i++){
	Plots[i].Theta->SetStats(false);
	Plots[i].Theta->SetLineColor( palette[i] ); Plots[i].Theta->Scale( 1.0 / Plots[i].Theta->Integral() );
	if( i == 0 ) Plots[i].Theta->Draw("hist");
	else Plots[i].Theta->Draw("hist same");
    } index++;
    AllPlots->cd(index);
    for(size_t i=0; i<Plots.size(); i++){
	Plots[i].Phi->SetStats(false);
	Plots[i].Phi->SetLineColor( palette[i] ); Plots[i].Phi->Scale( 1.0 / Plots[i].Phi->Integral() );
	if( i == 0 ) Plots[i].Phi->Draw("hist");
	else Plots[i].Phi->Draw("hist same");
    } index++;
    AllPlots->cd(index);
    for(size_t i=0; i<Plots.size(); i++){
	Plots[i].Vx->SetStats(false);
	Plots[i].Vx->GetXaxis()->SetRangeUser(-4,4);
	Plots[i].Vx->SetLineColor( palette[i] ); Plots[i].Vx->Scale( 1.0 / Plots[i].Vx->Integral() );
	if( i == 0 ) Plots[i].Vx->Draw("hist");
	else Plots[i].Vx->Draw("hist same");
    } index++;
    AllPlots->cd(index);
    for(size_t i=0; i<Plots.size(); i++){
	Plots[i].Vy->SetStats(false);
	Plots[i].Vy->GetXaxis()->SetRangeUser(-4,4);
	Plots[i].Vy->SetLineColor( palette[i] ); Plots[i].Vy->Scale( 1.0 / Plots[i].Vy->Integral() );
	if( i == 0 ) Plots[i].Vy->Draw("hist");
	else Plots[i].Vy->Draw("hist same");
    } index++;
    AllPlots->cd(index);
    for(size_t i=0; i<Plots.size(); i++){
	Plots[i].Vz->SetStats(false);
	Plots[i].Vz->GetXaxis()->SetRangeUser(-10,5);
	Plots[i].Vz->SetLineColor( palette[i] ); Plots[i].Vz->Scale( 1.0 / Plots[i].Vz->Integral() );
	if( i == 0 ) Plots[i].Vz->Draw("hist");
	else Plots[i].Vz->Draw("hist same");
    } index++;
    AllPlots->cd(index);
    for(size_t i=0; i<Plots.size(); i++){
	Plots[i].Q2->SetStats(false);
	Plots[i].Q2->SetLineColor( palette[i] ); Plots[i].Q2->Scale( 1.0 / Plots[i].Q2->Integral() );
	if( i == 0 ) Plots[i].Q2->Draw("hist");
	else Plots[i].Q2->Draw("hist same");
    } index++;
    AllPlots->cd(index);
    for(size_t i=0; i<Plots.size(); i++){
	Plots[i].W->SetStats(false);
	Plots[i].W->GetXaxis()->SetRangeUser(1,5);
	Plots[i].W->SetLineColor( palette[i] ); Plots[i].W->Scale( 1.0 / Plots[i].W->Integral() );
	if( i == 0 ) Plots[i].W->Draw("hist");
	else Plots[i].W->Draw("hist same");
    } index++;
    AllPlots->cd(index);
    for(size_t i=0; i<Plots.size(); i++){
	Plots[i].X->SetStats(false);
	Plots[i].X->SetLineColor( palette[i] ); Plots[i].X->Scale( 1.0 / Plots[i].X->Integral() );
	if( i == 0 ) Plots[i].X->Draw("hist");
	else Plots[i].X->Draw("hist same");
    } index++;
    AllPlots->cd(index);
    for(size_t i=0; i<Plots.size(); i++){
	Plots[i].M2_Dist->SetStats(false);
	Plots[i].M2_Dist->GetXaxis()->SetRangeUser(0,100);
	Plots[i].M2_Dist->SetLineColor( palette[i] ); Plots[i].M2_Dist->Scale( 1.0 / Plots[i].M2_Dist->Integral() );
	if( i == 0 ) Plots[i].M2_Dist->Draw("hist");
	else Plots[i].M2_Dist->Draw("hist same");
    } index++;
    AllPlots->cd(index);
    for(size_t i=0; i<Plots.size(); i++){
	Plots[i].Nphe->SetStats(false);
	Plots[i].Nphe->SetLineColor( palette[i] ); Plots[i].Nphe->Scale( 1.0 / Plots[i].Nphe->Integral() );
	if( i == 0 ) Plots[i].Nphe->Draw("hist");
	else Plots[i].Nphe->Draw("hist same");
    }

    // Make a legend
    TLegend* leg = new TLegend(0.6,0.5,0.9,0.9);
    leg->AddEntry( Plots[0].Phi, "Summer", "lp" );
    leg->AddEntry( Plots[1].Phi, "Fall (Neg. Sol.)", "lp" );
    leg->AddEntry( Plots[2].Phi, "Fall (Pos. Sol.)", "lp" );
    leg->AddEntry( Plots[3].Phi, "Spring (Inb.)", "lp" );
    AllPlots->cd(1);
    leg->Draw("same");

    return AllPlots;
}

// Dumb struct that holds the fcup information from the input file
struct Input{

    Input() = default;
    Input( int run, double fcp, double fcm, double rp, double rm, double fcp_nc, double fcm_nc, double rp_nc, double rm_nc )
	: Run(run), FC_P(fcp), FC_M(fcm), ReadP(rp), ReadM(rm), 
	  FC_P_NoCut(fcp_nc), FC_M_NoCut(fcm_nc), ReadP_NoCut(rp_nc), ReadM_NoCut(rm_nc) {}

    int Run            = 0; // Run number
    double FC_P        = 0; // Helicity+1 charge (nC)
    double FC_M        = 0; // Helicity-1 charge (nC)
    double ReadP       = 0; // Number of Helicity+1 bank readouts
    double ReadM       = 0; // Number of Helicity-1 bank readouts
    double FC_P_NoCut  = 0; // Helicity+1 w/o additional cuts (nC)
    double FC_M_NoCut  = 0; // Helicity-1 w/o additional cuts (nC)
    double ReadP_NoCut = 0; // Number of Helicity+1 bank readouts w/o cuts
    double ReadM_NoCut = 0; // Number of Helicity-1 bank readouts w/o cuts

    void Print() const{
	cout <<"========= Run "<< Run <<" =========\n";
	cout <<"--> FC_P        = "<< FC_P << endl;
	cout <<"--> FC_M        = "<< FC_M << endl;
	cout <<"--> ReadP       = "<< ReadP << endl;
	cout <<"--> ReadM       = "<< ReadM << endl;
	cout <<"--> FC_P_NoCut  = "<< FC_P_NoCut << endl;
	cout <<"--> FC_M_NoCut  = "<< FC_M_NoCut << endl;
	cout <<"--> ReadP_NoCut = "<< ReadP_NoCut << endl;
	cout <<"--> ReadM_NoCut = "<< ReadM_NoCut << endl;
	cout <<"=============================\n";
    }

};

// Creates a vector that holds all the FCup information...
vector<Input> Create_FCData( string fcupPath ){

    // Get the FC data ready to write
    vector<Input> FCData;
    
    ifstream fcin( fcupPath.c_str() );
    string line;
    getline( fcin, line ); // Throw away header row with column ID's
    while( getline( fcin, line ) ){
	stringstream sin(line);
	int run, run2; double fcp, fcm, readp, readm, fcp_nc, fcm_nc, readp_nc, readm_nc;
	sin >> run >> fcp >> fcm >> readp >> readm >> run2 >> fcp_nc >> fcm_nc >> readp_nc >> readm_nc;
	
	Input thisRun( run, fcp, fcm, readp, readm, fcp_nc, fcm_nc, readp_nc, readm_nc );
	FCData.push_back( thisRun );
    }
    fcin.close();

    return FCData;
}

// Given a run, this function returns an Input struct based on run number
Input FindRun( vector<Input>& Data, int Run ){

    for(auto Datum : Data) if( Datum.Run == Run ) return Datum;

    cout <<"ERROR: Couldn't find run "<< Run <<". Returning empty Input...\n";
    Input fail;
    return fail; // Return this if run doesn't exist
}

// Used to set the min and max values
void FindBounds( double val, double& min, double& max, vector<double> Vals ){

    int bins = Vals.size()-1;

    for(int i=0; i<bins; i++){
	if( val >= Vals[i] && val < Vals[i+1] ){
	    min = Vals[i];
	    max = Vals[i+1];
	    break;
	}
    }

}

// Pasing in a run number, create an output text file in the target directory 
// Returns true if a file was created
//bool MakeTextFile( int Run, RunPeriod& Period, vector<Input>& FCData, string folderPath ){
bool MakeTextFile( int Run, RunPeriod& Period, string folderPath, DiagnosticPlots& Histos, bool useQA=true ){

    string Target = Period.getTargetType( Run );

    //Input thisData = FindRun( FCData, Run );
    //if( thisData.Run == 0 ) return false; // Didn't find the run

    string fileName;
    if( useQA ) fileName = folderPath + "/ROOT_Files/"+ Target +"_"+ to_string(Run) +"_Data.root";
    else fileName = folderPath + "/ROOT_Files_NoQA/"+ Target +"_"+ to_string(Run) +"_Data_NoQA.root";

    TFile* file = new TFile( fileName.c_str() );
    if( file->IsZombie() ) return false;

    //Helicity  Q2BinData  Q2BinData_HelMinus  Q2BinData_HelPlus  XBinData  XBinData_HelMinus  XBinData_HelPlus

    // Get the data bin information:
    auto X_Bin_Data  = (TProfile2D*)file->Get("XBinData");
    auto Q2_Bin_Data = (TProfile2D*)file->Get("Q2BinData");

    // Get the data bin information for the helicity plus states:
    auto X_Bin_Data_HelPlus  = (TProfile2D*)file->Get("XBinData_HelPlus");
    auto Q2_Bin_Data_HelPlus = (TProfile2D*)file->Get("Q2BinData_HelPlus");

    // Get the data bin information for the helicity minus states:
    auto X_Bin_Data_HelMinus  = (TProfile2D*)file->Get("XBinData_HelMinus");
    auto Q2_Bin_Data_HelMinus = (TProfile2D*)file->Get("Q2BinData_HelMinus");

    // Get the data bin information for the helicity zero states:
    auto X_Bin_Data_HelZero  = (TProfile2D*)file->Get("XBinData_HelZero");
    auto Q2_Bin_Data_HelZero = (TProfile2D*)file->Get("Q2BinData_HelZero");

    // Various kinematic plots used for diagnostics
    auto Theta   = (TH1D*)file->Get("Theta"); Histos.Theta->Add(Theta);
    auto Phi     = (TH1D*)file->Get("Phi"); Histos.Phi->Add(Phi);
    auto Nphe    = (TH1D*)file->Get("nphe"); Histos.Nphe->Add(Nphe);
    auto Vx      = (TH1D*)file->Get("Vx"); Histos.Vx->Add(Vx);
    auto Vy      = (TH1D*)file->Get("Vy"); Histos.Vy->Add(Vy);
    auto Vz      = (TH1D*)file->Get("Vz"); Histos.Vz->Add(Vz);
    auto Q2      = (TH1D*)file->Get("Q2"); Histos.Q2->Add(Q2);
    auto W       = (TH1D*)file->Get("W"); Histos.W->Add(W);
    auto X       = (TH1D*)file->Get("X"); Histos.X->Add(X);
    auto M2_Dist = (TH1D*)file->Get("M2_Dist"); Histos.M2_Dist->Add(M2_Dist);

/*
    auto phi      = (TH1D*)file->Get(""); Phi->Add(phi);
    auto nphe     = (TH1D*)file->Get(""); Nphe->Add(nphe);
    auto vx       = (TH1D*)file->Get(""); Vx->Add(Vx)
    auto vy       = (TH1D*)file->Get("");
    auto vz       = (TH1D*)file->Get("");
    auto beta     = (TH1D*)file->Get("");
    auto q2	  = (TH1D*)file->Get("");
    auto w	  = (TH1D*)file->Get("");
    auto x	  = (TH1D*)file->Get("");
    auto m2_dist  = (TH1D*)file->Get("");
    auto pid      = (TH1D*)file->Get("");
 */
    // At this stage, corrections are applied for the target polarization and the solenoid polarity. The HelP states will now correspond to
    // the beam polarity and target polarity being in opposite directions. Using bt for beam and target polarity, we have the following:
    // --> "HelP" == ++ or --
    // --> "HelN" == +- or -+
    // The FC are swapped based on this configuration as well.
    // "solPol" is the solenoid polarization; negative for negative polarization, and positive for pos. pol.
    /* Temporarily removed for now
    double solPol = Period.getSolenoidScale( Run );
    double tarPol = 1;
    if( Target == "NH3" || Target == "ND3" ) tarPol = Period.getTargetPolarization( Run ); // Only care about this for ammonia
    // Tells whether or not to flip the HelP and HelN counts and FC charges
    bool flipStates = ( tarPol * solPol ) > 0;
    */
    // Now loop over all bins in the profile objects, writing the data for any bins with non-zero counts

    if( Target == "Empty" ) Target = "ET";
    else if( Target == "Foil" ) Target = "F";

    string outName;
    if( useQA ) outName = folderPath +"/Text_Files/"+ Target +"_"+ to_string(Run) +"_DF_Data.txt";
    else outName = folderPath +"/Text_Files_NoQA/"+ Target +"_"+ to_string(Run) +"_DF_Data_NoQA.txt";

    ofstream fout( outName.c_str() );
    fout << "# Q2_Min   Q2_Max   X_Min   X_Max   Np_Counts   Nm_Counts    FCp_Charge   FCm_Charge   N0_Counts   FC0_Charge   X_Avg   Q2_Avg   FC_RUN\n";

    int nXBins  = X_Bin_Bounds.size()-1;  // Number of x bins
    int nQ2Bins = Q2_Bin_Bounds.size()-1; // Number of Q2 bins
    int binMax  = X_Bin_Data->GetBin( nQ2Bins, nXBins ); // The maximum global bin index

    // Get the FC charges
    //double FC_P = thisData.FC_P;
    //double FC_M = thisData.FC_M;

    ROOT::EnableImplicitMT();

    ROOT::RDataFrame df("FCup", fileName.c_str());

    double FC_P   = *df.Sum("fcupgated_HelP");
    double FC_N   = *df.Sum("fcupgated_HelN");
    double FC_0   = *df.Sum("fcupgated_Hel0");
    double FC_RUN = *df.Sum("RUN_fcupgated");

    for( int i=0; i<nQ2Bins; i++ ){
	for( int j=0; j<nXBins; j++ ){

	    int binIndex = X_Bin_Data->GetBin(i,j);

	    // Get the counts for this bin
	    int HelP_Counts = X_Bin_Data_HelPlus->GetBinEntries(  binIndex );
	    int HelM_Counts = X_Bin_Data_HelMinus->GetBinEntries( binIndex );
	    int Hel0_Counts = X_Bin_Data_HelZero->GetBinEntries(  binIndex );

	    // Get the mean values of x and Q2 for this bin
	    double xMean  = X_Bin_Data->GetBinContent(  binIndex );
	    double q2Mean = Q2_Bin_Data->GetBinContent( binIndex );

	    // Get the bin boundaries
	    double Q2Min = 0; double Q2Max = 0; FindBounds( q2Mean, Q2Min, Q2Max, Q2_Bin_Bounds );
	    double XMin = 0;  double XMax  = 0; FindBounds( xMean, XMin, XMax, X_Bin_Bounds );

/*
	    if( flipStates ){
		fout << Q2Min <<"   "<< Q2Max <<"   "<< XMin <<"   "<< XMax <<"   "<< HelM_Counts <<"   "<< HelP_Counts <<"   "<<
		FC_M <<"   "<< FC_P <<"   "<< 0 <<"   "<< 0 <<"   "<< xMean <<"   "<< q2Mean << endl;
	    } 
	    else{
*/
		fout << Q2Min <<"   "<< Q2Max <<"   "<< XMin <<"   "<< XMax <<"   "<< HelP_Counts <<"   "<< HelM_Counts <<"   "<<
		FC_P <<"   "<< FC_N <<"   "<< Hel0_Counts <<"   "<< FC_0 <<"   "<< xMean <<"   "<< q2Mean <<"   "<< FC_RUN << endl;
//	    }
	    
	
	}
	
    }

    fout.close();
    //cout << "Finished making text file for run "<< Run << endl;

    return true;

}

// This function reads in a file and calculates the BSA for the target. Right now, it just checks the 
// sign of the target polarization to see if it's above zero. The signs should already be corrected.
// Returns true if greater than zero.
bool Asymmetry_Sign_Checker( int Run, RunPeriod& Period, string folderPath, bool useQA=true ){

    string Target = Period.getTargetType( Run );

    double Target_Pol = Period.getTargetPolarization( Run );
    if( Run > 17185 && Run < 17450 ) Target_Pol *= -1.0; // Accounts for solenoid polarity flip...

    if( !(Target == "NH3" || Target == "ND3") ) return true;

    if( Target == "Empty" ) Target = "ET";
    else if( Target == "Foil" ) Target = "F";

    string inName;
    if( useQA ) inName = folderPath +"/Text_Files/"+ Target +"_"+ to_string(Run) +"_DF_Data.txt";
    else        inName = folderPath +"/Text_Files_NoQA/"+ Target +"_"+ to_string(Run) +"_DF_Data_NoQA.txt";

    ifstream fin( inName.c_str() );
    // "Q2_Min   Q2_Max   X_Min   X_Max   Np_Counts   Nm_Counts    FCp_Charge   FCm_Charge   N0_Counts   FC0_Charge   X_Avg   Q2_Avg"
    string line;
    getline( fin, line ); // Throw away header row; contents of this row listed above...

    double totalP = 0; double totalM = 0;
    double totFCP = 0; double totFCM = 0;
    while( getline( fin, line ) ){
	stringstream sin(line);
	double Q2_Min, Q2_Max, X_Min, X_Max, Np_Counts, Nm_Counts,  FCp_Charge, FCm_Charge, N0_Counts, FC0_Charge, X_Avg, Q2_Avg, FC_RUN;
	sin >> Q2_Min >> Q2_Max >> X_Min >> X_Max >> Np_Counts >> Nm_Counts >>  FCp_Charge >> FCm_Charge >> N0_Counts >> FC0_Charge >> X_Avg >> Q2_Avg >> FC_RUN;
	
	// Exclude kinematic bins that were outside of the range
	if( Q2_Min > 0 && Q2_Max > 0 && X_Min > 0 && X_Max > 0 ){
	    totalP += Np_Counts; totalM += Nm_Counts;
	    totFCP = FCp_Charge; totFCM = FCm_Charge; // Don't accumulate FC charge since it's common to each bin
	}
    }
    fin.close();

    if( totalP > 0 && totalM > 0 && totFCP > 0 && totFCM > 0 ){
	double normP = totalP / totFCP;
	double normM = totalM / totFCM;
	double rawAsym  = (totalP - totalM) / (totalP + totalM);
	double normAsym = (normP - normM) / (normP + normM);
	//cout << "Run "<< Run <<" raw Asym = "<< rawAsym << endl;
	return Target_Pol * rawAsym > 0 && Target_Pol * normAsym > 0;
    }

    return false;

}


void Make_Text_Files(){

    RunPeriod Period;

    DiagnosticPlots DP_Su22("Su22"), DP_Su22_NoQA("Su22_NoQA");
    DiagnosticPlots DP_Fa22Neg("Fa22Neg"), DP_Fa22Neg_NoQA("Fa22Neg_NoQA");
    DiagnosticPlots DP_Fa22Pos("Fa22Pos"), DP_Fa22Pos_NoQA("Fa22Pos_NoQA");
    DiagnosticPlots DP_Sp23("Sp23"), DP_Sp23_NoQA("Sp23_NoQA");

    // File path to the input root files
    string folderPath = "../Latest_Skims";
    // File path to the Fcup data

    // Now, read in from each of the ROOT files to extract the counts in each
    //bool MakeTextFile( int Run, RunPeriod& Period, vector<Input>& FCData, string folderPath ){
    //vector<int> RunList = { 16137, 16989, 17212, 17584 };
    vector<int> RunList = Period.getAllGoodRuns();
    //vector<int> RunList = Period.getRunEpoch("P03");

    vector<int> MissedRuns, WrongSignRuns;           // Checks if any runs were missed or had the wrong sign for the asymmetry...
    vector<int> MissedRuns_NoQA, WrongSignRuns_NoQA; // Does the same for the no QADB cut data...

    for(int Run : RunList){
/*
        // Set the references based on which run is being analyzed
        DiagnosticPlots& QA_Data, NoQA_Data;
        if( Run < 16800 ){ QA_Data = &DP_Su22; NoQA_Data = &DP_Su22_NoQA; }
        if( Run > 16850 && Run < 17185 ){ QA_Data = &DP_Fa22Neg; NoQA_Data = &DP_Fa22Neg_NoQA; }
        if( Run > 17185 && Run < 17450 ){ QA_Data = &DP_Fa22Pos; NoQA_Data = &DP_Fa22Pos_NoQA; }
        if( Run > 17450 && Run <=17768 ){ QA_Data = &DP_Sp23; NoQA_Data = &DP_Sp23_NoQA; }

	// Create the text files for the data with QADB cuts...
	if( !MakeTextFile( Run, Period, folderPath, QA_Data, true ) )     MissedRuns.push_back( Run );
	if( !Asymmetry_Sign_Checker( Run, Period, folderPath, true ) ) WrongSignRuns.push_back( Run );

	// Create the text files for the data without QADB cuts...
	if( !MakeTextFile( Run, Period, folderPath, NoQA_Data, false ) ) MissedRuns_NoQA.push_back( Run );
	if( !Asymmetry_Sign_Checker( Run, Period, folderPath, false ) )  WrongSignRuns_NoQA.push_back( Run );
*/
	// Create the text files for the data with QADB cuts...
	if( Run < 16850 ) if( !MakeTextFile( Run, Period, folderPath, DP_Su22, true ) )     MissedRuns.push_back( Run );
	if( Run > 16850 && Run < 17185 ) if( !MakeTextFile( Run, Period, folderPath, DP_Fa22Neg, true ) )     MissedRuns.push_back( Run );
	if( Run > 17185 && Run < 17450 ) if( !MakeTextFile( Run, Period, folderPath, DP_Fa22Pos, true ) )     MissedRuns.push_back( Run );
	if( Run > 17450 && Run <=17768 ) if( !MakeTextFile( Run, Period, folderPath, DP_Sp23, true ) )     MissedRuns.push_back( Run );
	if( !Asymmetry_Sign_Checker( Run, Period, folderPath, true ) ) WrongSignRuns.push_back( Run );

	if( Run < 16850 ) if( !MakeTextFile( Run, Period, folderPath, DP_Su22_NoQA, false ) )     MissedRuns_NoQA.push_back( Run );
	if( Run > 16850 && Run < 17185 ) if( !MakeTextFile( Run, Period, folderPath, DP_Fa22Neg_NoQA, false ) )     MissedRuns_NoQA.push_back( Run );
	if( Run > 17185 && Run < 17450 ) if( !MakeTextFile( Run, Period, folderPath, DP_Fa22Pos_NoQA, false ) )     MissedRuns_NoQA.push_back( Run );
	if( Run > 17450 && Run <=17768 ) if( !MakeTextFile( Run, Period, folderPath, DP_Sp23_NoQA, false ) )     MissedRuns_NoQA.push_back( Run );
	if( !Asymmetry_Sign_Checker( Run, Period, folderPath, false ) ) WrongSignRuns_NoQA.push_back( Run );

	cout <<"Finished run "<< Run << endl;
    }

    // Check for any missed runs that weren't opened for the QA cut data...
    if( MissedRuns.size() > 0 ){
	cout <<"Missed "<< MissedRuns.size() <<" run(s) in the file loop.\n";
	cout <<"Runs missed:\n";
	for(int Run : MissedRuns) cout <<"--> Run "<< Run << endl;
	cout << endl;
    }
    // Check for any missed runs that weren't opened for the data without QA cuts...
    if( MissedRuns_NoQA.size() > 0 ){
	cout <<"Missed "<< MissedRuns_NoQA.size() <<" run(s) in the file loop without QA cuts.\n";
	cout <<"Runs missed:\n";
	for(int Run : MissedRuns_NoQA) cout <<"--> Run "<< Run << endl;
	cout << endl;
    }

    // Check for any polarized runs with the wrong sign...
     
    if( WrongSignRuns.size() > 0 ){
	cout <<"There is(are) "<< WrongSignRuns.size() <<" run(s) with the wrong target polarization.\n";
	cout <<"Runs with wrong sign:\n";
	for(int Run : WrongSignRuns) cout <<"--> Run "<< Run << endl;
	cout << endl;
    }
     if( WrongSignRuns_NoQA.size() > 0 ){
	cout <<"There is(are) "<< WrongSignRuns_NoQA.size() <<" run(s) with the wrong target polarization for the data without QA cuts.\n";
	cout <<"Runs with wrong sign:\n";
	for(int Run : WrongSignRuns_NoQA) cout <<"--> Run "<< Run << endl;
	cout << endl;
    }

    vector<DiagnosticPlots> All_QA_Data   = {DP_Su22, DP_Fa22Neg, DP_Fa22Pos, DP_Sp23};
    vector<DiagnosticPlots> All_NoQA_Data = {DP_Su22_NoQA, DP_Fa22Neg_NoQA, DP_Fa22Pos_NoQA, DP_Sp23_NoQA};

    auto Plots     = DrawDiagnosticPlots( All_QA_Data, "QA" );   Plots->Print("PDF_Plots/DiagnosticPlotsQA.pdf");
    auto PlotsNoQA = DrawDiagnosticPlots( All_QA_Data, "NoQA" ); PlotsNoQA->Print("PDF_Plots/DiagnosticPlotsNoQA.pdf");

    cout <<"All done! :)\n";

}
