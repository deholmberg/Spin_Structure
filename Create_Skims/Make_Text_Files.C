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
bool MakeTextFile( int Run, RunPeriod& Period, vector<Input>& FCData, string folderPath ){

    string Target = Period.getTargetType( Run );

    Input thisData = FindRun( FCData, Run );
    if( thisData.Run == 0 ) return false; // Didn't find the run

    string fileName = folderPath + "/ROOT_Files/"+ Target +"_"+ to_string(Run) +"_Data.root";
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

    // At this stage, corrections are applied for the target polarization and the solenoid polarity. The HelP states will now correspond to
    // the beam polarity and target polarity being in opposite directions. Using bt for beam and target polarity, we have the following:
    // --> "HelP" == ++ or --
    // --> "HelN" == +- or -+
    // The FC are swapped based on this configuration as well.
    // "solPol" is the solenoid polarization; negative for negative polarization, and positive for pos. pol.
    double solPol = Period.getSolenoidScale( Run );
    double tarPol = 1;
    if( Target == "NH3" || Target == "ND3" ) tarPol = Period.getTargetPolarization( Run ); // Only care about this for ammonia
    // Tells whether or not to flip the HelP and HelN counts and FC charges
    bool flipStates = ( tarPol * solPol ) > 0;

    // Now loop over all bins in the profile objects, writing the data for any bins with non-zero counts

    if( Target == "Empty" ) Target = "ET";
    else if( Target == "Foil" ) Target = "F";

    string outName = folderPath +"/Text_Files/"+ Target +"_"+ to_string(Run) +"_DF_Data.txt";
    ofstream fout( outName.c_str() );
    fout << "Q2_Min   Q2_Max   X_Min   X_Max   Np_Counts   Nm_Counts    FCp_Charge   FCm_Charge   N0_Counts   FC0_Charge   X_Avg   Q2_Avg\n";

    int nXBins  = X_Bin_Bounds.size()-1;  // Number of x bins
    int nQ2Bins = Q2_Bin_Bounds.size()-1; // Number of Q2 bins
    int binMax  = X_Bin_Data->GetBin( nQ2Bins, nXBins ); // The maximum global bin index

    for( int i=0; i<nQ2Bins; i++ ){
	for( int j=0; j<nXBins; j++ ){

	    int binIndex = X_Bin_Data->GetBin(i,j);

	    // Get the counts for this bin
	    int HelP_Counts = X_Bin_Data_HelPlus->GetBinEntries(  binIndex );
	    int HelM_Counts = X_Bin_Data_HelMinus->GetBinEntries( binIndex );

	    // Get the mean values of x and Q2 for this bin
	    double xMean  = X_Bin_Data->GetBinContent(  binIndex );
	    double q2Mean = Q2_Bin_Data->GetBinContent( binIndex );

	    // Get the bin boundaries
	    double Q2Min = 0; double Q2Max = 0; FindBounds( q2Mean, Q2Min, Q2Max, Q2_Bin_Bounds );
	    double XMin = 0;  double XMax  = 0; FindBounds( xMean, XMin, XMax, X_Bin_Bounds );

	    // Get the FC charges
	    double FC_P = thisData.FC_P;
	    double FC_M = thisData.FC_M;

	    // If there's data, write the bin...
	    //cout << Q2Min <<"   "<< Q2Max <<"   "<< XMin <<"   "<< XMax <<"   "<< HelP_Counts <<"   "<< HelM_Counts <<"   "<<
	    //	    FC_P <<"   "<< FC_M <<"   "<< 0 <<"   "<< 0 <<"   "<< xMean <<"   "<< q2Mean << endl;    

	    //if( HelP_Counts > 0 && HelM_Counts > 0 && Q2Min > 0 && Q2Max > 0 && XMin > 0 && XMax > 0 ){

		if( flipStates ){
			fout << Q2Min <<"   "<< Q2Max <<"   "<< XMin <<"   "<< XMax <<"   "<< HelM_Counts <<"   "<< HelP_Counts <<"   "<<
			FC_M <<"   "<< FC_P <<"   "<< 0 <<"   "<< 0 <<"   "<< xMean <<"   "<< q2Mean << endl;
		} 
		else{
			fout << Q2Min <<"   "<< Q2Max <<"   "<< XMin <<"   "<< XMax <<"   "<< HelP_Counts <<"   "<< HelM_Counts <<"   "<<
			FC_P <<"   "<< FC_M <<"   "<< 0 <<"   "<< 0 <<"   "<< xMean <<"   "<< q2Mean << endl;
		}
	    //}
	
	}
	
    }

    fout.close();
    //cout << "Finished making text file for run "<< Run << endl;

    return true;

}

// This function reads in a file and calculates the BSA for the target. Right now, it just checks the 
// sign of the target polarization to see if it's above zero. The signs should already be corrected.
// Returns true if greater than zero.
bool Asymmetry_Sign_Checker( int Run, RunPeriod& Period, string folderPath ){

    string Target = Period.getTargetType( Run );

    if( !(Target == "NH3" || Target == "ND3") ) return true;

    if( Target == "Empty" ) Target = "ET";
    else if( Target == "Foil" ) Target = "F";

    string inName = folderPath +"/Text_Files/"+ Target +"_"+ to_string(Run) +"_DF_Data.txt";
    ifstream fin( inName.c_str() );
    // "Q2_Min   Q2_Max   X_Min   X_Max   Np_Counts   Nm_Counts    FCp_Charge   FCm_Charge   N0_Counts   FC0_Charge   X_Avg   Q2_Avg"
    string line;
    getline( fin, line ); // Throw away header row; contents of this row listed above...

    double totalP = 0; double totalM = 0;
    double totFCP = 0; double totFCM = 0;
    while( getline( fin, line ) ){
	stringstream sin(line);
	double Q2_Min, Q2_Max, X_Min, X_Max, Np_Counts, Nm_Counts,  FCp_Charge, FCm_Charge, N0_Counts, FC0_Charge, X_Avg, Q2_Avg;
	sin >> Q2_Min >> Q2_Max >> X_Min >> X_Max >> Np_Counts >> Nm_Counts >>  FCp_Charge >> FCm_Charge >> N0_Counts >> FC0_Charge >> X_Avg >> Q2_Avg;
	
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
	return rawAsym > 0 && normAsym > 0;
    }

    return false;

}


void Make_Text_Files(){

    RunPeriod Period;

    // File path to the input root files
    string folderPath = "../Latest_Skims/";
    // File path to the Fcup data
    string fcupPath = "../../Dilution_Factors/LatestSkim/FCup_Data/FCup_Data_WithQA.txt";


    // Test to see if runs were read in successfully
    vector<Input> FCData = Create_FCData( fcupPath );
    
    //for( auto Datum : FCData ) Datum.Print();

    // Now, read in from each of the ROOT files to extract the counts in each
    //bool MakeTextFile( int Run, RunPeriod& Period, vector<Input>& FCData, string folderPath ){
    //vector<int> RunList = { 16137 };
    vector<int> RunList = Period.getAllGoodRuns();

    vector<int> MissedRuns, WrongSignRuns;

    for(int Run : RunList){
	if( !MakeTextFile( Run, Period, FCData, folderPath ) ) MissedRuns.push_back( Run );
	if( !Asymmetry_Sign_Checker( Run, Period,  folderPath ) ) WrongSignRuns.push_back( Run );
    }

    // Check for any missed runs that weren't opened...
    if( MissedRuns.size() > 0 ){
	cout <<"Missed "<< MissedRuns.size() <<" run(s) in the file loop.\n";
	cout <<"Runs missed:\n";
	for(int Run : MissedRuns) cout <<"--> Run "<< Run << endl;
	cout << endl;
    }

    // Check for any polarized runs with 
    if( WrongSignRuns.size() > 0 ){
	cout <<"There is(are) "<< WrongSignRuns.size() <<" run(s) with the wrong target polarization.\n";
	cout <<"Runs with wrong sign:\n";
	for(int Run : WrongSignRuns) cout <<"--> Run "<< Run << endl;
	cout << endl;
    }

    cout <<"All done! :)\n";

}
