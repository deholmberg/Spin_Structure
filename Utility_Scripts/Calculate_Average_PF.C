/***********************************************************************************
 *
 * Author: Derek Holmberg
 *
 * Date Created: 9/30/2026
 *
 * Last Modified: 10/6/2026
 *
 * Purpose:
 * This program calculates the statistics-weighted average of all packing fractions
 * from the "../Analysis_Channel/Output_Data/All_NH3_PF_Data.txt" file. The outputs
 * from this program are used in a modified version of Noemie's code to calculate
 * dilution factors for the elastic electron-proton scattering channel. Calculations
 * are broken down based on these divisions:
 * 
 * 1) Run period (summer, fallneg, etc.)
 * 2) Target polarizations (positive or negative)
 *
 * There are eight sub-periods in all (four inbending run periods with two target
 * polarizations). Due to the tight exclusivity cuts in the elastic analysis, the
 * cell value of the PF is the relevant quantity to calculate.
 *
***********************************************************************************/

#include "../Header_Files/Binning_Classes.h"
#include "../Header_Files/Background_Runs.h"

using namespace std;

// For a given run, it opens the associated ROOT file with all DIS/fiducial cuts applied 
// and returns the total counts from helicity plus or minus states.
double ReturnRunCounts( string Target, int Run ){

    string fileName = "../Latest_Skims/ROOT_Files/"+ Target +"_"+ to_string(Run) +"_Data.root";
    TFile* file = new TFile( fileName.c_str() );
    if( file->IsZombie() ) return 0;

    //Helicity  Q2BinData  Q2BinData_HelMinus  Q2BinData_HelPlus  XBinData  XBinData_HelMinus  XBinData_HelPlus
    auto Helicity = (TH1D*)file->Get("Helicity");
    // Get the counts for each helicity state
    double Helicity_Plus  = Helicity->GetBinContent(7);
    double Helicity_Minus = Helicity->GetBinContent(3);

    return Helicity_Plus + Helicity_Minus;

}

// Just a dumb struct to hold the PF info for each epoch
struct PF{

    string Epoch        ="";
    double PF_Bath      = 0;
    double PF_Bath_Err  = 0;
    double PF_Bath_Chi2 = 0;
    double PF_Cell      = 0;
    double PF_Cell_Err  = 0;
    double PF_Cell_Chi2 = 0;

    // Used in the script to set the values
    void SetPFs( stringstream& s, string Target="NH3" ){
	int epoch;
	s >> epoch >> PF_Bath >> PF_Bath_Err >> PF_Bath_Chi2
	  >> PF_Cell >> PF_Cell_Err >> PF_Cell_Chi2;
	// Convert "epoch" to a string
	if( epoch < 10 && Target == "NH3" ) Epoch = "P0" + to_string(epoch);
	else if( Target == "NH3" )          Epoch = "P"  + to_string(epoch);
	else if( epoch < 10 && Target == "ND3" ) Epoch = "D0" + to_string(epoch);
	else if( Target == "ND3" )               Epoch = "D"  + to_string(epoch);
    }
    void Print() const{
	cout <<"========== Epoch "<< Epoch <<" Info ==========\n";
	cout <<"--> PF_Bath = "<< PF_Bath <<" +- "<< PF_Bath_Err <<", Chi2_red = "<< PF_Bath_Chi2 << endl;
	cout <<"--> PF_Cell = "<< PF_Cell <<" +- "<< PF_Cell_Err <<", Chi2_red = "<< PF_Cell_Chi2 << endl;
	cout <<"===================================\n";
    }

};

// Gets the PF data and returns a vector holding the PF info
vector<PF> GetPFs( string TargetType="NH3" ){

    vector<PF> PFs;

    // Read in the PF data
    ifstream fin("../Analysis_Channel/Output_Data/All_NH3_PF_Data.txt");

    if( fin.fail() ){ cout <<"ERROR: Couldn't find PF file. Check inputs.\n"; return PFs; }

    string line;
    while( getline( fin, line ) ){
	stringstream sin(line);
	PF thisPF; thisPF.SetPFs( sin, TargetType );
	PFs.push_back(thisPF);
	// Test if read-in was successful
	//thisPF.Print();
    }    
    fin.close();

    return PFs;

}

// This function calculates the weighted average of all PF across a given run period for 
// each DIS epoch, using the statistical error in PF as the weighted value. Currently, this is 
// only calculated for the cell value
vector<double> CalculateWeightedAvgPF( vector<PF>& PFs, string Period ){

    // Set the upper and lower epoch IDs to loop over for the averaging calculation...
    string lowBound = ""; string upperBound = ""; 
    if( Period == "Su22"    ){ lowBound = "P01"; upperBound = "P10"; }
    if( Period == "Fa22Neg" ){ lowBound = "P11"; upperBound = "P15"; }
    if( Period == "Fa22Pos" ){ lowBound = "P16"; upperBound = "P19"; }
    if( Period == "Sp23"    ){ lowBound = "P20"; upperBound = "P23"; }

    // Numerator and denominator terms for calculating the weighted average
    double pfCalcNum = 0; double pfCalcDen = 0;
    for(auto pf : PFs){
	if( pf.Epoch >= lowBound && pf.Epoch <= upperBound ){
	    cout <<"Added epoch "<< pf.Epoch << endl;
	    pfCalcNum += pf.PF_Cell / pow( pf.PF_Cell_Err, 2 );
	    pfCalcDen += 1.0 / pow( pf.PF_Cell_Err, 2 );
	}
    }

    vector<double> PF_Vals = {0,0};
    if( pfCalcDen > 0 ){
	PF_Vals[0] = pfCalcNum / pfCalcDen;
	PF_Vals[1] = 1.0 / sqrt( pfCalcDen );
	cout <<"For PFs in the "<< Period <<" run period, avg. PF_cell = "<< PF_Vals[0] <<" +- "<< PF_Vals[1] << endl;
    }

    return PF_Vals;

}

// This function calculates the statistics-weighted average of all the PF for a given elastic
// epoch, weighting the PFs by the total counts from the runs
vector<double> CalculateStatAvgPF( RunPeriod& RP, vector<PF>& PFs, string Period, string Target, int targetPol ){

    // Get the runs for this period
    auto Runs = RP.getElasticEpoch( Period, Target, targetPol );

    // Tracks the missed runs
    vector<int> MissedRuns;

    // For each run, get the total counts and start computing the stat-weighted average
    double totalCounts = 0; double avgPF = 0;

    for( int run : Runs ){
	// Find the PF for this run
	double thisPF = 0;
	string thisEp = RP.getEpoch( run ); // Returns the string of the epoch number
	// Make sure to use the cell PF for the elastic calculation!
	for( auto PF : PFs ) if( PF.Epoch == thisEp ) thisPF = PF.PF_Cell;
	
	double runCounts = ReturnRunCounts( Target, run );
	if( runCounts > 0 && thisPF > 0 ){
	    totalCounts += runCounts;
	    avgPF += runCounts * thisPF;
	}
	else MissedRuns.push_back( run );
    }

    // Print missed runs
    if( MissedRuns.size() > 0 ){
	cout <<"Missed runs or runs with zero events:\n";
	for(int r : MissedRuns) cout <<"--> "<< r << endl;
    }

    vector<double> PF_Vals = {0,0};

    if( totalCounts > 0 ){
	PF_Vals[0] = avgPF / totalCounts;
	PF_Vals[1] = 1.0 / sqrt( totalCounts );
	cout <<"For "<< Period <<" in the "<< targetPol <<" polarization, avg. PF_cell = "<< PF_Vals[0] <<" +- "<< PF_Vals[1] << endl;
    }

    return PF_Vals; 

}
//Bin_Q2 Bin_X DF_NH3 Err_DF_NH3 PF_bath_NH3 ErrPF_bath_NH3 PF_cell_NH3 ErrPF_cell_NH3 FC_total N_NH3 N_C N_CH2 N_ET N_F N_CD2 FC_NH3 FC_C FC_CH2 FC_ET FC_F FC_CD2
//1.86303 0.0875 0.869179 0.00832395 0.490623 0.0178347 0.57501 0.0209023 2.99457e+06 13552 270 1031 102.402 1 0 2.05138e+06 214938 148785 116298 463163 0

// Calculate the average value of the elastic PF for a given run period and target polarization
vector<double> CalculateAvgElasticPF( RunPeriod& RP, string Period, string PosOrNeg ){

    vector<double> PF_Vals = {0,0};

    // Open the elastic data file
    //Elastic_DF_Data_Neg_Su22.txt
    ifstream fin( string("../Analysis_Channel/Output_Data/Elastic_DF_Data_"+PosOrNeg+"_"+Period+".txt") );
    if( fin.fail() ){ cout <<"ERROR: Couldn't find elastic PF file. Check inputs.\n"; return PF_Vals; }

    // Declare the values for calculating the weighted average for the cell
    double pfCalcNum = 0; double pfCalcDen = 0;

    string line;
    getline( fin, line ); // Throw away header row
    while( getline( fin, line ) ){
	stringstream sin(line);
	double Bin_Q2, Bin_X, DF_NH3, Err_DF_NH3, PF_bath_NH3, ErrPF_bath_NH3, PF_cell_NH3, ErrPF_cell_NH3;
	sin >> Bin_Q2 >> Bin_X >> DF_NH3 >> Err_DF_NH3 >> PF_bath_NH3 >> ErrPF_bath_NH3 >> PF_cell_NH3 >> ErrPF_cell_NH3;

	if( PF_cell_NH3 > 0 && ErrPF_cell_NH3 > 0 ){
	    pfCalcNum += PF_cell_NH3 / pow( ErrPF_cell_NH3, 2 );
	    pfCalcDen += 1.0 / pow( ErrPF_cell_NH3, 2 );
	}

    }
    fin.close();

    if( pfCalcDen > 0 ){
	PF_Vals[0] = pfCalcNum / pfCalcDen;
	PF_Vals[1] = 1.0 / sqrt( pfCalcDen );
	cout <<"Elastic average PF for "<< Period <<" "<< PosOrNeg <<" is "<< PF_Vals[0] <<" +- "<< PF_Vals[1] << endl;
    }

    return PF_Vals;

}

void Calculate_Average_PF(){

    //double test = ReturnRunCounts( "NH3", 16137 );
    //cout << setprecision(10) << test << endl;

    // Get a vector of the DIS PF data
    auto PFs = GetPFs();

    RunPeriod Period;

    // Statistical weighted method

    auto Su22_Pos = CalculateStatAvgPF( Period, PFs, "Su22", "NH3", 1 );
    auto Su22_Neg = CalculateStatAvgPF( Period, PFs, "Su22", "NH3",-1 );

    // Weighted average method
    auto Su22_Stat_Pos = CalculateWeightedAvgPF( PFs, "Su22" );

    // Average elastic PF
    auto Su22_Pos_Elastic = CalculateAvgElasticPF( Period, "Su22", "Pos" );
    auto Su22_Neg_Elastic = CalculateAvgElasticPF( Period, "Su22", "Neg" );

//vector<double> CalculateAvgElasticPF( RunPeriod& RP, string Period, string PosOrNeg ){

}
