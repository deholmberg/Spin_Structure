/***********************************************************************************
 *
 * Author: Derek Holmberg
 *
 * Date Created: 9/30/2026
 *
 * Last Modified: 9/30/2026
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

//#include "../Header_Files/Binning_Classes.h"
//#include "../Header_Files/Background_Runs.h"

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

void Calculate_Average_PF(){

    double test = ReturnRunCounts( "NH3", 16137 );
    cout << setprecision(10) << test << endl;

    // Read in the PF data

}
