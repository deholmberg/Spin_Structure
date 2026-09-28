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

void Calculate_PbPt(){

    RunPeriod Period; // Holds the run info

    DataSet Test;
    Test.ReadDFfromTXT("Output_Data/All_DF_Data_Epoch_10.txt");
    //Test.Print(true);
    if( Test.AddToBins("../Latest_Skims/Text_Files/NH3_16772_DF_Data.txt", "NH3") ) cout << "Opened run 16772!\n";
    Test.MaxLikelihoodPbPt( 16772, Period, true );

}
