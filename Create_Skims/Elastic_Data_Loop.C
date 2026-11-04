/***********************************************************************************
 *
 * Author: Derek Holmberg
 *
 * Date Created: 10/2/2026
 *
 * Last Modified: 11/1/2026
 *
 * Purpose:
 * This program parses the ROOT file skims created from the gmn trains from Noemie's
 * elastic analysis code. This program creates a set of ROOT files containing
 * histograms of certain kinematic quantities. These are read-in by the program
 * "Elastic_Diagnostic_Plots.C" to look for differences in kinematic quantities 
 * between run periods. Probably not the best way to do stuff but oh well lol. 
 *
***********************************************************************************/

#include "/w/hallb-scshelf2102/clas12/holmberg/Spin_Structure/Header_Files/Binning_Classes.h"

using namespace std;

string thisDir = "/w/hallb-scshelf2102/clas12/holmberg/Spin_Structure/";

// Loops over a set of input files and writes the histograms to a ROOT file
void ElasticFileLoop( string FilePath, string Filter, int runToAnalyze ){

    ROOT::EnableImplicitMT();

    ROOT::RDataFrame df("Elastic", string(FilePath+".root").c_str());
   
    df.GetColumnNames();

/*
  tree->Branch("Q2", &Q2);
  tree->Branch("delta_Q2", &delta_Q2);
  tree->Branch("W2", &W2); 
  tree->Branch("Emiss", &Emiss);
  tree->Branch("p_perp", &p_perp);
  tree->Branch("p_long", &p_long); 
  tree->Branch("pe_calc", &pe_calc);
  tree->Branch("Ee_calc", &Ee_calc);
  tree->Branch("Ee_calc_purelastic", &Ee_calc_purelastic);
  tree->Branch("pp_calc", &pp_calc);
  tree->Branch("Q2_calc", &Q2_calc);
  tree->Branch("dPhi", &dPhi);
  tree->Branch("delta_pe", &delta_pe);
  tree->Branch("delta_pp", &delta_pp);
  tree->Branch("Ath", &Ath);
*/
    //cout << Filter << endl;
    auto df_filtered = df.Filter( "Helicity==1 || Helicity==-1" );

    // Various kinematic plots used for diagnostics
    auto Ee_calc = df_filtered.Histo1D({"Ee_calc", "Ee_calc Distribution; E_{Beam} [GeV];",2000,0,20},"Ee_calc");
    auto mm2_total = df_filtered.Histo1D({"mm2_total", "mm2_total Distribution; MM^{2} [GeV^{2}];",1000,-0.5,0.5},"mm2_total");
    auto dPhi = df_filtered.Histo1D({"dPhi", "#Delta #phi Distribution; #Delta #phi [deg];",400,170,190},"dPhi");
    auto delta_pp = df_filtered.Histo1D({"delta_pp", "#delta P_{p} Distribution; #delta P_{p} [GeV/c];",2000,-1,1},"delta_pp");
    auto delta_pe = df_filtered.Histo1D({"delta_pe", "#delta P_{e-} Distribution; #delta P_{e-} [GeV/c];",2000,-1,1},"delta_pe");
    auto Q2 = df_filtered.Histo1D({"Q2", "Q^{2} Distribution; Q^{2} [GeV^{2}];",500,0,5},"Q2");
    auto W2 = df_filtered.Histo1D({"W2", "W^{2} Distribution; W^{2} [GeV^{2}];",400,0,4},"W2");

    // Make 2D histos as well
    auto dPhi_vs_dpp = df_filtered.Histo2D({"dPhi_vs_dpp", "#Delta#phi vs. #delta P_{p} Distribution; #Delta P_{p} [GeV/c]; #Delta #phi [deg]",2000,-1,1,400,170,190},"delta_pp","dPhi");
    auto dPhi_vs_dpe = df_filtered.Histo2D({"dPhi_vs_dpe", "#Delta#phi vs. #delta P_{e-} Distribution; #Delta P_{e-} [GeV/c]; #Delta #phi [deg]",2000,-1,1,400,170,190},"delta_pe","dPhi");

    Ee_calc->Write();
    mm2_total->Write();
    dPhi->Write();
    delta_pp->Write();
    delta_pe->Write();
    Q2->Write();
    W2->Write();
    dPhi_vs_dpp->Write();
    dPhi_vs_dpe->Write();

    cout <<"Finished loop on run "<< runToAnalyze <<"...\n";
}

// Takes a vector of runs and writes them to a file
void Write_Data( vector<int>& Runs, string fileExt="_" ){

    for( int Run : Runs ){
        string outFileName, inputFile;
	if( fileExt == "_" ){
	    outFileName = thisDir + "Latest_Skims/Diagnostic_Elastic/Elastic_Diagnostic"+ fileExt + to_string( Run ) +".root";
            inputFile = thisDir + "Latest_Skims/Elastic_ROOT_Files/fiducial_analysis_elastic_"+ to_string( Run ) ;
	}
	else{
	    outFileName = thisDir + "Latest_Skims/Diagnostic_Elastic/Elastic_Diagnostic"+ fileExt + to_string( Run ) +".root";
            inputFile = thisDir + "Latest_Skims/Elastic_ROOT_Files_Pcorr/fiducial"+ fileExt +"elastic_"+ to_string( Run ) ;
	}

	string TotalCut = ""; // No extra cuts for now...
        //string inputFile = thisDir + "Latest_Skims/Elastic_ROOT_Files_Pcorr/fiducial"+ fileExt +"elastic_"+ to_string( Run ) ;
        TFile* outFileTest = new TFile( outFileName.c_str() , "RECREATE" );
        ElasticFileLoop( inputFile, TotalCut, Run );
        outFileTest->Close();
    }

}

void Elastic_Data_Loop(){

  int run;
  cin >> run;
  vector<int> Runs = {run};

  string fileExt = "_correctedP_";

  //Write_Data( Runs, fileExt );
  Write_Data( Runs );

  // Start a loop over all runs:
  /*
  RunPeriod Period;

  
  for( int i=1; i<2; i++ ){

    string epoch = "P";
    if( i<10 ) epoch += "0"+to_string(i);
    else epoch += to_string(i);

    vector<int> Epoch_Runs = Period.getRunEpoch( epoch );

    Write_Data( Epoch_Runs );

  } // End of file loop

  // Now loop over the background runs
  
  Write_Data( Su22_CH2_Runs );
  Write_Data( Su22_F_Runs );
  Write_Data( Su22_ET_Runs );
  Write_Data( Su22_C_Runs );

  Write_Data( Fa22Neg_CH2_Runs );
  Write_Data( Fa22Neg_ET_Runs );
  Write_Data( Fa22Neg_F_Runs );
  Write_Data( Fa22Neg_C_Runs );

  Write_Data( Fa22Pos_CH2_Runs );
  Write_Data( Fa22Pos_C_Runs );

  Write_Data( Sp23Inb_CH2_Runs );
  Write_Data( Sp23Inb_ET_Runs );
  Write_Data( Sp23Inb_F_Runs );
  Write_Data( Sp23Inb_C_Runs );
  Write_Data( Sp23Inb_CD2_Runs );
  */
  cout <<"All done! :)\n";

}
