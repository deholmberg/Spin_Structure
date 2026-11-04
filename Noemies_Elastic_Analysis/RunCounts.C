#include <dirent.h>
#include <boost/filesystem.hpp>
//If you don't have boost, the code runs fine
//with std::filesystem instead
//#include <filesystem>
//#include "../utils/utils.C"
//#include "../utils/cosmetics.C"
#include "../Header_Files/Binning_Classes.h"

double Ebeam = 10.6;
double Pmass = 938.272081/1000.;

string dirPath     = "/w/hallb-scshelf2102/clas12/holmberg/Spin_Structure/Latest_Skims/";
string dirPathCuts = "/w/hallb-scshelf2102/clas12/holmberg/Spin_Structure/Noemies_Elastic_Analysis/";

//-------------------------------------------------------------------------------Define cuts---------------------------------------------------------------------------------------//
// This function returns an integer containing all the error bits to be set for a
// given run.
int ReturnErrorBits( RunPeriod& RP, int Run ){

    int ErrorBits = 0;

    bool SetMiscBit = RP.applyMiscBit( Run );

    // Binary holding all the error bits. Flip the misc bit if it's applicable
    bitset<21> SetBits(0b011110000110000011111);
    if( SetMiscBit ) SetBits.flip(5); // Flip bit 5 for misc

    ErrorBits = SetBits.to_ulong(); // Truncated to an int

    return ErrorBits;

}

TCut exclu_cuts(string target, int ErrorBits)
{
  target = "NH3";
  //This function reads exclusivity cuts
  //Input: string for target type
  //Ouput: TCut holding the exclusivity cuts
//  cout << "Exclusivity cuts" << endl;                                                                                                                                                                    
  TCut cut_all_topo;
  // Open file streams for CD and FD cuts files
  string path_CD = dirPathCuts + target+"_CD_cuts.txt";
  string path_FD = dirPathCuts + target+"_FD_cuts.txt";
  fstream file_CD, file_FD;
  file_CD.open(path_CD, ios::in);
  file_FD.open(path_FD, ios::in);
  //Check files sanity

  if( file_CD.fail() ) cout <<"OOHHHHH NNNNNOOOOOOOO!!!!!\n";

  if (file_CD.bad())
    {
      cerr << "Unable to open " + path_CD + "\n";
      exit(8);
    }
  if (file_FD.bad())
    {
      cerr << "Unable to open " + path_FD + "\n";
      exit(8);
    }
//  cout << ">>> Cut files opened " << path_CD << " and " << path_FD << endl;
  string CD,FD;
  TCut cut_CD, cut_FD;
  // Define topological cuts for CD and FD
  TCut topo_CD = "Nuc_status >= 4000 && Nuc_status < 8000";
  TCut topo_FD = "Nuc_status < 4000 && Nuc_status >=2000";
  //Read exclu cuts
  getline(file_CD, CD);
  cut_CD = CD.c_str();
  file_CD.close();
  getline(file_FD, FD);
  cut_FD = FD.c_str();
  file_FD.close();
  //Build final cut to be applied

  // if it's equal to zero, it passes the check
  string S1_QADB  = "(qadb_mask_s1 & "+to_string(ErrorBits)+") == 0";
  string S2_QADB  = "(qadb_mask_s2 & "+to_string(ErrorBits)+") == 0";
  string S3_QADB  = "(qadb_mask_s3 & "+to_string(ErrorBits)+") == 0";
  string S4_QADB  = "(qadb_mask_s4 & "+to_string(ErrorBits)+") == 0";
  string S5_QADB  = "(qadb_mask_s5 & "+to_string(ErrorBits)+") == 0";
  string S6_QADB  = "(qadb_mask_s6 & "+to_string(ErrorBits)+") == 0";
  TCut s1_qadb = S1_QADB.c_str();
  TCut s2_qadb = S2_QADB.c_str();
  TCut s3_qadb = S3_QADB.c_str();
  TCut s4_qadb = S4_QADB.c_str();
  TCut s5_qadb = S5_QADB.c_str();
  TCut s6_qadb = S6_QADB.c_str();
  //string All_QADB = "("+ S1_QADB +")&&("+ S2_QADB +")&&("+ S3_QADB +")&&("+ S4_QADB +")&&("+ S5_QADB +")&&("+ S6_QADB +")";
  TCut All_QADB = s1_qadb + s2_qadb + s3_qadb + s4_qadb + s5_qadb + s6_qadb;

  cut_all_topo = (cut_CD && topo_CD) || (cut_FD && topo_FD);
  cut_all_topo =  cut_all_topo + "El_fiducial && Nuc_chi2pid<3.4 && Nuc_chi2pid>-2.6 && TMath::Abs(El_chi2pid)<3"+ All_QADB;
//cout << ">>> Cuts read" << cut_all_topo << endl;
  return cut_all_topo;
}

vector<int> runs= {16194, 16872, 16976, 16978, 16979, 17482, 17763, 17764, 17765,

16309, 16970, 16975, 17766, 17767, 17768,

17514, 17516, 17517, 17518, 17520, 17668, 17669, 17670, 17671, 17672, 17673, 17674, 17675, 17676,

16128, 16134, 16290, 16291, 16292, 16293, 16296, 16297, 16700, 16701, 16702, 16704, 16859, 16860, 16874, 16877, 16878, 17138, 17139, 17140, 17170, 17171, 17173, 17177, 17178, 17179, 17180, 17181, 17182, 17183, 17304, 17305, 17306, 17307, 17308, 17310, 17312, 17313, 17314, 17316, 17317, 17319, 17320, 17321, 17322, 17323, 17324, 17325, 17485, 17486, 17487, 17562, 17566, 17567, 17568, 17569, 17570, 17742, 17744, 17745, 17746, 17747, 17805, 17806, 17807, 17808, 17809, 17810, 17811, 

16298, 16299, 16300, 16301, 16302, 16303, 17118, 17119, 17122, 17123, 17125, 17129, 17130, 17131, 17132, 17133, 17134, 17135, 17328, 17329, 17330, 17331, 17332, 17333, 17334, 17340, 17348, 17387, 17388, 17393, 17394, 17397, 17398, 17399, 17400, 17401, 17402, 17403, 17404, 17405, 17406, 17407, 17408, 17624, 17625, 17626, 17627, 17628, 17629, 17630, 17631, 17632, 17633, 17635, 17636, 17637, 17638, 17639, 17640, 17641, 17642, 17643, 17644,

16137, 16138, 16144, 16145, 16146, 16148, 16156, 16157, 16158, 16166, 16167, 16168, 16169, 16170, 16178, 16213, 16214, 16221,
16222, 16223, 16224, 16225, 16226, 16228, 16231, 16232, 16233, 16238, 16243, 16244, 16245, 16246, 16248, 16249, 16252, 16253, 16256, 16257, 16259,
16260, 16318,
16320, 16321, 16322, 16323, 16326, 16327, 16328, 16329, 16330, 16331, 16332, 16333, 16335, 16336, 16337, 16338, 16339, 16341, 16343, 16345, 16346,
16348, 16350, 16352, 16353, 16354, 16355, 16356, 16357, 16658,
16659, 16660, 16664, 16665, 16666, 16671, 16672, 16673, 16674, 16675, 16676, 16678, 16679, 16681, 16682, 16683, 16685, 16686, 16687, 16688, 16689,
16690, 16692, 16693, 16695, 16709, 16710, 16711, 16712, 16713, 16715, 16716, 16717, 16718, 16719, 16720, 16723, 16726,
16727, 16728, 16729, 16730, 16731, 16732, 16733, 16734, 16736, 16738, 16742, 16743, 16744, 16746, 16747, 16748, 16749, 16750, 16751, 16752, 16753,
16754, 16755, 16756, 16757, 16758, 16759, 16761, 16762, 16763, 16765, 16766, 16767, 16768, 16769, 16770, 16771, 16772, 16983,
16984, 16986, 16988, 16989, 16990, 16991, 16992, 16993, 16995, 16996, 17000, 17001, 17002, 17003, 17008, 17010, 17011, 17014, 17015, 17016, 17019,
17020, 17022, 17025, 17026, 17027, 17028, 17032, 17067, 

17069, 17070, 17071, 17073, 17074, 17075, 17076, 17077, 17078, 17084, 17085, 17086, 17087, 17088, 17092, 17093, 17095, 17096, 17097, 17098, 17100,
17101, 17102, 17144, 17146, 17147, 17149,
17150, 17151, 17153, 17155, 17156, 17157, 17158, 17159, 17160, 17161, 17162, 17163, 17164, 17165, 17166, 17167, 17168, 17169, 17188,
17189, 17191, 17192, 17195, 17197, 17199, 17203, 17204, 17206, 17207, 17211, 17212, 17214, 17215, 17216, 17221, 17222, 17223, 17224, 17225, 17353,
17355, 17356, 17357, 17358, 17359, 17360, 17361, 17362, 17363, 17364, 17366, 17368, 17371, 17373, 17374, 17375, 17376, 17379, 17381, 17382, 17575,
17576, 17577, 17581, 17583, 17584, 17585, 17586, 17587, 17588, 17589, 17593, 17594, 17595, 17596, 17597, 17598, 17599, 17600, 17601, 17602, 17603,
17604, 17605, 17606, 17607, 17608, 17609, 17610, 17611, 17720,
17721, 17724, 17725, 17726, 17727, 17728, 17732, 17733, 17734, 17737, 17738, 17740, 17741, 17748, 17750, 17751,
17752, 17753, 17754, 17755, 17756, 17757, 17758, 17759, 17760, 17762, 17769, 17770, 17771, 17772, 17773,
17774, 17775, 17776, 17777, 17778, 
};


//-----------------------------------------------------------------------------Propagating stat errors-------------------------------------------------------------------//
double error_asym(double Np, double Nm, double dNp, double dNm, double Df=1, double dDf=0)
{
  //This function computes the stat error including the error on the dilution factor A = (Np-Nm)/(Df*(Np+Nm))
  //Input: (double) FCup Normalized yields and their errors, Dilution Factor and its error
  //Ouput: (double) stat error
  double den = Df*pow((Np + Nm),2); 
  double term_Np = 2*Nm / den;        
  double term_Nm = 2*Np / den;
  double asym = (Np-Nm)/(Df*(Np+Nm));
  double term_Df = asym/Df;
  double err_squared = term_Np*term_Np*dNp*dNp + term_Nm*term_Nm*dNm*dNm + term_Df*term_Df*dDf*dDf;
  return TMath::Sqrt(err_squared);
}

//-----------------------------------------------------------------------------------------MAIN-------------------------------------------------------------------------------------//
//int RunCounts(string analysis_folder, string filename_C, string filename_signal, string target, string outFileName="")
//int RunCounts(string analysis_folder, int runToAnalyze, string target, string outFileName="")
void RunCounts()
{
  // This is a very slightly modified version of Noemie's script that returns the total number of elastic counts on a run-by-run basis.
  // Everything else besides that is identical, just removing stuff for calculating the dilution factor, PbPt, etc.

  //Main function to compute PbPt
  //Input: (string) filenames of the files containing the list of runs, (string) target type
  //Output: (int) sanity check
  //Results from the analysis are saved in the /results folder
  //It contains: - PbPt plot with the Q^2 dependence of the asymmetries and their ratio
  //             - Dilution factor plot
  //             - Text file with the PbPt result

  string target, outFileName;
  //int runToAnalyze;
  //cin >> runToAnalyze;

  RunPeriod RP;

 for(int runToAnalyze : runs ){

  target = RP.getTargetType( runToAnalyze );
  if( target == "Empty" ) target = "ET";
  else if( target == "Foil" ) target = "F";

  // Update the beam energy based on the run number
  if( runToAnalyze <= 17065 ) Ebeam = beam_energy1;
  else if( runToAnalyze >= 17067 && runToAnalyze <= 17704 ) Ebeam = beam_energy2;
  else if( runToAnalyze >= 17720 ) Ebeam = beam_energy3;

  // Set the error bits for this run:
  int ErrorBits = ReturnErrorBits( RP, runToAnalyze );

  // Create an output text file with the mean Q2 values, counts, FC charge, elastic DFs, theoretical asymmetries:
  //string dirPath = "/w/hallb-scshelf2102/clas12/holmberg/PairSymBkg/shared_for_RGC/PbPt_elastic/";
  //string dirPath = "/w/hallb-scshelf2102/clas12/holmberg/Spin_Structure/Latest_Skims/";

  ofstream fout(string(dirPath + "Elastic_Text_Files/"+ target +"_"+to_string(runToAnalyze)+"_Elastic.txt"));
  //ofstream fout(string(dirPath + "Elastic_Text_Files_Pcorr/"+ target +"_"+to_string(runToAnalyze)+"_Elastic.txt"));

  fout << "Q2_Mean   Theta_Mean   A_el   N+   N-   FC+   FC-\n";
  
  vector<int> run_signal = { runToAnalyze };
    
  //-------------------------------------------------------------------------------Open files and get trees-------------------------------------------------------------------------//
  cout << "Opening input files for run " << runToAnalyze << endl;

  vector<int> valid_runs_ND3, valid_runs_C;
  //SIGNAL
  //TChain *REC_signal = new TChain("Elastic_loose");
  //string analysis_folder = "/volatile/clas12/holmberg/ElasticData";
  string analysis_folder = dirPath + "Elastic_ROOT_Files"; 
  //string analysis_folder = dirPath + "Elastic_ROOT_Files_Pcorr";
 
  TChain *REC_signal = new TChain("Elastic");
  TChain *scaler_signal = new TChain("Scaler info");
  for (int i_files = 0; i_files < run_signal.size(); i_files++)       
    {                     
      int run_number = run_signal.at(i_files);
      if(run_number == 0) {cerr<<"RUNS READ INCORRECTLY!"<<endl; return;}

      string filename_signal = analysis_folder + "/fiducial_analysis_elastic_" + to_string(run_number) + ".root";
      //string filename_signal = analysis_folder + "/fiducial_correctedP_elastic_" + to_string(run_number) + ".root";

      //Check if the file exists
      if (boost::filesystem::exists(filename_signal))  
        {                         
          REC_signal->Add(filename_signal.c_str());    
          scaler_signal->Add(filename_signal.c_str());
	  valid_runs_ND3.push_back(run_number);
        }       
      else {cout << filename_signal << " does not exist, continuing" << endl;}     
    }

  //---------------------------------------------------------------------------Get FCup normalizations-------------------------------------------------------------------------------//
//  cout << "FCup normalization" << endl;
  //Adding normalizations between runs
  double FCup_p_signal = 0;
  double FCup_n_signal = 0;
  double FCup_run_signal = 0;
  double FCup_hel_p_signal, FCup_hel_n_signal, FCup_read_run_signal;
  scaler_signal->SetBranchAddress("FCup_hel_p", &FCup_hel_p_signal);
  scaler_signal->SetBranchAddress("FCup_hel_n", &FCup_hel_n_signal);
  scaler_signal->SetBranchAddress("FCup_run", &FCup_read_run_signal);
  
  for (int iEntry = 0; scaler_signal->LoadTree(iEntry) >= 0; ++iEntry)
    {
      scaler_signal->GetEntry(iEntry);
      if(FCup_read_run_signal==0) {cerr << "At least one signal file has 0 FCup normalization." << endl; return;}
      FCup_p_signal+=FCup_hel_p_signal;
      FCup_n_signal+=FCup_hel_n_signal;
      FCup_run_signal+=FCup_read_run_signal;
    }
  //-------------------------------------------------------------------------Get exclu cuts------------------------------------------------------------------------------------------//
  TCut total = exclu_cuts(target, ErrorBits);
  TCut hel_plus("(Helicity == 1)");
  TCut hel_minus("(Helicity == -1)");
  TCut total_hel_plus  = total + hel_plus ;
  TCut total_hel_minus = total + hel_minus;
    
  //------------------------------------------------------------------------Q2 binning----------------------------------------------------------------------------------------//
  //To be changed with whatever is needed
  //Q2 binning for: - max likelihood estimator
  //                - checking TPol indep of it
//  cout << endl;
//  cout << "Q2 binning" << endl;
  double Q2_low = 1.8;//1.8;//1.2;
  double Q2_high = 3.5;//3.5;//5;
  const int nbins = 10;

  // Create logarithmically spaced bin edges
  double logQ2_low = TMath::Log10(Q2_low);
  double logQ2_high = TMath::Log10(Q2_high);
  double step = (logQ2_high - logQ2_low) / nbins;
    
//  cout << ">>> Q2 from " << Q2_low << " GeV^2 to " << Q2_high << " GeV^2 in " << nbins << " bins" << endl;
//  cout << endl;
  
  //Binned values for asym computation
  double centers_Q2[nbins], centers_theta[nbins];
  double A_in_bins[nbins];
  double err_Q2[nbins];
  //Binned asym and PbPt are a sanity check
  double N_plus[nbins], N_minus[nbins], asym[nbins], PbPt[nbins];
  double N_plus_C[nbins], N_minus_C[nbins], asym_C[nbins];
  double err_N_plus[nbins], err_N_minus[nbins], err_asym[nbins], err_PbPt[nbins];   
  double err_N_plus_C[nbins], err_N_minus_C[nbins], err_asym_C[nbins];
  double Df_vect[nbins], err_Df[nbins];
  int count_bins = -1;
  double Np_C=0, Nm_C=0; //Integrated yields for C
  
  //Max likelihood binned value
  double numerator = 0;
  double numerator_for_error = 0;
  double denominator = 0;
  double denominator_for_error = 0;
  double numerator_error_DF = 0;
  
  for(double i_log = logQ2_low; TMath::Abs(i_log-logQ2_high) > 0.000001; i_log+=step){    
    
    count_bins++;
    
    double i = TMath::Power(10, i_log);
    double i_plus_step = TMath::Power(10, (i_log+step));

    //Binning cut and display
    TCut Q2_bin = ("Q2>=" + to_string(i) + "&&Q2<" + to_string(i_plus_step)).c_str();
    auto Q2_bin_fortext = ("Q^{2} >= " + to_string(i) + " Q^{2} < " + to_string(i_plus_step));

//    cout << "Q2>=" + to_string(i) + "&&Q2<" + to_string(i_plus_step) << endl;

    //Mean Q2 of the bin
    string name_Q2_distrib = ("h_Q2_distrib"+to_string(count_bins));
    auto h_Q2_distrib = new TH1F(name_Q2_distrib.c_str(),"",80,Q2_low,Q2_high);
    string todo_Q2_distrib = "Q2>>" + name_Q2_distrib;
    REC_signal->Draw(todo_Q2_distrib.c_str(),total+Q2_bin,"goff");
    centers_Q2[count_bins] = h_Q2_distrib->GetMean();
    err_Q2[count_bins] = 0;

    //Mean theta_e of the bin
    string name_theta_distrib = ("h_theta_distrib"+to_string(i));         
    auto h_theta_distrib = new TH1F(name_theta_distrib.c_str(),"",80,0,50);                                                                                                                            
    string todo_theta_distrib = "El_Theta>>" + name_theta_distrib;
    REC_signal->Draw(todo_theta_distrib.c_str(),total+Q2_bin,"goff");                                                                                                                             
    centers_theta[count_bins] = h_theta_distrib->GetMean();

    //Theoretical asymmetry. arXiv:0707.1861
    double tau =centers_Q2[count_bins]/(4*Pmass*Pmass);
    double theta_rad = centers_theta[count_bins]*M_PI/180.;
    double tan2 = tan(theta_rad/2)*tan(theta_rad/2);
    double epsilon = 1/(1+2*(1+tau)*tan2);
    double G_E = (1+3.439*tau-1.602*pow(tau,2)+0.068*pow(tau,3))/(1+15.055*tau+48.061*pow(tau,2)+99.304*pow(tau,3)+0.012*pow(tau,4)+8.65*pow(tau,5));
    double G_M = 2.7928473*(1-1.465*tau+1.26*pow(tau,2)+0.262*pow(tau,3))/(1+9.627*tau+11.179*pow(tau,4)+13.245*pow(tau,5));
    double G = G_M/G_E;
    double A_num = 2*tau*G*(Pmass/Ebeam + G*(tau*Pmass/Ebeam + (1+tau)*tan2));
    double A_den = 1 + G*G*tau/epsilon;
    A_in_bins[count_bins] = A_num/A_den;
    //Dilution Factor
    /*
    vector<double> Df_and_errors = dilution_factor(REC_C, REC_signal, total, Q2_bin, FCup_run_C, FCup_run_signal, target);
    Df_vect[count_bins] = Df_and_errors[0];
    err_Df[count_bins] = Df_and_errors[1]; 
    */
    //Signal
    N_plus[count_bins] = REC_signal->GetEntries(total_hel_plus+Q2_bin); 
    //N_plus[count_bins] = REC_signal->GetEntries(); 
    err_N_plus[count_bins] = TMath::Sqrt(N_plus[count_bins]);
    N_minus[count_bins] = REC_signal->GetEntries(total_hel_minus+Q2_bin);
    err_N_minus[count_bins] = TMath::Sqrt(N_minus[count_bins]);

    // Print out raw counts before FCup normalization
//    cout <<"Raw counts before FCup normalization:\n";
//    cout <<"--> N_plus  = "<< N_plus[count_bins]  << endl;
//    cout <<"--> N_minus = "<< N_minus[count_bins] << endl;
    // Write the counts for the run
    fout << centers_Q2[count_bins] <<"  "<< centers_theta[count_bins] <<"  "<< A_in_bins[count_bins]<<"  "<<
	     N_plus[count_bins] <<"  "<< N_minus[count_bins] <<"  "<< FCup_p_signal <<"  "<< FCup_n_signal << endl;

    // Clear out the histograms already in memory...
// h_Q2_distrib0 (Potential memory leak).
// Warning in <TROOT::Append>: Replacing existing TH1: h_theta_distrib1.800000
    h_Q2_distrib->Delete();
    h_theta_distrib->Delete();

  }
  fout.close();

 } // End of loop on "runs" vector

}



