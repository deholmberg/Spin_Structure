/***********************************************************************************
 *
 * Author: Derek Holmberg
 *
 * Date Created: 10/2/2026
 *
 * Last Modified: 10/2/2026
 *
 * Purpose:
 * This program parses the ROOT file skims created from the gmn trains from Noemie's
 * elastic analysis code. Using the FD and CD cuts defined by her code, this
 * program creates TProfile2D objects that are compatible with my data structure,
 * making sure that all the same fiducial/exlusivity cuts Noemie uses are also 
 * applied to my code.
 *
***********************************************************************************/

#include "Input_Functions.h"

using namespace std;

// This calculates the FC asymmetry with an error.
// First entry is the FC asymmetry, the second is the error.
vector<double> FCup_Asymmetry( double Fp, double Fn, double ReadP, double ReadN ){

    vector<double> Asyms = {0,0};

    double d_Fp = ReadP / 906.2; // Uncertainties on each FC charge
    double d_Fn = ReadN / 906.2;

    if( Fp + Fn > 0 && ReadP > 0 && ReadN > 0 ){
	Asyms[0] = (Fp - Fn) / (Fp + Fn);
	Asyms[1] = (2.0/pow((Fp+Fn),2)) * sqrt( Fn*Fn*d_Fp + Fp*Fp*d_Fn );
    }

    return Asyms;

}

// This function gets the Faraday cup charge for a given run
void FCupLoop( string FilePath, string Filter, int runToAnalyze, ofstream& fout ){

    cout << "Opening file for run "<< runToAnalyze << endl;

    ROOT::EnableImplicitMT();

    ROOT::RDataFrame df("FaradayCup", string(FilePath+".root").c_str());
   
    df.GetColumnNames();

    // Apply filters and get the FCup charge for each helicity state
    string helP_Filter = "("+ Filter +")&&( helicity == 1 )";
    auto df_helP = df.Filter( helP_Filter.c_str() );
    double fcup_HelP = *df_helP.Sum("fcupgated");
    auto HelP_Counts = *df_helP.Count();

    string helN_Filter = "("+ Filter +")&&( helicity == -1 )";
    auto df_helN = df.Filter( helN_Filter.c_str() );
    double fcup_HelN = *df_helN.Sum("fcupgated");
    auto HelN_Counts = *df_helN.Count();
/*
    cout << "For run "<< runToAnalyze <<":\n";
    cout << "--> FCup Hel +1 = "<< fcup_HelP <<", "<< HelP_Counts <<" counts\n";
    cout << "--> FCup Hel -1 = "<< fcup_HelN <<", "<< HelN_Counts <<" counts\n";
    vector<double> FC_Asyms = FCup_Asymmetry( fcup_HelP, fcup_HelN, HelP_Counts, HelN_Counts );
    cout << FC_Asyms[0] <<" +- "<< FC_Asyms[1] << endl;
*/
    //vector<double> FCup_Values = { fcup_HelP, fcup_HelN, HelP_Counts, HelN_Counts };
    //return FCup_Values;
    fout << runToAnalyze <<"   "<< fcup_HelP <<"   "<< fcup_HelN <<"   "<< HelP_Counts <<"   "<< HelN_Counts <<"   ";

}

// Loops over a set of input files and writes the histograms to a ROOT file
void FileLoop( string FilePath, string Filter, int runToAnalyze ){

    ROOT::EnableImplicitMT();

    ROOT::RDataFrame df("Particle", string(FilePath+".root").c_str());
    //ROOT::RDataFrame df("Particle", string( FilePath ).c_str());
   
    //cout << df.Count().GetValue() << endl;
    df.GetColumnNames();

    // Define strings that are used in the Define commands to create relevant variables...
    //string nucMass = format( "{:.9f}", nucleon_mass );
    stringstream sm; sm << fixed << setprecision(9) << nucleon_mass;
    string nucMass = sm.str();

    string beamEng; // Declare the energy based on run period
/*
    if( runToAnalyze < 17066 ) beamEng = format( "{:.4f}", beam_energy1 );
    else if( runToAnalyze >= 17067 && runToAnalyze <= 17704 ) beamEng = format( "{:.4f}", beam_energy2 );
    else if( runToAnalyze >= 17720 ) beamEng = format( "{:.4f}", beam_energy3 );
*/
    stringstream s;
    if( runToAnalyze < 17066 ){
	s << fixed << setprecision(4) << beam_energy1;
	beamEng = s.str();	
    }
    else if( runToAnalyze >= 17067 && runToAnalyze <= 17704 ){
	s << fixed << setprecision(4) << beam_energy2;
	beamEng = s.str();
    }
    else if( runToAnalyze >= 17720 ){
	s << fixed << setprecision(4) << beam_energy3;
	beamEng = s.str();
    }

    string Q2_Variable = "2 * "+ beamEng +"* Pt * (1.0 - TMath::Cos( theta * M_PI/180.0 ))";
    string W_Variable  = "TMath::Sqrt("+ nucMass+"*"+nucMass+" + 2 * "+nucMass+"*("+beamEng+" - Pt) - Q2 )";
    string X_Variable  = "Q2 / (2 * "+ nucMass +"*("+beamEng+" - Pt) )";

    //cout << Filter << endl;

    cout <<"Starting input loop on run "<< runToAnalyze <<"...\n";
/*
    cout << nucMass <<"  "<< beamEng << endl;
    cout << Q2_Variable << endl;
    cout << W_Variable << endl;
    cout << X_Variable << endl;
*/
    // Define new branches that weren't saved in the files for space reasons...
    auto df_extra =   df.Define("Pt", "TMath::Sqrt(Px*Px+Py*Py+Pz*Pz)") \
		        .Define("theta", "TMath::ACos(Pz/Pt)*180.0/M_PI") \
		        .Define("phi", "TMath::ATan2(Py,Px)*180.0/M_PI" ) \
		        .Define("SF", "(PcalE + EcinE + EoutE)/Pt" ) \
		        .Define("Einner","(EcinE + EoutE)") \
			.Define("Etotal","(PcalE + EcinE + EoutE)") \
		        .Define("SFpcal","PcalE/Pt") \
			.Define("SFecin","EcinE/Pt") \
			.Define("SFeout","EoutE/Pt") \
			.Define("M2","(M2_PCAL + M2_ECIN + M2_EOUT)/3.0") \
			.Define("Q2", Q2_Variable.c_str() ) \
			.Define("W", W_Variable.c_str() ) \
			.Define("X", X_Variable.c_str() );

    auto df_filtered = df_extra.Filter( Filter.c_str() );

    auto Helicity = df_filtered.Histo1D({"Helicity","Helicity Distribution; Helicity;",8,-2,2},"Helicity");

    // Number of x and Q2 bins
    int nQ2Bins = Q2_Bin_Bounds.size()-1;
    int nXBins  = X_Bin_Bounds.size()-1;

    // This stores the average value of Q2 for each of the x, Q2 bins
    auto Q2BinData = df_filtered.Profile2D({"Q2BinData","Q2 Bin Data; Q2 [GeV^{2}]; X", nQ2Bins, Q2_Bin_Bounds.data(), nXBins, X_Bin_Bounds.data()},"Q2","X","Q2");
    // This does the same but for the X data
    auto XBinData = df_filtered.Profile2D({"XBinData","X Bin Data; Q2 [GeV^{2}]; X", nQ2Bins, Q2_Bin_Bounds.data(), nXBins, X_Bin_Bounds.data()},"Q2","X","X");

    // Now create the same bins as above but for each beam helicity state
    // Helicity +1
    auto df_helPlus = df_filtered.Filter("Helicity == 1");
    auto Q2BinData_HelPlus = df_helPlus.Profile2D({"Q2BinData_HelPlus","Q2 Bin Data; Q2 [GeV^{2}]; X", nQ2Bins, Q2_Bin_Bounds.data(), nXBins, X_Bin_Bounds.data()},"Q2","X","Q2");
    auto XBinData_HelPlus = df_helPlus.Profile2D({"XBinData_HelPlus"  ,"X Bin Data; Q2 [GeV^{2}]; X", nQ2Bins, Q2_Bin_Bounds.data(), nXBins, X_Bin_Bounds.data()},"Q2","X","X");
    // Helicity -1
    auto df_helMinus = df_filtered.Filter("Helicity == -1");
    auto Q2BinData_HelMinus = df_helMinus.Profile2D({"Q2BinData_HelMinus","Q2 Bin Data; Q2 [GeV^{2}]; X",nQ2Bins, Q2_Bin_Bounds.data(), nXBins, X_Bin_Bounds.data()},"Q2","X","Q2");
    auto XBinData_HelMinus = df_helMinus.Profile2D({"XBinData_HelMinus"  ,"X Bin Data; Q2 [GeV^{2}]; X", nQ2Bins, Q2_Bin_Bounds.data(), nXBins, X_Bin_Bounds.data()},"Q2","X","X");

/*
    auto Theta    = df_filtered.Histo1D({"Theta","Theta Distribution; #theta [deg];",450,0,45},"theta");
    auto Phi      = df_filtered.Histo1D({"Phi","Phi Distribution; #phi [deg];",3600,-180,180},"phi");
    auto Nphe     = df_filtered.Histo1D({"nphe","Photoelectrons; Nphe;;",50,0,50},"nphe");
    auto Vx       = df_filtered.Histo1D({"Vx","X-Vertex Distribution; V_{X} [cm];",2000,-15,5},"Vx");
    auto Vy       = df_filtered.Histo1D({"Vy","Y-Vertex Distribution; V_{Y} [cm];",2000,-15,5},"Vy");
    auto Vz       = df_filtered.Histo1D({"Vz","Z-Vertex Distribution; V_{Z} [cm];",2000,-15,5},"Vz");
    auto Beta     = df_filtered.Histo1D({"Beta","#beta Distribution; #beta;",1200,-100,2},"Beta");
    auto Q2	  = df_filtered.Histo1D({"Q2","Q^{2} Distribution; Q^{2} [GeV^{2}];",1100,0,11},"Q2");
    auto W	  = df_filtered.Histo1D({"W","W Distribution; W [GeV];",500,0,5},"W");
    auto X	  = df_filtered.Histo1D({"X","X Distribution; X;",100,0,1},"X");
    auto M2_Dist  = df_filtered.Histo1D({"M2_Dist","ECAL Second Moments Distribution; M2 [cm^{2}];",4500,0,450},"M2");
    auto PID      = df_filtered.Histo1D({"pid","PID Distribution; PID;",6000,-3000,3000},"pid");
  
    auto M2_PCAL  = df_filtered.Histo1D({"M2_PCAL","Second Moment of PCAL; M2_{PCAL};",450,0,450},"M2_PCAL");
    auto M2_ECIN  = df_filtered.Histo1D({"M2_ECIN","Second Moment of ECIN; M2_{ECIN};",450,0,450},"M2_ECIN");
    auto M2_EOUT  = df_filtered.Histo1D({"M2_EOUT","Second Moment of EOUT; M2_{EOUT};",450,0,450},"M2_EOUT");

    auto SFpcal   = df_filtered.Histo1D({"SFpcal","Partial SF for PCAL; E_{PCAL}/P;",50,0,0.25},"SFpcal");
    auto SFecin   = df_filtered.Histo1D({"SFecin","Partial SF for ECIN; E_{ECIN}/P;",50,0,0.25},"SFecin");
    auto SFeout   = df_filtered.Histo1D({"SFeout","Partial SF for EOUT; E_{EOUT}/P;",50,0,0.25},"SFeout");

    auto DC_R1 = df_filtered.Histo2D({"DC_R1","DC Reduced #chi^{2}/NDF vs Edge Region 1; Edge [cm]; #chi^{2}/NDF",500,0,50,100,0,10},"EdgeR1","TrackChi2NDF");
    auto DC_R2 = df_filtered.Histo2D({"DC_R2","DC Reduced #chi^{2}/NDF vs Edge Region 2; Edge [cm]; #chi^{2}/NDF",800,0,80,100,0,10},"EdgeR2","TrackChi2NDF");
    auto DC_R3 = df_filtered.Histo2D({"DC_R3","DC Reduced #chi^{2}/NDF vs Edge Region 3; Edge [cm]; #chi^{2}/NDF",1000,0,100,100,0,10},"EdgeR3","TrackChi2NDF");

    auto SFvsP       = df_filtered.Histo2D({"SFvsP","Sampling Fraction vs. Momentum; P [GeV]; SF",110,0,11,40,0,0.4},"Pt","SF");
    auto PcalVsSF    = df_filtered.Histo2D({"PcalVsSF","E_{PCAL} vs. SF; SF; E_{PCAL} [GeV]",40,0,0.4,200,0,2},"SF","PcalE");
    auto EtotalVsSF  = df_filtered.Histo2D({"EtotalVsSF","E_{Total} vs. SF; SF; E_{PCAL}+E_{ECIN}+E_{EOUT} [GeV]",40,0,0.4,400,0,4},"SF","Etotal");
    auto PcalEcin    = df_filtered.Histo2D({"PcalEcin","E_{PCAL} vs. E_{ECIN}; E_{ECIN} [GeV]; E_{PCAL} [GeV]",200,0,2,200,0,2},"EcinE","PcalE");
    auto PartialSF   = df_filtered.Histo2D({"PartialSF","E_{PCAL}/P vs. E_{ECIN}/P; E_{ECIN}/P; E_{PCAL}/P",40,0,0.4,40,0,0.4},"SFecin","SFpcal");
    auto SF1D	     = df_filtered.Histo1D({"SF1D","Sampling Fraction Distribution; SF;",400,0,0.4},"SF");
    auto NpheVsP     = df_filtered.Histo2D({"NpheVsP","Photoelectrons vs. Momentum; P [GeV]; Nphe",110,0,11,50,0,50},"Pt","nphe");
    auto BetaVsP     = df_filtered.Histo2D({"BetaVsP","#beta vs. Momentum; P [GeV]; #beta",110,0,11,110,0,1.1},"Pt","Beta");
    auto SFvsM2      = df_filtered.Histo2D({"SFvsM2","Sampling Fraction vs. M2; M2 [cm^{2}]; SF",450,0,450,40,0,0.4},"M2","SF");
    auto NphevsM2    = df_filtered.Histo2D({"NphevsM2","Photoelectrons vs. M2; M2 [cm^{2}]; Nphe",450,0,450,50,0,50},"M2","nphe");   
*/
    //auto PcalSFvsPcalLv  = df_filtered.Histo2D({"PcalPcalLv","PCAL SF vs. PCAL Lv; Lv [cm]; SF",100,0,450,40,0,0.4},"PcalLv","SF");
    //auto PcalSFvsPcalLw  = df_filtered.Histo2D({"PcalPcalLw","PCAL SF vs. PCAL Lw; Lw [cm]; SF",100,0,450,40,0,0.4},"PcalLw","SF");

/*
    // Save a snapshot with the newly created variables for the second moments (M2's) and the partial sampling fractions
    if( snapName != "None" ){
        vector<string> SnapshotVars = {"M2_PCAL","M2_ECIN","M2_EOUT","SFpcal","SFecin","SFeout","nphe","Pt","theta"};
        string snapPath = "ParticleTrees/Snapshots/"+snapName+"_Snapshot.root";
        df_filtered.Snapshot( "Particle", snapPath.c_str(), SnapshotVars );
    }

    // Create a Profile2D object for calculating the average values of Pt and theta within
    // each theta and Pt bin...
    vector<double> ThetaBinBounds = {8, 11, 14, 17, 20, 23, 26, 29, 32, 40};
    int nThetaBins = ThetaBinBounds.size()-1;

    int nMomentumBins = 110;
    double pLow = 0; double pBig = 11;

    auto MomentumBinData = df_filtered.Profile2D({"MomentumBinData","Momentum Bin Data; P [GeV]; #theta [deg]", nMomentumBins, pLow, pBig, nThetaBins, ThetaBinBounds.data()},"Pt","theta","Pt");
    auto ThetaBinData    = df_filtered.Profile2D({"ThetaBinData"   ,"#theta Bin Data; P [GeV]; #theta [deg]"  , nMomentumBins, pLow, pBig, nThetaBins, ThetaBinBounds.data()},"Pt","theta","theta");

    string sector1Filter = "("+ Filter +") && SectorECAL == 1";
    auto df_sector1 = df_filtered.Filter( sector1Filter.c_str() );
    auto PcalSF_S1   = df_sector1.Histo2D({"PcalSF_S1","PCAL SF vs. Lv Sector 1; Lv [cm]; SF",450,0,450,40,0,0.4},"PcalLv","SF");

    string sector2Filter = "("+ Filter +") && SectorECAL == 2";
    auto df_sector2 = df_filtered.Filter( sector2Filter.c_str() );
    auto PcalSF_S2   = df_sector2.Histo2D({"PcalSF_S2","PCAL SF vs. Lv Sector 2; Lv [cm]; SF",450,0,450,40,0,0.4},"PcalLv","SF");
*/


    Helicity->Write();
    Q2BinData->Write();
    XBinData->Write();
    Q2BinData_HelPlus->Write();
    XBinData_HelPlus->Write();
    Q2BinData_HelMinus->Write();
    XBinData_HelMinus->Write();

/*
    Nphe->Write();
    Theta->Write();
    Phi->Write();
    SFvsP->Write();
    PcalVsSF->Write();
    EtotalVsSF->Write();
    PcalEcin->Write();
    PartialSF->Write();
    NpheVsP->Write();
    DC_R1->Write();
    DC_R2->Write();
    DC_R3->Write();
    Vx->Write();
    Vy->Write();
    Vz->Write();
    Beta->Write();
    PID->Write();

    Q2->Write();
    X->Write();
    W->Write();

    BetaVsP->Write();
    SFpcal->Write();
    SFecin->Write();
    SFeout->Write();

    M2_PCAL->Write();
    M2_ECIN->Write();
    M2_EOUT->Write();
    M2_Dist->Write();
    SFvsM2->Write();
    NphevsM2->Write();

    PcalSF_S1->Write();
    PcalSF_S2->Write();
*/
    cout <<"Finished loop on run "<< runToAnalyze <<"...\n";
}

void Data_Loop(string runperiod){

/* Things to add:
 * Run-dependent PCAL fiducial cut
 * Run-dependent vertex cut
 * Run-dependent SF cut
 */

  // Start a loop over all runs:
  RunPeriod Period;

  // Output file for all the fcup values
  ofstream fout( string("LatestSkim/Text_Files/"+ runperiod +"_FCup_Data.txt") );
  
  // Get all good runs to read in:
  //vector<int> All_Runs = Period.getAllGoodRuns();
  vector<int> All_Runs = Period.getDataSetGoodRuns( runperiod );

  //vector<int> All_Runs = { 16137 };//, 17796, 17797 };
 
  for( int Run : All_Runs ){

    string DC_Cuts     = "EdgeR1 > 4 && EdgeR2 > 5 && EdgeR3 > 8";
    //string PCAL_Cuts   = "(SectorECAL == 1 && (PcalLv > 22.5 && PcalLw > 22.5)) || ( (SectorECAL >= 2) && (PcalLv > 13.5 && PcalLw > 13.5) )";
    //string PCAL_Cuts   = " (SectorECAL >= 1) && (PcalLv > 13.5 && PcalLw > 13.5) ";
    string VtxCut_Su22 = "Vz > -9   && Vz < 1   && Vx < 2 && Vx > -2 && Vy < 2 && Vy > -2"; // Vertex cuts for the summer
    string VtxCut_FaSp = "Vz > -7.5 && Vz < 2.5 && Vx < 2 && Vx > -2 && Vy < 2 && Vy > -2"; // Vertex cuts for the fall and spring

    string Energy_Cuts = "PcalE > 0.06 && EcinE > 0.01";
    string Nphe_LowCut = "nphe < 2";
    string Nphe_BigCut = "nphe > 7";
    string Nphe_AvgCut = "nphe > 2";
    string Beta_Cut    = "Beta > 0.997 && Beta < 1.003";
    string Kinematic_Cuts = "Q2 > 2.2 && Pt > 2.6 && theta > 7 && theta < 40 && W > 2";

    string SFCut_Su22_S1 = "(SectorECAL == 1 && SF > (0.190837 + 0.005999*Pt - 0.000627*Pt*Pt) && SF < (0.303454 - 0.005738*Pt + 0.000352*Pt*Pt) )";
    string SFCut_Su22_S2 = "(SectorECAL == 2 && SF > (0.187330 + 0.006956*Pt - 0.000670*Pt*Pt) && SF < (0.305030 - 0.000522*Pt - 0.000462*Pt*Pt) )";
    string SFCut_Su22_S3 = "(SectorECAL == 3 && SF > (0.186403 + 0.008207*Pt - 0.000806*Pt*Pt) && SF < (0.304494 + 0.001433*Pt - 0.000819*Pt*Pt) )";
    string SFCut_Su22_S4 = "(SectorECAL == 4 && SF > (0.176055 + 0.011733*Pt - 0.001074*Pt*Pt) && SF < (0.309051 - 0.001400*Pt - 0.000397*Pt*Pt) )";
    string SFCut_Su22_S5 = "(SectorECAL == 5 && SF > (0.178196 + 0.009725*Pt - 0.000872*Pt*Pt) && SF < (0.311447 - 0.004086*Pt - 0.000196*Pt*Pt) )";
    string SFCut_Su22_S6 = "(SectorECAL == 6 && SF > (0.183186 + 0.008936*Pt - 0.000870*Pt*Pt) && SF < (0.306155 - 0.000247*Pt - 0.000602*Pt*Pt) )";
    string Total_Su22_SF_Cut = "("+ SFCut_Su22_S1 +")||("+SFCut_Su22_S2 +")||("+ SFCut_Su22_S3 +")||("+ SFCut_Su22_S4 +")||("+ SFCut_Su22_S5 +")||("+ SFCut_Su22_S6 +")";
	    
    string SFCut_Fa22_S1 = "(SectorECAL == 1 && SF > (0.173820 + 0.011899*Pt - 0.001259*Pt*Pt) && SF < (0.322944 - 0.012585*Pt + 0.001100*Pt*Pt) )";
    string SFCut_Fa22_S2 = "(SectorECAL == 2 && SF > (0.180660 + 0.008456*Pt - 0.000790*Pt*Pt) && SF < (0.311690 - 0.004250*Pt - 0.000068*Pt*Pt) )";
    string SFCut_Fa22_S3 = "(SectorECAL == 3 && SF > (0.190434 + 0.005721*Pt - 0.000539*Pt*Pt) && SF < (0.311784 - 0.003025*Pt - 0.000286*Pt*Pt) )";
    string SFCut_Fa22_S4 = "(SectorECAL == 4 && SF > (0.182037 + 0.008472*Pt - 0.000775*Pt*Pt) && SF < (0.316576 - 0.003039*Pt - 0.000191*Pt*Pt) )";
    string SFCut_Fa22_S5 = "(SectorECAL == 5 && SF > (0.181404 + 0.007415*Pt - 0.000624*Pt*Pt) && SF < (0.313592 - 0.005845*Pt + 0.000040*Pt*Pt) )";
    string SFCut_Fa22_S6 = "(SectorECAL == 6 && SF > (0.179893 + 0.009279*Pt - 0.000846*Pt*Pt) && SF < (0.318022 - 0.007117*Pt + 0.000108*Pt*Pt) )";
    string Total_Fa22_SF_Cut = "("+ SFCut_Fa22_S1 +")||("+SFCut_Fa22_S2 +")||("+ SFCut_Fa22_S3 +")||("+ SFCut_Fa22_S4 +")||("+ SFCut_Fa22_S5 +")||("+ SFCut_Fa22_S6 +")";
    
    string SFCut_Sp23_S1 = "(SectorECAL == 1 && SF > (0.188243 + 0.004928*Pt - 0.000498*Pt*Pt) && SF < (0.308340 - 0.007150*Pt + 0.000487*Pt*Pt) )";
    string SFCut_Sp23_S2 = "(SectorECAL == 2 && SF > (0.183701 + 0.008171*Pt - 0.000765*Pt*Pt) && SF < (0.310391 - 0.003058*Pt - 0.000184*Pt*Pt) )";
    string SFCut_Sp23_S3 = "(SectorECAL == 3 && SF > (0.185788 + 0.007974*Pt - 0.000773*Pt*Pt) && SF < (0.310277 - 0.004496*Pt - 0.000111*Pt*Pt) )";
    string SFCut_Sp23_S4 = "(SectorECAL == 4 && SF > (0.177398 + 0.010710*Pt - 0.000954*Pt*Pt) && SF < (0.314889 - 0.003913*Pt - 0.000101*Pt*Pt) )";
    string SFCut_Sp23_S5 = "(SectorECAL == 5 && SF > (0.181192 + 0.007219*Pt - 0.000601*Pt*Pt) && SF < (0.312211 - 0.007023*Pt + 0.000140*Pt*Pt) )";
    string SFCut_Sp23_S6 = "(SectorECAL == 6 && SF > (0.187191 + 0.005478*Pt - 0.000496*Pt*Pt) && SF < (0.309825 - 0.004853*Pt - 0.000123*Pt*Pt) )";
    string Total_Sp23_SF_Cut = "("+ SFCut_Sp23_S1 +")||("+SFCut_Sp23_S2 +")||("+ SFCut_Sp23_S3 +")||("+ SFCut_Sp23_S4 +")||("+ SFCut_Sp23_S5 +")||("+ SFCut_Sp23_S6 +")";

    // Apply run-dependent cuts
    string PCAL_Cuts, SF_Cuts, Vertex_Cuts, inputFile, fcupFile;

    string Target = Period.getTargetType( Run );
    if( Target == "Empty" or Target == "Foil" ) Target = "ET";

    // Summer Data
    if( Run < 16156 ){ // Summer data before the ECAL sector 1 issue
	PCAL_Cuts   = " (SectorECAL >= 1) && (PcalLv > 13.5 && PcalLw > 13.5) ";
	SF_Cuts     = Total_Su22_SF_Cut;
	Vertex_Cuts = VtxCut_Su22;
	inputFile   = "/volatile/clas12/holmberg/Summer_22_Data/Electron_"+ Target +"_"+ to_string( Run ) +"_Data";
	fcupFile    = "../FCup_Analysis/Summer/Electron_"+ Target +"_"+ to_string( Run ) +"_FCup";
    }
    else if( Run >= 16156 && Run <= 16772 ){ // Summer data after the ECAL sector 1 issue appeared
	PCAL_Cuts   = "(SectorECAL == 1 && (PcalLv > 22.5 && PcalLw > 22.5)) || ( (SectorECAL >= 2) && (PcalLv > 13.5 && PcalLw > 13.5) )";
	SF_Cuts     = Total_Su22_SF_Cut;
	Vertex_Cuts = VtxCut_Su22;
	inputFile   = "/volatile/clas12/holmberg/Summer_22_Data/Electron_"+ Target +"_"+ to_string( Run ) +"_Data";
	fcupFile    = "../FCup_Analysis/Summer/Electron_"+ Target +"_"+ to_string( Run ) +"_FCup";
    }
    else if( Run > 16800 && Run <= 17183 ){ // Fall data (negative solenoid for run below and including 17183; pos. sol. for runs above and including 17188)
	PCAL_Cuts   = "(SectorECAL == 1 && (PcalLv > 22.5 && PcalLw > 22.5)) || ( (SectorECAL >= 2) && (PcalLv > 13.5 && PcalLw > 13.5) )";
	SF_Cuts     = Total_Fa22_SF_Cut;
	Vertex_Cuts = VtxCut_FaSp;
	inputFile   = "/volatile/clas12/holmberg/Fall_22_Data/Electron_"+ Target +"_"+ to_string( Run ) +"_Data";
	fcupFile    = "../FCup_Analysis/FallNeg/Electron_"+ Target +"_"+ to_string( Run ) +"_FCup";
    }
    else if( Run > 17185 && Run <= 17408 ){ // Fall data (negative solenoid for run below and including 17183; pos. sol. for runs above and including 17188)
	PCAL_Cuts   = "(SectorECAL == 1 && (PcalLv > 22.5 && PcalLw > 22.5)) || ( (SectorECAL >= 2) && (PcalLv > 13.5 && PcalLw > 13.5) )";
	SF_Cuts     = Total_Fa22_SF_Cut;
	Vertex_Cuts = VtxCut_FaSp;
	inputFile   = "/volatile/clas12/holmberg/Fall_22_Data/Electron_"+ Target +"_"+ to_string( Run ) +"_Data";
	fcupFile    = "../FCup_Analysis/FallPos/Electron_"+ Target +"_"+ to_string( Run ) +"_FCup";
    }
    else if( Run > 17450 && Run <= 17768 ){ // Spring inbending
	PCAL_Cuts   = "(SectorECAL == 1 && (PcalLv > 22.5 && PcalLw > 22.5)) || ( (SectorECAL >= 2) && (PcalLv > 13.5 && PcalLw > 13.5) )";
	SF_Cuts     = Total_Sp23_SF_Cut;
	Vertex_Cuts = VtxCut_FaSp;
	inputFile   = "/volatile/clas12/holmberg/Spring_23_Data/Electron_"+ Target +"_"+ to_string( Run ) +"_Data";
	fcupFile    = "../FCup_Analysis/SpringInb/Electron_"+ Target +"_"+ to_string( Run ) +"_FCup";
    }
    else if( Run >= 17769 ){ // Spring outbending
	PCAL_Cuts   = " (SectorECAL >= 1) && ( PcalLv > 18 && PcalLw > 18 ) ";
	SF_Cuts     = Total_Sp23_SF_Cut;
	Vertex_Cuts = VtxCut_Su22; // Reuses the summer vertex cuts
	inputFile   = "/volatile/clas12/holmberg/Spring_23_Data/Electron_"+ Target +"_"+ to_string( Run ) +"_Data";
	fcupFile    = "../FCup_Analysis/SpringOutb/Electron_"+ Target +"_"+ to_string( Run ) +"_FCup";
    }

    //string TotalCutAllSu22 = "("+ DC_Cuts +")&&("+ PCAL_Cuts +")&&("+ Vertex_Cuts + ")&&("+ Total_Su22_SF_Cut +")&&("+ Kinematic_Cuts +")&&( pid == 11 )";

    string TotalCut = "("+ DC_Cuts +")&&("+ PCAL_Cuts +")&&("+ Vertex_Cuts + ")&&("+ SF_Cuts +")&&("+ Kinematic_Cuts +")&&( pid == 11 )&&( Helicity != 0 )";

    string outFileName = "LatestSkim/" + Period.getTargetType( Run ) +"_"+ to_string( Run ) +"_Data.root";

    TFile* outFileTest = new TFile( outFileName.c_str() , "RECREATE" );
    FileLoop( inputFile, TotalCut, Run );
    outFileTest->Close();

    FCupLoop( fcupFile, "(fcupgated > 0) && (clockgated > 1000)", Run, fout );
    FCupLoop( fcupFile, "(clockgated > -100000)", Run, fout );
    fout << endl;

    //cout << Run << endl;

  } // End of file loop

  fout.close();

}
