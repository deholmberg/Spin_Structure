/***********************************************************************************
 *
 * Author: Derek Holmberg
 *
 * Date Created: 8/24/2026
 *
 * Last Modified: 9/18/2026
 *
 * Purpose:
 * The purpose of this program is to create TTrees of all available runs
 * from SIDIS skims in the RGC data set, looking at only trigger electrons.
 * This program creates root files in my /volatile/ directory, along with
 * text files tracking the integrated FC charge. 
 *
***********************************************************************************/

#include "/w/hallb-scshelf2102/clas12/holmberg/Spin_Structure/Header_Files/Binning_Classes.h" 
#include "/w/hallb-scshelf2102/clas12/holmberg/Spin_Structure/Header_Files/Functions_RGC.h"
#include "QADB.h"

using namespace std;
using namespace QA;


// This is the main loop that goes through the data
void GetPairTree(){

   string path1, path2, thisDir, outDir, fileExt;
   double percentCut; // Percentage of the raster to keep
   int runToAnalyze; // The starting run to look at
   int PID; // The particle ID to look at
   //cin >> path_to_sidis >> TARGET_TYPE >> runStart >> runEnd >> thisDir;
   cin >> path1 >> path2 >> runToAnalyze >> PID >> fileExt;

   // The following cuts off the "/slurm" part of the string
   int event_counter = 0; // Counts total number of events

   // These are used to subtract off from the "AccumulateCharge" QADB functions so the total
   // charges can be found on a run by run basis
   double totalRUNcharge = 0.0;
   double totalHLMcharge = 0.0; // Charge for HEL +1 states (Tpol sign-corrected)
   double totalHLPcharge = 0.0; // Charge for HEL -1 states (Tpol sign-corrected)
   double totalHL0charge = 0.0; // Charge for HEL  0 states
   double totalHLMug_charge = 0.0; // Ungated charge for HEL +1 states (Tpol sign-corrected)
   double totalHLPug_charge = 0.0; // Ungated charge for HEL -1 states (Tpol sign-corrected)
   double totalHL0ug_charge = 0.0; // Ungated charge for HEL  0 states

   // Vector that holds all the possible target configurations
   vector<string> Targets = {"ND3","NH3","C","CH2","ET","CD2"};

   // Read in the information for the runs
   RunPeriod Period; //Period.SetRunPeriod();
   cout << "read in run period.\n";

   // instantiate QADB
   QADB* qa = new QADB("latest");

   // custom QA cut definition
   qa->SetMaskBit("TotalOutlier",false);
   qa->SetMaskBit("TerminalOutlier",false);
   qa->SetMaskBit("MarginalOutlier",false);
   qa->SetMaskBit("SectorLoss",false); 
   qa->SetMaskBit("LowLiveTime",false);
   qa->SetMaskBit("TotalOutlierFT",false);
   qa->SetMaskBit("TerminalOutlierFT",false);
   qa->SetMaskBit("MarginalOutlierFT",false);
   qa->SetMaskBit("LossFT",false);
   qa->SetMaskBit("BSAWrong",false);
   //qa->SetMaskBit("BSAUnknown",false);
   qa->SetMaskBit("ChargeHigh",false);
   qa->SetMaskBit("ChargeNegative",false);
   qa->SetMaskBit("ChargeUnknown",false);
   qa->SetMaskBit("PossiblyNoBeam",false);

 //for(int runToAnalyze = runStart; runToAnalyze <= runEnd; runToAnalyze++){
 //for( string TARGET_TYPE : Targets ){
 
  if( Period.isGoodRun( runToAnalyze ) ){
  
   //qa->SetMaskBit("Misc", Period.applyMiscBit(runToAnalyze) );

   if( runToAnalyze == 17482 ){
	qa->SetMaskBit("TotalOutlier",false);
	qa->SetMaskBit("TotalOutlierFT",false);
   }

   // Set the run period;
   string DataPeriod;
   if(      runToAnalyze > 16100 && runToAnalyze < 16800 ) DataPeriod = "Su22";
   else if( runToAnalyze > 16800 && runToAnalyze < 17450 ) DataPeriod = "Fa22";
   else if( runToAnalyze > 17480 && runToAnalyze < 17812 ) DataPeriod = "Sp23";

   // Set the target type
   string TitlePID;
   if( PID == 11 ) TitlePID = "e^{-}";
   else if( PID == -11 ) TitlePID = "e^{+}";
   else if( PID == 211 ) TitlePID = "#pi^{+}";
   else if( PID == -211) TitlePID = "#pi^{-}";
   else if( PID == 2212) TitlePID = "P";
   else if( PID == 1 ) TitlePID = "Positive_Track";
   else if( PID ==-1 ) TitlePID = "Negative_Track";

   cout << "Get target polarization.\n";
   double thisTPol = Period.getTargetPolarization( runToAnalyze );
   string TARGET_TYPE = Period.getTargetType( runToAnalyze );

   if( TARGET_TYPE == "Empty" || TARGET_TYPE == "Foil" ) TARGET_TYPE = "ET";

   // This is a correction factor for the target polarization; right now, it's set to 1 since it's unused...
   double corr_factor = 1;

   // Gated FCup counts from scalers
   double FCup_hel_n = 0; // gated counts
   double FCup_hel_p = 0;
   double FCup_hel_0 = 0;

   double FCup_ughel_n = 0; // un-gated counts
   double FCup_ughel_p = 0;
   double FCup_ughel_0 = 0;

   double FCup_run = 0; // RUN::scaler value

   //string fullFilePath = path1 + TARGET_TYPE + path2 + to_string(runToAnalyze) + "/*.hipo";
   string fullFilePath = path1 + TARGET_TYPE + path2 + to_string(runToAnalyze) + ".hipo";

   cout << fullFilePath << endl;
 
   TChain *chain = new TChain("","");
   chain->Add( fullFilePath.c_str() );

   hipo::reader  reader;


    // Create the ROOT file for storing the diagnostic plots
    cout << "Create the root file.\n";
    string outPID;
    if( PID == 11 ) outPID = "Electron";
    else if( PID == -11 ) outPID = "Positron";
    else if( PID == 1 ) outPID = "Positive_Track";
    else if( PID ==-1 ) outPID = "Negative_Track";

    //string tfile_path = thisDir+"/"+outDir+"/ROOT_Files/"+outPID+"_"+TARGET_TYPE+"_"+to_string(runToAnalyze)+"_"+fileExt+"_EVIO.root";

    string tfile_path = "/volatile/clas12/holmberg/";
    //string tfile_path = "/w/hallb-scshelf2102/clas12/holmberg/Dilution_Factors/";

    if( runToAnalyze < 16800 ) tfile_path += "Summer_22_NoQA/"+ outPID+"_"+TARGET_TYPE+"_"+to_string(runToAnalyze)+"_"+fileExt+".root";
    else if( runToAnalyze > 16800 && runToAnalyze < 17410 ) tfile_path += "Fall_22_NoQA/"+ outPID+"_"+TARGET_TYPE+"_"+to_string(runToAnalyze)+"_"+fileExt+".root";
    else if( runToAnalyze > 17450 ) tfile_path += "Spring_23_NoQA/"+ outPID+"_"+TARGET_TYPE+"_"+to_string(runToAnalyze)+"_"+fileExt+".root";

    TFile* file = new TFile(tfile_path.c_str(), "RECREATE");

   //Reading the input
   TObjArray *files = chain->GetListOfFiles();
   cout << "size of files   " << files->GetLast() + 1 << endl;

    TTree* tree = new TTree("Particle","Particle");

    // Variables to store in the particle tree
    int pid = 0; // Monte-carlo ID of particle
    int Status = 0; // Location of hit within CLAS12 geometry
    double Chi2PID = 0; // How closely PID matches time-of-flight information based on particle momentum
    int Helicity = 0; // HWP and target polarization-corrected helicity state of particle
    int Charge = 0; // Charge of the particle
    double Px = 0; double Py = 0; double Pz = 0; // Particle momenta
    double Vx = 0; double Vy = 0; double Vz = 0; // Z-vertex location in CLAS12 geometry
    double TrackChi2NDF = 0; // Quality of track reconstruction in the DC
    double EdgeR1 = 0; double EdgeR2 = 0; double EdgeR3 = 0; // Edge variable in each of the three regions of DC
    int SectorDC = 0; int SectorECAL = 0; int SectorHTCC = 0; // Location of track/hit in DC, ECAL, and HTCC, respectively (should be equal)
    double PcalX = 0; double PcalY = 0; // Hit locations in the Pcal
    double PcalE = 0; // Energy deposited in Pcal
    double PcalLu = 0; double PcalLv = 0; double PcalLw = 0; // Lu, Lv, Lw hit positions in the Pcal
    double EcinX = 0; double EcinY = 0; // Hit locations in the Ecin
    double EcinE = 0; // Energy deposited in Ecin
    double EcinLu = 0; double EcinLv = 0; double EcinLw = 0; // Lu, Lv, Lw hit positions in the Ecin
    double EoutX = 0; double EoutY = 0; // Hit locations in the Eout
    double EoutE = 0; // Energy deposited in Eout
    double EoutLu = 0; double EoutLv = 0; double EoutLw = 0; // Lu, Lv, Lw hit positions in the Eout
    double Edep = 0; // Total energy deposited into the ECAL system
    double nphe = 0; // Number of photoelectrons produced in the HTCC
    double Beta = 0; // Velocity divided by light speed (taken from REC::Particle bank)
    double RasterX = 0; // Raster x-position
    double RasterY = 0; // Raster y-position
    double M2_PCAL = 0; // Second moment of the shower in the PCAL
    double M2_ECIN = 0; // Second moment of the shower in the ECIN
    double M2_EOUT = 0; // Second moment of the shower in the EOUT

    // Create branches on the tree to store all of the variables
    tree->Branch("pid", &pid, "pid/I");
    tree->Branch("Status", &Status, "Status/I");
    tree->Branch("Chi2PID", &Chi2PID, "Chi2PID/D");
    tree->Branch("Helicity", &Helicity, "Helicity/I");
    tree->Branch("Charge", &Charge, "Charge/I");
    tree->Branch("SectorDC", &SectorDC, "SectorDC/I");
    tree->Branch("SectorECAL", &SectorECAL, "SectorECAL/I");
    tree->Branch("SectorHTCC", &SectorHTCC, "SectorHTCC/I");
    tree->Branch("Px", &Px, "Px/D"); tree->Branch("Py", &Py, "Py/D"); tree->Branch("Pz", &Pz, "Pz/D");
    tree->Branch("Vx", &Vx, "Vx/D"); tree->Branch("Vy", &Vy, "Vy/D"); tree->Branch("Vz", &Vz, "Vz/D");
    tree->Branch("TrackChi2NDF", &TrackChi2NDF, "TrackChi2NDF/D");
    tree->Branch("EdgeR1", &EdgeR1, "EdgeR1/D"); tree->Branch("EdgeR2", &EdgeR2, "EdgeR2/D"); tree->Branch("EdgeR3", &EdgeR3, "EdgeR3/D");

    tree->Branch("PcalX", &PcalX, "PcalX/D"); tree->Branch("PcalY", &PcalY, "PcalY/D"); tree->Branch("PcalE", &PcalE, "PcalE/D");
    tree->Branch("PcalLu", &PcalLu, "PcalLu/D"); tree->Branch("PcalLv", &PcalLv, "PcalLv/D"); tree->Branch("PcalLw", &PcalLw, "PcalLw/D");
    tree->Branch("EcinX", &EcinX, "EcinX/D"); tree->Branch("EcinY", &EcinY, "EcinY/D"); tree->Branch("EcinE", &EcinE, "EcinE/D");
    tree->Branch("EcinLu", &EcinLu, "EcinLu/D"); tree->Branch("EcinLv", &EcinLv, "EcinLv/D"); tree->Branch("EcinLw", &EcinLw, "EcinLw/D");
    tree->Branch("EoutX", &EoutX, "EoutX/D"); tree->Branch("EoutY", &EoutY, "EoutY/D"); tree->Branch("EoutE", &EoutE, "EoutE/D");
    tree->Branch("EoutLu", &EoutLu, "EoutLu/D"); tree->Branch("EoutLv", &EoutLv, "EoutLv/D"); tree->Branch("EoutLw", &EoutLw, "EoutLw/D");

    tree->Branch("nphe", &nphe, "nphe/D");
    tree->Branch("Beta", &Beta, "Beta/D");

    tree->Branch("RasterX", &RasterX, "RasterX/D");
    tree->Branch("RasterY", &RasterY, "RasterY/D");

    tree->Branch("M2_PCAL", &M2_PCAL, "M2_PCAL/D");
    tree->Branch("M2_ECIN", &M2_ECIN, "M2_ECIN/D");
    tree->Branch("M2_EOUT", &M2_EOUT, "M2_EOUT/D");

    //Start loop on files
    for (int noffiles = 0; noffiles < files->GetEntries(); noffiles++){

        reader.open(files->At(noffiles)->GetTitle());
	cout << files->At(noffiles)->GetTitle() << endl;

        // This is the main analysis loop
        hipo::dictionary  factory;
        reader.readDictionary(factory);
        hipo::structure  particles;
        hipo::structure  detectors;
        hipo::event      event;
        hipo::bank  dataPART;

        hipo::bank PART(factory.getSchema("REC::Particle"));
        hipo::bank HEL(factory.getSchema("REC::Event"));
        hipo::bank TRACK(factory.getSchema("REC::Track"));
        hipo::bank TRAJ(factory.getSchema("REC::Traj"));
        hipo::bank CAL(factory.getSchema("REC::Calorimeter"));
        hipo::bank CHKV(factory.getSchema("REC::Cherenkov"));
        hipo::bank HEL_SCALER(factory.getSchema("HEL::scaler"));
        hipo::bank RUN_SCALER(factory.getSchema("RUN::scaler"));
        hipo::bank CONF(factory.getSchema("RUN::config"));
        hipo::bank RASTER(factory.getSchema("RASTER::position"));

      while(reader.next()==true){ // #1
	reader.read(event);
	event.getStructure(PART);
	event.getStructure(HEL);
	event.getStructure(TRACK);
	event.getStructure(TRAJ);
	event.getStructure(CAL);
	event.getStructure(CHKV);
	event.getStructure(CONF);
        event.getStructure(HEL_SCALER);
        event.getStructure(RUN_SCALER);
	event.getStructure(RASTER);
	int evnum = CONF.getInt("event",0);
	int runnum= CONF.getInt("run",0); // Already read in, but this is a double-check

      if( qa->Pass( runnum, evnum ) ) {
	
	int helRows = HEL_SCALER.getRows(); 
	if( helRows > 0 ){
          //for(int row=0; row<HEL_SCALER.getRows(); row++){
	      int hel = corr_factor * HEL_SCALER.getByte("helicity",0);
	      // I swapped the helicity signs to match what I'm using for e- spins
	   for(int h=0; h<helRows; h++){
	      if(hel == 1){
		     FCup_hel_n += HEL_SCALER.getFloat("fcupgated",h);
		     FCup_ughel_n += HEL_SCALER.getFloat("fcup",h);
	      }
	      else if(hel == -1){
		     FCup_hel_p += HEL_SCALER.getFloat("fcupgated",h);
		     FCup_ughel_p += HEL_SCALER.getFloat("fcup",h);
	      }
	      else if(hel == 0){
		     FCup_hel_0 += HEL_SCALER.getFloat("fcupgated",h);
		     FCup_ughel_0 += HEL_SCALER.getFloat("fcup",h);
	      }
	  }
        }

	qa->AccumulateCharge();   // Accumulate gated FC charge from QADB (from RUN::scaler I think??)
	qa->AccumulateChargeHL(); // Accumulate gated, helicity-latched FC charge (from HEL::scaler)

        vector<Particle> Part = SetParticleData( PID, qa, corr_factor, CONF, PART, HEL, TRACK, TRAJ, CAL, CHKV, RASTER );
	//   if( theta >= 5.0 && theta < 40.0 && w > 2.0 && Q2 > 1.0 && E > 2.6 && passedDCcut  ){

	//if( Part.size() == 1 ){
	for(auto P : Part ){

	    pid = P.PID; 
	    Status = P.Status;
	    Chi2PID = P.Chi2PID;
	    Helicity = P.Helicity;
	    Charge = P.Charge;
	    Px = P.Px; Py = P.Py; Pz = P.Pz;
	    Vx = P.Vx; Vy = P.Vy; Vz = P.Vz;
	    TrackChi2NDF = P.TrackChi2NDF;
	    EdgeR1 = P.EdgeR1; EdgeR2 = P.EdgeR2; EdgeR3 = P.EdgeR3;
	    SectorDC = P.SectorDC; SectorECAL = P.SectorECAL; SectorHTCC = P.SectorHTCC;
	    PcalX = P.PcalX; PcalY = P.PcalY;
	    PcalE = P.PcalE;
	    PcalLu = P.PcalLu; PcalLv = P.PcalLv; PcalLw = P.PcalLw;
	    EcinX = P.EcinX; EcinY = P.EcinY;
	    EcinE = P.EcinE;
	    EcinLu = P.EcinLu; EcinLv = P.EcinLv; EcinLw = P.EcinLw;
	    EoutX = P.EoutX; EoutY = P.EoutY;
	    EoutE = P.EoutE;
	    EoutLu = P.EoutLu; EoutLv = P.EoutLv; EoutLw = P.EoutLw;

	    nphe = P.nphe;
	    Beta = P.Beta;

	    RasterX = P.RasterX;
	    RasterY = P.RasterY;

	    M2_PCAL = P.M2_PCAL;
	    M2_ECIN = P.M2_ECIN;
	    M2_EOUT = P.M2_EOUT;

	    // Only select particles in the forward detector
	    if( abs(Status) < 4000 ){
		tree->Fill();
   	        event_counter++;
	    }
	}
	//event_counter++;
      } // End of QA pass "if" statement
    } // End of "reader.next()" while loop #1

   } // end of loop over all files on TChain

   // Write the raster constants
   string rasterData = "/volatile/clas12/holmberg/";
   //string rasterData = "/w/hallb-scshelf2102/clas12/holmberg/Dilution_Factors/";
   if( runToAnalyze < 16800 ) rasterData += "Summer_22_Data/Text_Files/"+ outPID+"_"+TARGET_TYPE+"_"+to_string(runToAnalyze)+"_"+fileExt+"_FC_Info.txt";
   else if( runToAnalyze > 16800 && runToAnalyze < 17410 ) rasterData += "Fall_22_Data/Text_Files/"+ outPID+"_"+TARGET_TYPE+"_"+to_string(runToAnalyze)+"_"+fileExt+"_FC_Info.txt";
   else if( runToAnalyze > 17450 ) rasterData += "Spring_23_Data/Text_Files/"+ outPID+"_"+TARGET_TYPE+"_"+to_string(runToAnalyze)+"_"+fileExt+"_FC_Info.txt";

   ofstream rasterOut( rasterData.c_str() );
   rasterOut << "UG_HelN   UG_HelP   UG_Hel0   G_HelN   G_HelP   G_Hel0   QA_HelN   QA_HelP   QA_Hel0   RUN_Hel\n";
   rasterOut << FCup_ughel_n <<"	"<< FCup_ughel_p <<"	"<< FCup_ughel_0 <<"	"<< FCup_hel_n <<"	"<< FCup_hel_p <<"	"<< FCup_hel_0 <<"	"
	     << qa->GetAccumulatedChargeHL(1) <<"	"<< qa->GetAccumulatedChargeHL(-1) <<"	"<< qa->GetAccumulatedChargeHL(0) <<"	"<< qa->GetAccumulatedCharge() << endl;
   rasterOut.close();

   // Update charges
   totalRUNcharge += qa->GetAccumulatedCharge()     - totalRUNcharge;
   totalHLMcharge += qa->GetAccumulatedChargeHL(1)  - totalHLMcharge;
   totalHLPcharge += qa->GetAccumulatedChargeHL(-1) - totalHLPcharge;
   totalHL0charge += qa->GetAccumulatedChargeHL(0)  - totalHL0charge;

   tree->Write("Particle",TObject::kOverwrite); // Removes the autosave branches for simplicity
   
   file->Close();
   cout << "Total events = "<< event_counter << endl;
   event_counter = 0;


  } // End of if statement checking if this is a good run to look at
 //} // End of for loop on target types
 //} // End of for loop on input files

}
