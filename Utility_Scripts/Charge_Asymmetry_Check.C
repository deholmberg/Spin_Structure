#include "../Header_Files/Binning_Classes.h"

using namespace std;

// This struct holds the information for each individual run
struct RunFCData{

	int Run = 0; // Run Number
	string Target = "N/A"; // Target type
	double Np    = 0; // helicity +1 counts
	double Nm    = 0; // helicity -1 counts
	double Fp_ug = 0; // FC +1 charge (ungated)
	double Fm_ug = 0; // FC -1 charge (ungated)
	double Fp_g  = 0; // FC +1 charge (gated)
	double Fm_g  = 0; // FC -1 charge (gated)
	double ReadP = 0; // Number of helicity +1 readouts of HEL::scaler from the updated list
	double ReadN = 0; // Number of helicity -1 readouts of HEL::scaler from the updated list

	// This holds the updated FC asymmetries
	double BCA_EPICS_FCUP     = 0; // BCA from the Fcup EPICS
	double BCA_EPICS_SLM      = 0; // BCA from the Slm EPICS
	double BCA_STRUCK_gated   = 0; // BCA from the gated STRUCK readouts
	double BCA_STRUCK_ungated = 0; // BCA from the ungated STRUCK readouts

	void SetValues( int run, string targ, double np, double nm, double fp_ug, double fm_ug, double fp_g, double fm_g ){
		Run = run; Target = targ;
		Np = np; Nm = nm; Fp_ug = fp_ug; Fm_ug = fm_ug;
		Fp_g = fp_g; Fm_g = fm_g;
	}
	
	void SetBCA( double epicsFC, double epicsSLM, double gStruck, double ugStruck, double countsP, double countsN ){
		BCA_EPICS_FCUP = epicsFC; BCA_EPICS_SLM = epicsSLM;
		BCA_STRUCK_gated = gStruck; BCA_STRUCK_ungated = ugStruck;
		ReadP = countsP; ReadN = countsN;
	}

	// This opens a run file for a given run and target type and sets the relevant info
	void SetCounts( int run, string targ, bool useQA=true ){
	    string filePath;
	    if( useQA ) filePath = "../Latest_Skims/Text_Files/"+ targ +"_"+ to_string(run) +"_DF_Data.txt";
	    else filePath = "../Latest_Skims/Text_Files_NoQA/"+ targ +"_"+ to_string(run) +"_DF_Data_NoQA.txt";

	    ifstream fin( filePath.c_str() );
	    if( fin.fail() ){ cout << "ERROR: Couldn't open file "<< filePath <<" to set updated BCA info.\n"; return; }

	    string line;
	    getline( fin, line ); // Throw away header row
	    double NP_total = 0;
	    double NM_total = 0;
	    double FP_total = 0;
	    double FM_total = 0;
	    while( getline( fin, line ) ){
		stringstream sin(line);
		double Q2_Min, Q2_Max, X_Min, X_Max, Np_Counts, Nm_Counts,  FCp_Charge, FCm_Charge, N0_Counts, FC0_Charge, X_Avg, Q2_Avg;
		sin >> Q2_Min >> Q2_Max >> X_Min >> X_Max >> Np_Counts >> Nm_Counts >>  FCp_Charge >> FCm_Charge >> N0_Counts >> FC0_Charge >> X_Avg >> Q2_Avg;
		// Accumulate the counts
		NP_total += Np_Counts; NM_total += Nm_Counts;
		// The FC charge is the same for each bin, so just override which is dumb but it works
		FP_total = FCp_Charge; FM_total = FCm_Charge;
	    }
	    fin.close();

	    Np = NP_total; Nm = NM_total;
	    Fp_g = FP_total; Fm_g = FM_total; // Only the gated values are set
	    Run = run;
	    Target = targ;
	
	//Q2_Min   Q2_Max   X_Min   X_Max   Np_Counts   Nm_Counts    FCp_Charge   FCm_Charge   N0_Counts   FC0_Charge   X_Avg   Q2_Avg
	//0   0   0.125   0.15   120   114   2759.29   2747.69   0   0   0.148904   2.21101

	}

	vector<double> FCAsymmetry( string gate, bool normToReadouts=false ){

		vector<double> FC_Asyms {0,0};

		double Fp = 0; double Fm = 0; // The FC charges
		double d_Fp = ReadP / 906.2; // Uncertainties on each FC charge
		double d_Fm = ReadN / 906.2;

		if( gate == "Ungated" ){ Fp = Fp_ug; Fm = Fm_ug; }
		else if( gate == "Gated" ){ Fp = Fp_g; Fm = Fm_g; }
		else if( gate == "EPICS_FCUP"      ){ FC_Asyms[0] = BCA_EPICS_FCUP; return FC_Asyms; }
		else if( gate == "EPICS_SLM"       ){ FC_Asyms[0] = BCA_EPICS_SLM; return FC_Asyms; }
		else if( gate == "STRUCK_gated"    ){ FC_Asyms[0] = BCA_STRUCK_gated; return FC_Asyms; }
		else if( gate == "STRUCK_ungated" ){ FC_Asyms[0] = BCA_STRUCK_ungated; return FC_Asyms; }

		if( Fp + Fm > 0 && !normToReadouts ){
			FC_Asyms[0] = (Fp - Fm) / (Fp + Fm);
			FC_Asyms[1] = (2.0/pow((Fp+Fm),2)) * sqrt( Fm*Fm*d_Fp + Fp*Fp*d_Fm );
		}
		else if( Fp +Fm > 0 && ReadP > 0 && ReadN > 0 && normToReadouts ){
			FC_Asyms[0] = (1.0*Fp/ReadP - 1.0*Fm/ReadN) / (1.0*Fp/ReadP + 1.0*Fm/ReadN);
			// Check: does dividing by ReadP and ReadN introduce an extra error that must be propagated?
			FC_Asyms[1] = (2.0/pow( ( (Fp/ReadP) +(Fm/ReadN) ),2)) * sqrt( Fm*Fm*d_Fp/(ReadN*ReadN*ReadP) + Fp*Fp*d_Fm/(ReadP*ReadP*ReadN) );
		}

		return FC_Asyms;
	}

	vector<double> RawAsymmetry() const{

		vector<double> Asyms = {0,0};

		if( Np + Nm > 0 ){
			Asyms[0] = 1.0*(Np - Nm)/(Np + Nm);
			Asyms[1] = (2.0/pow(1.0*(Np+Nm),2)) * sqrt( 1.0*Nm*Nm*Np + 1.0*Np*Np*Nm );
		}
		
		return Asyms;
	}

	vector<double> NormAsymmetry(string normType="Gated") const{

		// The first is the value of the asym, the second the stat error of counts and FC charge, and the third
		// is the stat error of the counts only
		vector<double> Asyms = {0,0,0};

		double Fp = 0; double Fm = 0;
		double d_Fp = ReadP / 906.2; // Uncertainties on each FC charge
		double d_Fm = ReadN / 906.2;

		if( normType == "Gated" ){ Fp = Fp_g; Fm = Fm_g; }
		else if( normType == "Ungated" ){ Fp = Fp_ug; Fm = Fm_ug; }

		if( (Np + Nm) > 0 && Fp > 0 && Fm > 0 ){
			double np = 1.0 * Np / Fp;
			double nm = 1.0 * Nm / Fm;

			Asyms[0] = 1.0*(np - nm)/(np + nm);
			//Asyms[1] = (2.0/pow(1.0*(np+nm),2)) * sqrt( 1.0*nm*nm*np/Fp + 1.0*np*np*nm/Fm );
			//double asymRawStat = (2.0/pow(1.0*(Np+Nm),2)) * sqrt( 1.0*Nm*Nm*Np + 1.0*Np*Np*Nm );
			double asymRawStat = (2.0/pow(1.0*(np+nm),2)) * sqrt( 1.0*nm*nm*np/Fp + 1.0*np*np*nm/Fm );
			double asymFCStat  = (2.0/pow((Fp+Fm),2)) * sqrt( Fm*Fm*d_Fp + Fp*Fp*d_Fm );
			Asyms[1] = sqrt( pow(asymRawStat,2) + pow(asymFCStat,2) );
			Asyms[2] = asymRawStat;
		}

		return Asyms;
	}

	void Print() const{
		cout <<"======== Run: "<< Run <<" ========\n";
		cout <<"---> Target: "<< Target << endl;
		cout <<"---> Nm = "<< Nm << endl;
		cout <<"---> Np = "<< Np << endl;
		cout <<"---> Gated Fm   = "<< Fm_g << endl;
		cout <<"---> Gated Fp   = "<< Fp_g << endl;
		cout <<"---> Scaler Readouts -1 = "<< ReadN << endl;
		cout <<"---> Scaler Readouts +1 = "<< ReadP << endl;
		cout <<"---> BCA_EPICS_FCUP     = "<< BCA_EPICS_FCUP << endl;
		cout <<"---> BCA_EPICS_SLM      = "<< BCA_EPICS_SLM << endl;
		cout <<"---> BCA_STRUCK_gated   = "<< BCA_STRUCK_gated << endl;
		cout <<"---> BCA_STRUCK_ungated = "<< BCA_STRUCK_ungated << endl;
		cout <<"============================\n";
	}
};

// This struct holds all the FC info for a set of runs
struct AllFCRunInfo{

	// Holds all the run info
	vector<RunFCData> AllRuns;

	// Using an input vector of runs, sets the FC info
	void SetRunInfo( RunPeriod& RP, vector<int>& Runs, bool useQA=true ){
	
	    for( int run : Runs ){
	        string Target = RP.getTargetType( run ); 
		if( Target == "Foil" ) Target = "F";
		else if( Target == "Empty" ) Target = "ET";
		RunFCData thisRun;
		thisRun.SetCounts( run, Target, useQA );
		AllRuns.push_back( thisRun );
	    }

	    int totalRuns = AllRuns.size();

	    // Now set the updated FC asym values
	    ifstream fin("/w/hallb-scshelf2102/clas12/holmberg/bcadb/db/rgc_su22.txt");
	    string line;
	    getline( fin, line ); // Throw away header row
	    while( getline( fin, line ) ){
//#runnum       BCA_EPICS_FCUP        BCA_EPICS_SLM     BCA_STRUCK_gated   BCA_STRUCK_ungated  count_STRUCK_helPos  count_STRUCK_helNeg
//16042                0.00185              0.00177              0.00000              0.00000                    0                    0
 		stringstream sin(line);
		double epicsfc, epicsslm, gated, ungated, countsPos, countsNeg;
		int run;
		sin >> run >> epicsfc >> epicsslm >> gated >> ungated >> countsPos >> countsNeg;
		// This is inefficient, but loop through the "AllRuns" vector and set the info
		if( run >= 16128 ){
		    for( int i=0; i<totalRuns; i++){
			if( AllRuns[i].Run == run ) AllRuns[i].SetBCA( epicsfc, epicsslm, gated, ungated, countsPos, countsNeg );
		    }
		}
	    }
	    fin.close(); 

	}

	// Prints all the info
	void Print() const{
	    for(auto r : AllRuns) r.Print();
	}

	// Returns a TGraphErrors object with the desired double-spin asymmetry plotted against the desired type of FC asymmetry
	TGraphErrors* MakePlot( string Target, string Period, string BCA_Type, int Color=1, bool UseNorm=false ){
	    int RunMin = 0; int RunMax = 0; // Bounds of potential run ranges to use for plotting
	    if( Period == "Summer" ){ RunMin = 16128; RunMax = 16772; }

	    vector<double> BCA, ALL, ALL_Err;
	    // Loop over and get the desired info
	    for( auto Run : AllRuns ){
	        vector<double> FC_Asyms = Run.FCAsymmetry( BCA_Type );
		vector<double> ALL_Asym = Run.RawAsymmetry();

		// Check the target with a target pass
		bool TargetCheck = false;
		string targ = Run.Target;

		//if( Target == "Unpol" && ( targ == "F" || targ == "ET" || targ == "C" || targ == "CH2" || targ == "CD2" ) ) TargetCheck = true;
		if( Target == "Unpol" && ( targ == "ET" || targ == "C" || targ == "CH2" || targ == "CD2" ) ) TargetCheck = true; // removed foil targets for now...
		else TargetCheck = Target == targ; 

		if( TargetCheck && Run.Run >= RunMin && Run.Run <= RunMax && FC_Asyms[0] != 0 ){
		    // This vector holds the FC info. Ignore the second entry (the error) for now...
		    BCA.push_back( FC_Asyms[0] );

		    double all = ALL_Asym[0];
		    if( UseNorm ) all = all - FC_Asyms[0]; 

		    ALL.push_back( all );
		    // FIXME: The error isn't scaled when subtracting the FC asym
		    ALL_Err.push_back( ALL_Asym[1] );
		    cout <<"Run "<< Run.Run << endl;
		}
	    }

	    TGraphErrors* gr = new TGraphErrors( BCA.size(), BCA.data(), ALL.data(), nullptr, ALL_Err.data() );
	    string norm = "LU";
	    if( UseNorm ) norm = "LU,norm";
	    gr->SetTitle( string( Period + " A_{"+norm+"} vs. A_{FC} for "+ Target +" Targets; A_{FC}; A_{"+norm+"}" ).c_str() );
	    gr->SetMarkerStyle(kFullCircle);
	    gr->SetMarkerColor( Color );

	    return gr;
	}

	// Creates a plot of all the FCup types drawn together
	TCanvas* MakeAllFCPlots( string Target, string Period, string Version="All", string TitleExt="", double yMin=-0.01, double yMax=0.01 ){

	    // Used to track what the axis ranges should be
	    //bool useUnpolTargets = !( Target == "NH3" or Target == "ND3" );	    

	    // Make multigraph and legend to hold everything
	    TMultiGraph* mg = new TMultiGraph();
	    TLegend* leg = new TLegend(0.18,0.1, 0.98, 0.3); leg->SetNColumns(2);
	    string canName = Target + Period +"_can";
	    TCanvas* c = new TCanvas( canName.c_str(), canName.c_str(), 1600, 600 ); c->Divide(2,1);
	    
	    // Plot the unnormalized graphs
	    TFitResultPtr r;
	    c->cd(1);
	    c->GetPad(1)->SetLeftMargin(0.18);
	    c->GetPad(1)->SetRightMargin(0.02);

	    auto gr_Gatedz = this->MakePlot( Target, Period, "Gated", 1, false );
	    r = gr_Gatedz->Fit("pol1","S","Q");
	    gr_Gatedz->GetFunction("pol1")->SetLineColor(1);
	    mg->Add( gr_Gatedz, "p" );
	    stringstream leg_Gatedz; leg_Gatedz << "Gated A_{FC}, Slope = "<< setprecision(3) << r->Value(1) <<" #pm "<< r->Error(1);
	    leg->AddEntry( gr_Gatedz, leg_Gatedz.str().c_str(), "p" );

	    if( Version == "All" or Version == "EPICS_FCUP" ){
		auto gr_EPICS_FCUPz = this->MakePlot( Target, Period, "EPICS_FCUP", 2, false );
	        r = gr_EPICS_FCUPz->Fit("pol1","S","Q");
	        gr_EPICS_FCUPz->GetFunction("pol1")->SetLineColor(2);
	        mg->Add( gr_EPICS_FCUPz, "p" );
	        stringstream leg_EPICS_FCUPz; leg_EPICS_FCUPz << "EPICS_FCUP A_{FC}, Slope = "<< setprecision(3) << r->Value(1) <<" #pm "<< r->Error(1);
	        leg->AddEntry( gr_EPICS_FCUPz, leg_EPICS_FCUPz.str().c_str(), "p" );
	    }
	    if( Version == "All" or Version == "EPICS_SLM" ){
	        auto gr_EPICS_SLMz = this->MakePlot( Target, Period, "EPICS_SLM", 3, false );
	        r = gr_EPICS_SLMz->Fit("pol1","S","Q");
	        gr_EPICS_SLMz->GetFunction("pol1")->SetLineColor(3);
	        mg->Add( gr_EPICS_SLMz, "p" );
	        stringstream leg_EPICS_SLMz; leg_EPICS_SLMz << "EPICS_SLM A_{FC}, Slope = "<< setprecision(3) << r->Value(1) <<" #pm "<< r->Error(1);
	        leg->AddEntry( gr_EPICS_SLMz, leg_EPICS_SLMz.str().c_str(), "p" );
	    }
	    if( Version == "All" or Version == "STRUCK_gated" ){
	        auto gr_STRUCK_gatedz = this->MakePlot( Target, Period, "STRUCK_gated", 4, false );
	        r = gr_STRUCK_gatedz->Fit("pol1","S","Q");
	        gr_STRUCK_gatedz->GetFunction("pol1")->SetLineColor(4);
	        mg->Add( gr_STRUCK_gatedz, "p" );
	        stringstream leg_STRUCK_gatedz; leg_STRUCK_gatedz << "STRUCK_gated A_{FC}, Slope = "<< setprecision(3) << r->Value(1) <<" #pm "<< r->Error(1);
	        leg->AddEntry( gr_STRUCK_gatedz, leg_STRUCK_gatedz.str().c_str(), "p" );
	    }
	    if( Version == "All" or Version == "STRUCK_ungated" ){
	        auto gr_STRUCK_ungatedz = this->MakePlot( Target, Period, "STRUCK_ungated", 6, false );
	        r = gr_STRUCK_ungatedz->Fit("pol1","S","Q");
	        gr_STRUCK_ungatedz->GetFunction("pol1")->SetLineColor(6);
	        mg->Add( gr_STRUCK_ungatedz, "p" );
	        stringstream leg_STRUCK_ungatedz; leg_STRUCK_ungatedz << "STRUCK_ungated A_{FC}, Slope = "<< setprecision(3) << r->Value(1) <<" #pm "<< r->Error(1);
	        leg->AddEntry( gr_STRUCK_ungatedz, leg_STRUCK_ungatedz.str().c_str(), "p" );
	    }
	    mg->SetTitle( string( Period +" A_{LU} vs. A_{FC} for "+ Target +" "+TitleExt+"; A_{FC}; A_{LU}").c_str() ); 
	    mg->GetYaxis()->SetRangeUser( yMin, yMax );
	    mg->Draw("ap"); leg->Draw("same");

	    c->cd(2);
	    c->GetPad(2)->SetLeftMargin(0.18);
	    c->GetPad(2)->SetRightMargin(0.02);

	    TMultiGraph* mgNorm = new TMultiGraph();
	    TLegend* legNorm = new TLegend(0.18,0.7, 0.98, 0.9); legNorm->SetNColumns(2);

	    auto gr_GatedNorm = this->MakePlot( Target, Period, "Gated", 1, true );
	    r = gr_GatedNorm->Fit("pol1","S","Q");
	    gr_GatedNorm->GetFunction("pol1")->SetLineColor(1);
	    mgNorm->Add( gr_GatedNorm, "p" );
	    stringstream legNorm_GatedNorm; legNorm_GatedNorm << "Norm Gated A_{FC}, Slope = "<< setprecision(3) << r->Value(1) <<" #pm "<< r->Error(1);
	    legNorm->AddEntry( gr_GatedNorm, legNorm_GatedNorm.str().c_str(), "p" );

	    if( Version == "All" or Version == "EPICS_FCUP" ){
	        auto gr_EPICS_FCUPNorm = this->MakePlot( Target, Period, "EPICS_FCUP", 2, true );
	        r = gr_EPICS_FCUPNorm->Fit("pol1","S","Q");
	        gr_EPICS_FCUPNorm->GetFunction("pol1")->SetLineColor(2);
	        mgNorm->Add( gr_EPICS_FCUPNorm, "p" );
	        stringstream legNorm_EPICS_FCUPNorm; legNorm_EPICS_FCUPNorm << "EPICS_FCUP A_{FC}, Slope = "<< setprecision(3) << r->Value(1) <<" #pm "<< r->Error(1);
	        legNorm->AddEntry( gr_EPICS_FCUPNorm, legNorm_EPICS_FCUPNorm.str().c_str(), "p" );
	    }
	    if( Version == "All" or Version == "EPICS_SLM" ){
	        auto gr_EPICS_SLMNorm = this->MakePlot( Target, Period, "EPICS_SLM", 3, true );
	        r = gr_EPICS_SLMNorm->Fit("pol1","S","Q");
	        gr_EPICS_SLMNorm->GetFunction("pol1")->SetLineColor(3);
	        mgNorm->Add( gr_EPICS_SLMNorm, "p" );
	        stringstream legNorm_EPICS_SLMNorm; legNorm_EPICS_SLMNorm << "EPICS_SLM A_{FC}, Slope = "<< setprecision(3) << r->Value(1) <<" #pm "<< r->Error(1);
	        legNorm->AddEntry( gr_EPICS_SLMNorm, legNorm_EPICS_SLMNorm.str().c_str(), "p" );
	    }
	    if( Version == "All" or Version == "STRUCK_gated" ){
	        auto gr_STRUCK_gatedNorm = this->MakePlot( Target, Period, "STRUCK_gated", 4, true );
	        r = gr_STRUCK_gatedNorm->Fit("pol1","S","Q");
	        gr_STRUCK_gatedNorm->GetFunction("pol1")->SetLineColor(4);
	        mgNorm->Add( gr_STRUCK_gatedNorm, "p" );
	        stringstream legNorm_STRUCK_gatedNorm; legNorm_STRUCK_gatedNorm << "STRUCK_gated A_{FC}, Slope = "<< setprecision(3) << r->Value(1) <<" #pm "<< r->Error(1);
	        legNorm->AddEntry( gr_STRUCK_gatedNorm, legNorm_STRUCK_gatedNorm.str().c_str(), "p" );
	    }
	    if( Version == "All" or Version == "STRUCK_ungated" ){
	        auto gr_STRUCK_ungatedNorm = this->MakePlot( Target, Period, "STRUCK_ungated", 6, true );
	        r = gr_STRUCK_ungatedNorm->Fit("pol1","S","Q");
	        gr_STRUCK_ungatedNorm->GetFunction("pol1")->SetLineColor(6);
	        mgNorm->Add( gr_STRUCK_ungatedNorm, "p" );
	        stringstream legNorm_STRUCK_ungatedNorm; legNorm_STRUCK_ungatedNorm << "STRUCK_ungated A_{FC}, Slope = "<< setprecision(3) << r->Value(1) <<" #pm "<< r->Error(1);
	        legNorm->AddEntry( gr_STRUCK_ungatedNorm, legNorm_STRUCK_ungatedNorm.str().c_str(), "p" );
	    }

	    mgNorm->SetTitle( string( Period +" A_{norm} vs. A_{FC} for "+ Target +" "+TitleExt+"; A_{FC}; A_{norm}").c_str() ); 
	    mgNorm->GetYaxis()->SetRangeUser( yMin, yMax );
	    mgNorm->Draw("ap"); legNorm->Draw("same");

	    return c;
	}
};

void Charge_Asymmetry_Check( ){

//getDataSet(string Dataset, string TorPol, string SolPol)
	
	RunPeriod RP;
	vector<int> AllSummerRuns = RP.getDataSet("Su22","Neg","Neg");
        //vector<int> getRunsForHWPandTPol( string DataSet, string Target, int HWPstatus, int TPolSign ){
        vector<int> AllHWPInRuns = RP.getRunsForHWPandTPol( "Su22", "Unpol", 1, 0 );
	//for(int r : AllHWPInRuns ) cout << r << endl;

	AllFCRunInfo AllData, AllDataNoQA; 
	AllData.SetRunInfo( RP, AllSummerRuns );
	AllDataNoQA.SetRunInfo( RP, AllSummerRuns, false );

	// Make a plot of the summer carbon data with QA cuts applied
	auto Plt_C_v1 = AllData.MakeAllFCPlots( "C", "Summer", "EPICS_FCUP" );     Plt_C_v1->Print("PDF_Plots/Summer_Carbon_FCup.pdf(");
	auto Plt_C_v2 = AllData.MakeAllFCPlots( "C", "Summer", "EPICS_SLM" );      Plt_C_v2->Print("PDF_Plots/Summer_Carbon_FCup.pdf");
	auto Plt_C_v3 = AllData.MakeAllFCPlots( "C", "Summer", "STRUCK_gated" );   Plt_C_v3->Print("PDF_Plots/Summer_Carbon_FCup.pdf");
	auto Plt_C_v4 = AllData.MakeAllFCPlots( "C", "Summer", "STRUCK_ungated" ); Plt_C_v4->Print("PDF_Plots/Summer_Carbon_FCup.pdf)");

	// Do the same without QA cuts
	auto NoQA_C_v1 = AllDataNoQA.MakeAllFCPlots( "C", "Summer", "EPICS_FCUP", "(No QA)" );     NoQA_C_v1->Print("PDF_Plots/Summer_Carbon_FCup_NoQA.pdf(");
	auto NoQA_C_v2 = AllDataNoQA.MakeAllFCPlots( "C", "Summer", "EPICS_SLM", "(No QA)" );      NoQA_C_v2->Print("PDF_Plots/Summer_Carbon_FCup_NoQA.pdf");
	auto NoQA_C_v3 = AllDataNoQA.MakeAllFCPlots( "C", "Summer", "STRUCK_gated", "(No QA)" );   NoQA_C_v3->Print("PDF_Plots/Summer_Carbon_FCup_NoQA.pdf");
	auto NoQA_C_v4 = AllDataNoQA.MakeAllFCPlots( "C", "Summer", "STRUCK_ungated", "(No QA)" ); NoQA_C_v4->Print("PDF_Plots/Summer_Carbon_FCup_NoQA.pdf)");

	// Make a plot of the summer carbon data with QA cuts applied
	auto Plt_Unpol_v1 = AllData.MakeAllFCPlots( "Unpol", "Summer", "EPICS_FCUP" );     Plt_Unpol_v1->Print("PDF_Plots/Summer_Unpol_FCup.pdf(");
	auto Plt_Unpol_v2 = AllData.MakeAllFCPlots( "Unpol", "Summer", "EPICS_SLM" );      Plt_Unpol_v2->Print("PDF_Plots/Summer_Unpol_FCup.pdf");
	auto Plt_Unpol_v3 = AllData.MakeAllFCPlots( "Unpol", "Summer", "STRUCK_gated" );   Plt_Unpol_v3->Print("PDF_Plots/Summer_Unpol_FCup.pdf");
	auto Plt_Unpol_v4 = AllData.MakeAllFCPlots( "Unpol", "Summer", "STRUCK_ungated" ); Plt_Unpol_v4->Print("PDF_Plots/Summer_Unpol_FCup.pdf)");

	// Do the same without QA cuts
	auto NoQA_Unpol_v1 = AllDataNoQA.MakeAllFCPlots( "Unpol", "Summer", "EPICS_FCUP", "(No QA)" );     NoQA_Unpol_v1->Print("PDF_Plots/Summer_Unpol_FCup_NoQA.pdf(");
	auto NoQA_Unpol_v2 = AllDataNoQA.MakeAllFCPlots( "Unpol", "Summer", "EPICS_SLM", "(No QA)" );      NoQA_Unpol_v2->Print("PDF_Plots/Summer_Unpol_FCup_NoQA.pdf");
	auto NoQA_Unpol_v3 = AllDataNoQA.MakeAllFCPlots( "Unpol", "Summer", "STRUCK_gated", "(No QA)" );   NoQA_Unpol_v3->Print("PDF_Plots/Summer_Unpol_FCup_NoQA.pdf");
	auto NoQA_Unpol_v4 = AllDataNoQA.MakeAllFCPlots( "Unpol", "Summer", "STRUCK_ungated", "(No QA)" ); NoQA_Unpol_v4->Print("PDF_Plots/Summer_Unpol_FCup_NoQA.pdf)");


	// Loop over all HWP and Tpol configs for NH3
	vector<int> Target_Pols = {-1,1};
	vector<int> HWP_Status  = { 0,1};
/*
	// NH3 Runs
	for(size_t i=0; i<HWP_Status.size(); i++){
	    for(size_t j=0; j<Target_Pols.size(); j++){
	        vector<int> Runs = RP.getRunsForHWPandTPol( "Su22", "NH3", HWP_Status[i], Target_Pols[j] );
	        AllFCRunInfo thisData; thisData.SetRunInfo( RP, Runs );

		string TitleExt = "HWP = "+to_string( HWP_Status[i] )+", P_{t} = "+to_string( Target_Pols[j] );
		double yMin = -0.01; double yMax = 0.01; // Bounds for the plots
		if( Target_Pols[j] < 0 ){ yMin = -0.04; yMax = -0.01; }
		else if( Target_Pols[j] > 0 ){ yMin = 0.01; yMax = 0.04; }

	        auto Plt_NH3_v1 = thisData.MakeAllFCPlots( "NH3", "Summer", "EPICS_FCUP", TitleExt, yMin, yMax );
		Plt_NH3_v1->Print( string("PDF_Plots/Summer_NH3_FCup_HWP"+to_string( HWP_Status[i] )+"_Tpol"+to_string( Target_Pols[j] )+".pdf(").c_str() );

	        auto Plt_NH3_v2 = thisData.MakeAllFCPlots( "NH3", "Summer", "EPICS_SLM", TitleExt, yMin, yMax );
		Plt_NH3_v2->Print( string("PDF_Plots/Summer_NH3_FCup_HWP"+to_string( HWP_Status[i] )+"_Tpol"+to_string( Target_Pols[j] )+".pdf").c_str() );

	        auto Plt_NH3_v3 = thisData.MakeAllFCPlots( "NH3", "Summer", "STRUCK_gated", TitleExt, yMin, yMax );
		Plt_NH3_v3->Print( string("PDF_Plots/Summer_NH3_FCup_HWP"+to_string( HWP_Status[i] )+"_Tpol"+to_string( Target_Pols[j] )+".pdf").c_str() );

	        auto Plt_NH3_v4 = thisData.MakeAllFCPlots( "NH3", "Summer", "STRUCK_ungated", TitleExt, yMin, yMax );
		Plt_NH3_v4->Print( string("PDF_Plots/Summer_NH3_FCup_HWP"+to_string( HWP_Status[i] )+"_Tpol"+to_string( Target_Pols[j] )+".pdf)").c_str() );
	    }
	}

	// ND3 Runs
	for(size_t i=0; i<HWP_Status.size(); i++){
	    for(size_t j=0; j<Target_Pols.size(); j++){
	        vector<int> Runs = RP.getRunsForHWPandTPol( "Su22", "ND3", HWP_Status[i], Target_Pols[j] );
	        AllFCRunInfo thisData; thisData.SetRunInfo( RP, Runs );

		string TitleExt = "HWP = "+to_string( HWP_Status[i] )+", P_{t} = "+to_string( Target_Pols[j] );
		double yMin = -0.01; double yMax = 0.01; // Bounds for the plots
		if( Target_Pols[j] < 0 ){ yMin = -0.03; yMax = 0.03; }
		else if( Target_Pols[j] > 0 ){ yMin = -0.03; yMax = 0.03; }

	        auto Plt_ND3_v1 = thisData.MakeAllFCPlots( "ND3", "Summer", "EPICS_FCUP", TitleExt, yMin, yMax );
		Plt_ND3_v1->Print( string("PDF_Plots/Summer_ND3_FCup_HWP"+to_string( HWP_Status[i] )+"_Tpol"+to_string( Target_Pols[j] )+".pdf(").c_str() );

	        auto Plt_ND3_v2 = thisData.MakeAllFCPlots( "ND3", "Summer", "EPICS_SLM", TitleExt, yMin, yMax );
		Plt_ND3_v2->Print( string("PDF_Plots/Summer_ND3_FCup_HWP"+to_string( HWP_Status[i] )+"_Tpol"+to_string( Target_Pols[j] )+".pdf").c_str() );

	        auto Plt_ND3_v3 = thisData.MakeAllFCPlots( "ND3", "Summer", "STRUCK_gated", TitleExt, yMin, yMax );
		Plt_ND3_v3->Print( string("PDF_Plots/Summer_ND3_FCup_HWP"+to_string( HWP_Status[i] )+"_Tpol"+to_string( Target_Pols[j] )+".pdf").c_str() );

	        auto Plt_ND3_v4 = thisData.MakeAllFCPlots( "ND3", "Summer", "STRUCK_ungated", TitleExt, yMin, yMax );
		Plt_ND3_v4->Print( string("PDF_Plots/Summer_ND3_FCup_HWP"+to_string( HWP_Status[i] )+"_Tpol"+to_string( Target_Pols[j] )+".pdf)").c_str() );
	    }
	}
*/


	// Also generate plots for the NH3 runs of each targer polarization and HWP configuration...
        //vector<int> NH3_In_Pos_Runs = RP.getRunsForHWPandTPol( "Su22", "NH3", 1, 1 );
	//AllFCRunInfo Data_NH3_In_Pos; Data_NH3_In_Pos.SetRunInfo( RP, NH3_In_Pos_Runs );
	//auto Plot_NH3_In_Pos = Data_NH3_In_Pos.MakeAllFCPlots( "NH3", "Summer" );

	//cout << AllHWPInRuns.size() << endl;
/*
	AllFCRunInfo HWPInData;
	HWPInData.SetRunInfo( RP, AllHWPInRuns );
	auto Plt_HWP_In = HWPInData.MakeAllFCPlots( "Unpol", "Summer" );
*/

	//AllData.Print();
/*
	//TGraphErrors* MakePlot( string Target, string Period, string BCA_Type="Gated" ){
	auto gr_EPICS_FCUP     = AllData.MakePlot( "C", "Summer", "EPICS_FCUP", 1 );
	auto gr_EPICS_FCUP_norm = AllData.MakePlot( "C", "Summer", "EPICS_FCUP", 2, true );

	TCanvas* c = new TCanvas("c","c",1600,600); c->Divide(2,1);
	c->cd(1);
	gr_EPICS_FCUP->Draw("ap");
	c->cd(2);
	gr_EPICS_FCUP_norm->Draw("ap");
*/
/*
	// Set the file extension based on whether or not the user wants to use QA cuts or not
	string fileExt, normExt;
	if( noQA ) fileExt = "_NoQA";
	if( normToReadouts ) normExt = "_NormR";
	else if( plotVsRun ) normExt = "_VsRun";

	// Get the count readout info...
	// Right now, there's only readout data for the NoQA option...
	vector<int> RunList, RunListNoQA;
	vector<double> HelN_R, HelP_R, HelN_R_NoQA, HelP_R_NoQA;
	ifstream finR( string("Scratch/Scaler_Readouts"+ fileExt +".txt") );
	string line; getline( finR, line ); // Throw away header row
	while( getline( finR, line ) ){
		stringstream sin(line);
		int run; double heln, help;
		sin >> run >> heln >> help;
		RunListNoQA.push_back( run );
		HelN_R_NoQA.push_back( heln );
		HelP_R_NoQA.push_back( help );
	}
	finR.close();
	

	RunPeriod Period; // Holds run info for all RGC runs...
	
	//UG_HelN   UG_HelP   UG_Hel0   G_HelN   G_HelP   G_Hel0   QA_HelN   QA_HelP   QA_Hel0   RUN_Hel   Elec_HelN   Elec_HelP   Elec_Hel0

	// Read in the FC charge data with electron counts...
	ifstream fin( string("All_FC_Info"+ fileExt +".txt") );
	cout << string("All_FC_Info"+ fileExt +".txt") << endl; 

	// Holds the runs for different data sets
	vector<RunFCData> Summer, FallNeg, FallPos, SpringInb, SpringOutb;
	
	// Loop over the input files.
	// The "UG" values are ungated, the "G" values are gated, the QA values can be ignored for now, the RUN_Hel is the RUN::scaler bank value,
	// and the "Elec_x" values are the counts for each helicity state that passed the DIS analysis cuts
	while( getline( fin, line ) ){
		stringstream sin(line);
		int Run; string Target;
		double UG_HelN, UG_HelP, UG_Hel0, G_HelN, G_HelP, G_Hel0, QA_HelN, QA_HelP, QA_Hel0, RUN_Hel, Elec_HelN, Elec_HelP, Elec_Hel0;
		sin >> Run >> Target >> UG_HelN >> UG_HelP >> UG_Hel0 >> G_HelN >> G_HelP >> G_Hel0 >> QA_HelN >> QA_HelP >> QA_Hel0 >> RUN_Hel >> Elec_HelN >> Elec_HelP >> Elec_Hel0;

		RunFCData thisRun; thisRun.SetValues( Run, Target, Elec_HelN, Elec_HelP, UG_HelN, UG_HelP, G_HelN, G_HelP );
		// This is really sloppy, but it just finds the number of scaler readouts corresponding to this run since it's stored in a separate file...
		int thisRunIndex = FindRunIndex( Run, RunListNoQA );
		thisRun.ReadP = HelN_R_NoQA[ thisRunIndex ];
		thisRun.ReadN = HelP_R_NoQA[ thisRunIndex ];

		thisRun.Print();
		if( Run < 16800 ) Summer.push_back( thisRun );
		else if( Run > 16800 && Run < 17185 ) FallNeg.push_back( thisRun );
		else if( Run > 17185 && Run < 17450 ) FallPos.push_back( thisRun );
		else if( Run > 17450 && Run <=17768 ) SpringInb.push_back( thisRun );
		else if( Run >= 17769 ) SpringOutb.push_back( thisRun );

	}
	fin.close();
*/
/*
	vector<TGraphErrors*> InP_ND3Summer  = MakeGraph( Period, Summer, 1, 1, Target, "Ungated", false );
	vector<TGraphErrors*> OutP_ND3Summer = MakeGraph( Period, Summer, 1, 0, Target, "Ungated", false );
	vector<TGraphErrors*> InN_ND3Summer  = MakeGraph( Period, Summer,-1, 1, Target, "Ungated", false );
	vector<TGraphErrors*> OutN_ND3Summer = MakeGraph( Period, Summer,-1, 0, Target, "Ungated", false );

	//TCanvas* ND3_Summer = FC_Asyms_vs_Runs( InP_ND3Summer, OutP_ND3Summer, InN_ND3Summer, OutN_ND3Summer, "Summer" );
	TCanvas* ND3_Summer = All_FC_Asyms( InP_ND3Summer, OutP_ND3Summer, InN_ND3Summer, OutN_ND3Summer, "Summer" );
*/
/*
	vector<TGraphErrors*> InP_C_Summer  = MakeGraph( Period, Summer, 0, 1, Target, "Ungated" );
	vector<TGraphErrors*> OutP_C_Summer = MakeGraph( Period, Summer, 0, 0, Target, "Ungated" );

	TCanvas* C_Summer = All_FC_Asyms( InP_C_Summer, OutP_C_Summer, "Summer" );
*/
/*

    // These functions plot the FC histograms for each of the run periods...
    TCanvas* Summer_Hist = Plot_FC_Histos( Summer, "Summer" );
    Summer_Hist->Print( string("PDF_Files/All_FC_Histos"+ fileExt +".pdf(").c_str() );

    TCanvas* FallNeg_Hist = Plot_FC_Histos( FallNeg, "FallNeg" );
    FallNeg_Hist->Print( string("PDF_Files/All_FC_Histos"+ fileExt +".pdf").c_str() );

    TCanvas* FallPos_Hist = Plot_FC_Histos( FallPos, "FallPos" );
    FallPos_Hist->Print( string("PDF_Files/All_FC_Histos"+ fileExt +".pdf").c_str() );

    TCanvas* SpringInb_Hist = Plot_FC_Histos( SpringInb, "SpringInb" );
    SpringInb_Hist->Print( string("PDF_Files/All_FC_Histos"+ fileExt +".pdf").c_str() );

    TCanvas* SpringOutb_Hist = Plot_FC_Histos( SpringOutb, "SpringOutb" );
    SpringOutb_Hist->Print( string("PDF_Files/All_FC_Histos"+ fileExt +".pdf)").c_str() );


    string outFile; // Used for writing the plots...

    // Plotting options
    //bool plotVsRun = false;
    //bool normToReadouts = true;

    ofstream fout( string(Target+"_ALU_vs_AFC_Vals"+ fileExt +".txt").c_str() );

    if( Target == "Unpol" ){
	vector<TGraphErrors*> InSummer  = MakeGraph( Period, Summer, 0, 1, Target, "Ungated", plotVsRun, normToReadouts, plotFCerrs );
	vector<TGraphErrors*> OutSummer = MakeGraph( Period, Summer, 0, 0, Target, "Ungated", plotVsRun, normToReadouts, plotFCerrs );
	TCanvas* SummerPlot = All_FC_Asyms( fout, InSummer, OutSummer, "Summer" );
	SummerPlot->Print( string("PDF_Files/Unpolarized_FC_Asyms"+ fileExt + normExt +".pdf(").c_str() );

	vector<TGraphErrors*> InFallNeg  = MakeGraph( Period, FallNeg, 0, 1, Target, "Ungated", plotVsRun, normToReadouts, plotFCerrs );
	vector<TGraphErrors*> OutFallNeg = MakeGraph( Period, FallNeg, 0, 0, Target, "Ungated", plotVsRun, normToReadouts, plotFCerrs );
	TCanvas* FallNegPlot = All_FC_Asyms( fout, InFallNeg, OutFallNeg, "FallNeg" );
	FallNegPlot->Print( string("PDF_Files/Unpolarized_FC_Asyms"+ fileExt + normExt +".pdf").c_str() );

	vector<TGraphErrors*> InFallPos  = MakeGraph( Period, FallPos, 0, 1, Target, "Ungated", plotVsRun, normToReadouts, plotFCerrs );
	vector<TGraphErrors*> OutFallPos = MakeGraph( Period, FallPos, 0, 0, Target, "Ungated", plotVsRun, normToReadouts, plotFCerrs );
	TCanvas* FallPosPlot = All_FC_Asyms( fout, InFallPos, OutFallPos, "FallPos" );
	FallPosPlot->Print( string("PDF_Files/Unpolarized_FC_Asyms"+ fileExt + normExt +".pdf").c_str() );

	vector<TGraphErrors*> InSpringInb  = MakeGraph( Period, SpringInb, 0, 1, Target, "Ungated", plotVsRun, normToReadouts, plotFCerrs );
	vector<TGraphErrors*> OutSpringInb = MakeGraph( Period, SpringInb, 0, 0, Target, "Ungated", plotVsRun, normToReadouts, plotFCerrs );
	TCanvas* SpringInbPlot = All_FC_Asyms( fout, InSpringInb, OutSpringInb, "SpringInb" );
	SpringInbPlot->Print( string("PDF_Files/Unpolarized_FC_Asyms"+ fileExt + normExt +".pdf").c_str() );

	vector<TGraphErrors*> InSpringOutb  = MakeGraph( Period, SpringOutb, 0, 1, Target, "Ungated", plotVsRun, normToReadouts, plotFCerrs );
	vector<TGraphErrors*> OutSpringOutb = MakeGraph( Period, SpringOutb, 0, 0, Target, "Ungated", plotVsRun, normToReadouts, plotFCerrs );
	TCanvas* SpringOutbPlot = All_FC_Asyms( fout, InSpringOutb, OutSpringOutb, "SpringOutb" );
	SpringOutbPlot->Print( string("PDF_Files/Unpolarized_FC_Asyms"+ fileExt + normExt +".pdf)").c_str() );


	vector<TGraphErrors*> InSummerG  = MakeGraph( Period, Summer, 0, 1, Target, "Gated", plotVsRun, normToReadouts, plotFCerrs );
	vector<TGraphErrors*> OutSummerG = MakeGraph( Period, Summer, 0, 0, Target, "Gated", plotVsRun, normToReadouts, plotFCerrs );
	TCanvas* SummerPlotG = All_FC_Asyms( fout, InSummerG, OutSummerG, "Summer", "Gated" );
	SummerPlotG->Print( string("PDF_Files/Unpolarized_GatedFC_Asyms"+ fileExt + normExt +".pdf(").c_str() );

	vector<TGraphErrors*> InFallNegG  = MakeGraph( Period, FallNeg, 0, 1, Target, "Gated", plotVsRun, normToReadouts, plotFCerrs );
	vector<TGraphErrors*> OutFallNegG = MakeGraph( Period, FallNeg, 0, 0, Target, "Gated", plotVsRun, normToReadouts, plotFCerrs );
	TCanvas* FallNegPlotG = All_FC_Asyms( fout, InFallNegG, OutFallNegG, "FallNeg", "Gated" );
	FallNegPlotG->Print( string("PDF_Files/Unpolarized_GatedFC_Asyms"+ fileExt + normExt +".pdf").c_str() );

	vector<TGraphErrors*> InFallPosG  = MakeGraph( Period, FallPos, 0, 1, Target, "Gated", plotVsRun, normToReadouts, plotFCerrs );
	vector<TGraphErrors*> OutFallPosG = MakeGraph( Period, FallPos, 0, 0, Target, "Gated", plotVsRun, normToReadouts, plotFCerrs );
	TCanvas* FallPosPlotG = All_FC_Asyms( fout, InFallPosG, OutFallPosG, "FallPos", "Gated" );
	FallPosPlotG->Print( string("PDF_Files/Unpolarized_GatedFC_Asyms"+ fileExt + normExt +".pdf").c_str() );

	vector<TGraphErrors*> InSpringInbG  = MakeGraph( Period, SpringInb, 0, 1, Target, "Gated", plotVsRun, normToReadouts, plotFCerrs );
	vector<TGraphErrors*> OutSpringInbG = MakeGraph( Period, SpringInb, 0, 0, Target, "Gated", plotVsRun, normToReadouts, plotFCerrs );
	TCanvas* SpringInbPlotG = All_FC_Asyms( fout, InSpringInbG, OutSpringInbG, "SpringInb", "Gated" );
	SpringInbPlotG->Print( string("PDF_Files/Unpolarized_GatedFC_Asyms"+ fileExt + normExt +".pdf").c_str() );

	vector<TGraphErrors*> InSpringOutbG  = MakeGraph( Period, SpringOutb, 0, 1, Target, "Gated", plotVsRun, normToReadouts, plotFCerrs );
	vector<TGraphErrors*> OutSpringOutbG = MakeGraph( Period, SpringOutb, 0, 0, Target, "Gated", plotVsRun, normToReadouts, plotFCerrs );
	TCanvas* SpringOutbPlotG = All_FC_Asyms( fout, InSpringOutbG, OutSpringOutbG, "SpringOutb", "Gated" );
	SpringOutbPlotG->Print( string("PDF_Files/Unpolarized_GatedFC_Asyms"+ fileExt + normExt +".pdf)").c_str() );
    }
    else{
	// Ungated values
	vector<TGraphErrors*> InPSummer  = MakeGraph( Period, Summer, 1, 1, Target, "Ungated", plotVsRun, normToReadouts, plotFCerrs );
	vector<TGraphErrors*> OutPSummer = MakeGraph( Period, Summer, 1, 0, Target, "Ungated", plotVsRun, normToReadouts, plotFCerrs );
	vector<TGraphErrors*> InNSummer  = MakeGraph( Period, Summer,-1, 1, Target, "Ungated", plotVsRun, normToReadouts, plotFCerrs );
	vector<TGraphErrors*> OutNSummer = MakeGraph( Period, Summer,-1, 0, Target, "Ungated", plotVsRun, normToReadouts, plotFCerrs );

	TCanvas* SummerPlot = All_FC_Asyms( fout, InPSummer, OutPSummer, InNSummer, OutNSummer, "Summer" );
	outFile = "PDF_Files/"+Target+"_FC_Asyms"+ fileExt + normExt +".pdf(";
	SummerPlot->Print( outFile.c_str() );

	vector<TGraphErrors*> InPFallNeg  = MakeGraph( Period, FallNeg, 1, 1, Target, "Ungated", plotVsRun, normToReadouts, plotFCerrs );
	vector<TGraphErrors*> OutPFallNeg = MakeGraph( Period, FallNeg, 1, 0, Target, "Ungated", plotVsRun, normToReadouts, plotFCerrs );
	vector<TGraphErrors*> InNFallNeg  = MakeGraph( Period, FallNeg,-1, 1, Target, "Ungated", plotVsRun, normToReadouts, plotFCerrs );
	vector<TGraphErrors*> OutNFallNeg = MakeGraph( Period, FallNeg,-1, 0, Target, "Ungated", plotVsRun, normToReadouts, plotFCerrs );

	TCanvas* FallNegPlot = All_FC_Asyms( fout, InPFallNeg, OutPFallNeg, InNFallNeg, OutNFallNeg, "FallNeg" );
	outFile = "PDF_Files/"+Target+"_FC_Asyms"+ fileExt + normExt +".pdf";
	FallNegPlot->Print( outFile.c_str() );

	vector<TGraphErrors*> InPFallPos  = MakeGraph( Period, FallPos, 1, 1, Target, "Ungated", plotVsRun, normToReadouts, plotFCerrs );
	vector<TGraphErrors*> OutPFallPos = MakeGraph( Period, FallPos, 1, 0, Target, "Ungated", plotVsRun, normToReadouts, plotFCerrs );
	vector<TGraphErrors*> InNFallPos  = MakeGraph( Period, FallPos,-1, 1, Target, "Ungated", plotVsRun, normToReadouts, plotFCerrs );
	vector<TGraphErrors*> OutNFallPos = MakeGraph( Period, FallPos,-1, 0, Target, "Ungated", plotVsRun, normToReadouts, plotFCerrs );

	TCanvas* FallPosPlot = All_FC_Asyms( fout, InPFallPos, OutPFallPos, InNFallPos, OutNFallPos, "FallPos" );
	outFile = "PDF_Files/"+Target+"_FC_Asyms"+ fileExt + normExt +".pdf";
	FallPosPlot->Print( outFile.c_str() );

	vector<TGraphErrors*> InPSpringInb  = MakeGraph( Period, SpringInb, 1, 1, Target, "Ungated", plotVsRun, normToReadouts, plotFCerrs );
	vector<TGraphErrors*> OutPSpringInb = MakeGraph( Period, SpringInb, 1, 0, Target, "Ungated", plotVsRun, normToReadouts, plotFCerrs );
	vector<TGraphErrors*> InNSpringInb  = MakeGraph( Period, SpringInb,-1, 1, Target, "Ungated", plotVsRun, normToReadouts, plotFCerrs );
	vector<TGraphErrors*> OutNSpringInb = MakeGraph( Period, SpringInb,-1, 0, Target, "Ungated", plotVsRun, normToReadouts, plotFCerrs );

	TCanvas* SpringInbPlot = All_FC_Asyms( fout, InPSpringInb, OutPSpringInb, InNSpringInb, OutNSpringInb, "SpringInb" );
	outFile = "PDF_Files/"+Target+"_FC_Asyms"+ fileExt + normExt +".pdf";
	SpringInbPlot->Print( outFile.c_str() );

	vector<TGraphErrors*> InPSpringOutb  = MakeGraph( Period, SpringOutb, 1, 1, Target, "Ungated", plotVsRun, normToReadouts, plotFCerrs );
	vector<TGraphErrors*> OutPSpringOutb = MakeGraph( Period, SpringOutb, 1, 0, Target, "Ungated", plotVsRun, normToReadouts, plotFCerrs );
	vector<TGraphErrors*> InNSpringOutb  = MakeGraph( Period, SpringOutb,-1, 1, Target, "Ungated", plotVsRun, normToReadouts, plotFCerrs );
	vector<TGraphErrors*> OutNSpringOutb = MakeGraph( Period, SpringOutb,-1, 0, Target, "Ungated", plotVsRun, normToReadouts, plotFCerrs );

	TCanvas* SpringOutbPlot = All_FC_Asyms( fout, InPSpringOutb, OutPSpringOutb, InNSpringOutb, OutNSpringOutb, "SpringOutb" );
	outFile = "PDF_Files/"+Target+"_FC_Asyms"+ fileExt + normExt +".pdf)";
	SpringOutbPlot->Print( outFile.c_str() );

	// Gated values
	vector<TGraphErrors*> InPSummerG   = MakeGraph( Period, Summer, 1, 1, Target, "Gated", plotVsRun, normToReadouts, plotFCerrs );
	vector<TGraphErrors*> OutPSummerG  = MakeGraph( Period, Summer, 1, 0, Target, "Gated", plotVsRun, normToReadouts, plotFCerrs );
	vector<TGraphErrors*> InNSummerG   = MakeGraph( Period, Summer,-1, 1, Target, "Gated", plotVsRun, normToReadouts, plotFCerrs );
	vector<TGraphErrors*> OutNSummerG  = MakeGraph( Period, Summer,-1, 0, Target, "Gated", plotVsRun, normToReadouts, plotFCerrs );

	TCanvas* SummerPlotG  = All_FC_Asyms( fout, InPSummerG , OutPSummerG , InNSummerG , OutNSummerG , "Summer", "Gated" );
	outFile = "PDF_Files/"+Target+"_GatedFC_Asyms"+ fileExt + normExt +".pdf(";
	SummerPlotG ->Print( outFile.c_str() );

	vector<TGraphErrors*> InPFallNegG   = MakeGraph( Period, FallNeg, 1, 1, Target, "Gated", plotVsRun, normToReadouts, plotFCerrs );
	vector<TGraphErrors*> OutPFallNegG  = MakeGraph( Period, FallNeg, 1, 0, Target, "Gated", plotVsRun, normToReadouts, plotFCerrs );
	vector<TGraphErrors*> InNFallNegG   = MakeGraph( Period, FallNeg,-1, 1, Target, "Gated", plotVsRun, normToReadouts, plotFCerrs );
	vector<TGraphErrors*> OutNFallNegG  = MakeGraph( Period, FallNeg,-1, 0, Target, "Gated", plotVsRun, normToReadouts, plotFCerrs );

	TCanvas* FallNegPlotG  = All_FC_Asyms( fout, InPFallNegG , OutPFallNegG , InNFallNegG , OutNFallNegG , "FallNeg", "Gated"  );
	outFile = "PDF_Files/"+Target+"_GatedFC_Asyms"+ fileExt + normExt +".pdf";
	FallNegPlotG ->Print( outFile.c_str() );

	vector<TGraphErrors*> InPFallPosG   = MakeGraph( Period, FallPos, 1, 1, Target, "Gated", plotVsRun, normToReadouts, plotFCerrs );
	vector<TGraphErrors*> OutPFallPosG  = MakeGraph( Period, FallPos, 1, 0, Target, "Gated", plotVsRun, normToReadouts, plotFCerrs );
	vector<TGraphErrors*> InNFallPosG   = MakeGraph( Period, FallPos,-1, 1, Target, "Gated", plotVsRun, normToReadouts, plotFCerrs );
	vector<TGraphErrors*> OutNFallPosG  = MakeGraph( Period, FallPos,-1, 0, Target, "Gated", plotVsRun, normToReadouts, plotFCerrs );

	TCanvas* FallPosPlotG  = All_FC_Asyms( fout, InPFallPosG , OutPFallPosG , InNFallPosG , OutNFallPosG , "FallPos", "Gated"  );
	outFile = "PDF_Files/"+Target+"_GatedFC_Asyms"+ fileExt + normExt +".pdf";
	FallPosPlotG ->Print( outFile.c_str() );

	vector<TGraphErrors*> InPSpringInbG   = MakeGraph( Period, SpringInb, 1, 1, Target, "Gated", plotVsRun, normToReadouts, plotFCerrs );
	vector<TGraphErrors*> OutPSpringInbG  = MakeGraph( Period, SpringInb, 1, 0, Target, "Gated", plotVsRun, normToReadouts, plotFCerrs );
	vector<TGraphErrors*> InNSpringInbG   = MakeGraph( Period, SpringInb,-1, 1, Target, "Gated", plotVsRun, normToReadouts, plotFCerrs );
	vector<TGraphErrors*> OutNSpringInbG  = MakeGraph( Period, SpringInb,-1, 0, Target, "Gated", plotVsRun, normToReadouts, plotFCerrs );

	TCanvas* SpringInbPlotG  = All_FC_Asyms( fout, InPSpringInbG , OutPSpringInbG , InNSpringInbG , OutNSpringInbG , "SpringInb", "Gated"  );
	outFile = "PDF_Files/"+Target+"_GatedFC_Asyms"+ fileExt + normExt +".pdf";
	SpringInbPlotG ->Print( outFile.c_str() );

	vector<TGraphErrors*> InPSpringOutbG   = MakeGraph( Period, SpringOutb, 1, 1, Target, "Gated", plotVsRun, normToReadouts, plotFCerrs );
	vector<TGraphErrors*> OutPSpringOutbG  = MakeGraph( Period, SpringOutb, 1, 0, Target, "Gated", plotVsRun, normToReadouts, plotFCerrs );
	vector<TGraphErrors*> InNSpringOutbG   = MakeGraph( Period, SpringOutb,-1, 1, Target, "Gated", plotVsRun, normToReadouts, plotFCerrs );
	vector<TGraphErrors*> OutNSpringOutbG  = MakeGraph( Period, SpringOutb,-1, 0, Target, "Gated", plotVsRun, normToReadouts, plotFCerrs );

	TCanvas* SpringOutbPlotG  = All_FC_Asyms( fout, InPSpringOutbG , OutPSpringOutbG , InNSpringOutbG , OutNSpringOutbG , "SpringOutb", "Gated"  );
	outFile = "PDF_Files/"+Target+"_GatedFC_Asyms"+ fileExt + normExt +".pdf)";
	SpringOutbPlotG ->Print( outFile.c_str() );

    }
*/

}
