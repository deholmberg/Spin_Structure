/***********************************************************************************
 *
 * Author: Derek Holmberg
 *
 * Date Created: 9/18/2026
 *
 * Last Modified: 9/20/2026
 *
 * Purpose:
 * This program calculates the scaling factors required for the dilution factor
 * calculations in "Analysis_Channel/Dilution_Factors.C". This program creates 
 * factors for every x, Q2 bin that scales empty (ET), foil (F), and CH2 counts
 * to match the expected behavior in different solenoid configurations as well
 * as scaling CH2 counts to match CD2 behavior for run periods where no CD2 data
 * is available.
 *
***********************************************************************************/

#include "../Header_Files/Binning_Classes.h"
#include "../Header_Files/Background_Runs.h"

using namespace std;

// This is just a dumb class to hold the scaling data
class ScaleBin{

  public:
    double Q2_Min = 0; double Q2_Max = 0;
    double X_Min = 0;  double X_Max = 0;
    double ScaleFactor   = 0; // For CH2 and CD2
    double ScaleFactorET = 0; // For ET targets
    
    void SetBin( double q2min, double q2max, double xmin, double xmax, double sf){
	Q2_Min = q2min; Q2_Max = q2max; X_Min = xmin; X_Max = xmax; ScaleFactor = sf;
    }

};

// This function calculates the statistical error on a scaling factor given the normalized 
// counts and FC charges. Assumes ratio of the form f = A / B
double F_Error( double nA, double nB, double fcA, double fcB ){

    if( nB > 0 && fcA > 0 && fcB > 0 ){
	return sqrt( (nA/(nB*nB*fcA)) + (nA*nA)/(nB*nB*nB*fcB)  );
    }

    return 0;

}

// Used to plot the scaling factors after calculation...
TGraphErrors* Plot(string filePath, int color, string title="Scaling Factors", bool useRadLenCorr=false){

	// Holds the scaling factors, their statistical errors, and the corresponding x values.
	// Multiple scaling factors for each x because there's a Q2 dependence as well.
	vector<double> XBins, ScaleF, ScaleFerr;

	string line; // Used for reading in data
	ifstream fin( filePath );
	while( getline( fin, line ) ){
		stringstream sin(line);
		double qmin, qmax, xmin, xmax, sf, sferr, qavg, xavg;
		sin >> qmin >> qmax >> xmin >> xmax >> sf >> sferr >> qavg >> xavg;
		if( useRadLenCorr ){
		    double corr = ET_Scale_Factor( xavg, qavg );
		    sf = sf*corr;
		    sferr = sferr*corr;
		}

		if( sf != 0 ){ // If it's zero, there's no data there
			//cout << xavg <<" "<< sf << endl;
			XBins.push_back( xavg );
			ScaleF.push_back( sf );
			ScaleFerr.push_back( sferr );
		}
	}
	fin.close();

	TGraphErrors* SF_Plt = new TGraphErrors( XBins.size(), XBins.data(), ScaleF.data(), nullptr, ScaleFerr.data() );
	SF_Plt->SetTitle( title.c_str() );
	SF_Plt->SetMarkerStyle( kFullCircle );
	SF_Plt->SetMarkerColor( color );
	SF_Plt->GetYaxis()->SetRangeUser(0.8,1.2);
	//SF_Plt->SetLineColor( color );
	// Calculate the weighted average, if applicable
	cout << "Results for "<< filePath <<":\n";
	TFitResultPtr r;
	r = SF_Plt->Fit("pol0","S","Q"); // Fit to a constant
	SF_Plt->GetFunction("pol0")->SetLineColor( color );

	//vgVal = r->Value(0); avgValErr = r->Error(0);
	//SF_Plt->GetYaxis()->SetRangeUser(0.8,1.2);

	return SF_Plt;

	//TCanvas* c = new TCanvas("c","c",800,600);
	//c->cd();
	//SF_Plt->Draw("ap");
	
}

void Create_Scaling_Factors(){

    DataSet ScalingPolyEthData;   // This is used to create the CD2/CH2 scaling factors
    DataSet ScalingSolenoidData;  // This is used to create the C_pos/C_neg and CH2_pos/CH2_neg ratios for the
				  // solenoid scaling factors.

    // These vectors are used to tell the DataSet objects which target kinematics to average over
    vector<string> Spring_PE_Targets = {"CH2", "CD2"};

    // Start with generating the count ratios for the spring CD2/CH2 data:
    for(int run : Sp23Inb_CD2_Runs){
	stringstream s1; s1 << "../"<< "Latest_Skims/Text_Files" <<"/CD2_"<< run <<"_DF_Data.txt";
	ScalingPolyEthData.AddToBins( s1.str(), "CD2" );
    }
    for(int run : Sp23Inb_CH2_Runs){
	stringstream s1; s1 << "../"<< "Latest_Skims/Text_Files" <<"/CH2_"<< run <<"_DF_Data.txt";
	ScalingPolyEthData.AddToBins( s1.str(), "CH2" );
    }
   
    // Calculates the average x, Q2 values across both CH2 and CD2 target types and sets them to 
    // the "All" category.
    ScalingPolyEthData.CalculateAvgXQ2( "Sp23Inb", Spring_PE_Targets );

    // Ok, this next part is a little confusing, so I need to explain the madness. The DataSet object "ScalingSolenoidData"
    // holds the information for four sets of targets: C_pos, C_neg, CH2_pos, CH2_neg. I will take the average of the
    // ratios C_pos/C_neg, CH2_pos/CH2_neg as the solenoid scaling factors to generate the pseudo-data for ET and F for
    // the fall positive solenoid run period. As such, I need to store all these counts into one DataSet object. However,
    // the DataSets only contain storage for each target type once, so I need to store the neg. solenoid data in other 
    // target types. Here's how they're stored:
    // --> C_pos   = "C"
    // --> C_neg   = "ET"
    // --> CH2_pos = "CH2"
    // --> CH2_neg = "CD2"
    // I hope that kinda sorta maybe makes sense and won't confuse me even more later lol. I'll use these strings to be less
    // confusing:
    string C_pos   = "C";
    string C_neg   = "ET";
    string CH2_pos = "CH2";
    string CH2_neg = "CD2";
    vector<string> Solenoid_Targets = { C_pos, C_neg, CH2_pos, CH2_neg };

    // Now set the solenoid scaling ratio data
    // C_pos data
    for(int run : Fa22Pos_C_Runs){
	string path = "../Latest_Skims/Text_Files/C_"+ to_string(run) +"_DF_Data.txt";
	ScalingSolenoidData.AddToBins( path, C_pos );
    }
    // C_neg data
    for(int run : Fa22Neg_C_Runs){
	string path = "../Latest_Skims/Text_Files/C_"+ to_string(run) +"_DF_Data.txt";
	ScalingSolenoidData.AddToBins( path, C_neg );
    }
    // CH2_pos data
    for(int run : Fa22Pos_CH2_Runs){
	string path = "../Latest_Skims/Text_Files/CH2_"+ to_string(run) +"_DF_Data.txt";
	ScalingSolenoidData.AddToBins( path, CH2_pos );
    }
    // CH2_neg data
    for(int run : Fa22Neg_CH2_Runs){
	string path = "../Latest_Skims/Text_Files/CH2_"+ to_string(run) +"_DF_Data.txt";
	ScalingSolenoidData.AddToBins( path, CH2_neg );
    }

    // Calculate the average x, Q2 for this data:
    ScalingSolenoidData.CalculateAvgXQ2( "Solenoid_Scale", Solenoid_Targets );

    // Output file holding the scaling factor for each of the bins
    // NOTE: The N0 and FC0 counts are left unscaled for the time being
    ofstream fout(   "../Input_Text_Files/SuFa22_CH2_Scaling_Factors.txt");
    ofstream foutC(  "../Input_Text_Files/Fall_C_Ratios.txt");
    ofstream foutCH2("../Input_Text_Files/Fall_CH2_Ratios.txt");
    ofstream foutAvg("../Input_Text_Files/Solenoid_Scaling_Factors.txt");

    for(size_t i=0; i<Q2_Bin_Bounds.size()-1; i++){ // Loop over Q2 bins
	double qmin = Q2_Bin_Bounds[i];
	double qmax = Q2_Bin_Bounds[i+1];
	double qmid = (qmin+qmax)/2.0;
	for(size_t j=0; j<X_Bin_Bounds.size()-1; j++){ // Loop over X bins
	    double xmin = X_Bin_Bounds[j];
	    double xmax = X_Bin_Bounds[j+1];
	    double xmid = (xmin+xmax)/2.0;
	    Bin thisBin = ScalingPolyEthData.getThisBin(qmid, xmid);
	    // Get the total counts and total FC charges for all target types.
	    // The Nt and FCt values include only the helicity-latched counts and charges, not the
	    // undefined charges/counts like N0, FC0.
	    double Nt_CH2 = thisBin.getNt("CH2"); double FCt_CH2 = thisBin.getFCt("CH2");
	    double Nt_CD2 = thisBin.getNt("CD2"); double FCt_CD2 = thisBin.getFCt("CD2");
	    if( FCt_CH2 > 0 && FCt_CD2 > 0 && Nt_CH2 > 0 && Nt_CD2 > 0 ){
		// Statistical error on the scaling factor, not sure if this is needed...
		double nCD = Nt_CD2 / FCt_CD2; double nCH = Nt_CH2 / FCt_CH2;
		//double f_error = sqrt( (nCD / (nCH*nCH *FCt_CD2)) + ((nCD*nCD) / (nCH*nCH*nCH*FCt_CH2)) );
		double f_error = F_Error( nCD, nCH, FCt_CD2, FCt_CH2 );

		//cout << "Scale factor for CD2/CH2 bin X = " << xmid <<", Q2 = "<< qmid <<" is: "<< (Nt_CD2/FCt_CD2)/(Nt_CH2/FCt_CH2) <<" +- "<< f_error << endl;
		fout << qmin <<"	"<< qmax <<"	"<< xmin <<"	"<< xmax <<"	"<<  (Nt_CD2/FCt_CD2)/(Nt_CH2/FCt_CH2) <<"	"<< f_error <<"	";
		fout << thisBin.getAvgQ2("All") <<"	"<< thisBin.getAvgX("All") << endl;
	    }
	    else
		fout << qmin <<"	"<< qmax <<"	"<< xmin <<"	"<< xmax <<"	"<< 0 <<"	"<< 0 <<"	"<< 0 <<"	"<< 0 << endl;

	    // Now do the same for the carbon and CH2 solenoid ratios and their average...
	    Bin thisBinSol = ScalingSolenoidData.getThisBin( qmid, xmid );
	    // Get the counts and FC charge info
	    double Nt_Cpos   = thisBinSol.getNt( C_pos );   double FCt_Cpos   = thisBinSol.getFCt( C_pos );
	    double Nt_Cneg   = thisBinSol.getNt( C_neg );   double FCt_Cneg   = thisBinSol.getFCt( C_neg );
	    double Nt_CH2pos = thisBinSol.getNt( CH2_pos ); double FCt_CH2pos = thisBinSol.getFCt( CH2_pos );
	    double Nt_CH2neg = thisBinSol.getNt( CH2_neg ); double FCt_CH2neg = thisBinSol.getFCt( CH2_neg );

	    if( Nt_Cneg > 0 && FCt_Cneg > 0 && Nt_CH2neg > 0 && FCt_CH2neg > 0 && FCt_Cpos > 0 && FCt_CH2pos > 0 ){
		double n_Cpos   = Nt_Cpos / FCt_Cpos;
		double n_Cneg   = Nt_Cneg / FCt_Cneg;
		double n_CH2pos = Nt_CH2pos / FCt_CH2pos;
		double n_CH2neg = Nt_CH2neg / FCt_CH2neg;
		// calculate the scale factors
		double f_c       = n_Cpos / n_Cneg;
		double f_c_err   = F_Error( n_Cpos, n_Cneg, FCt_Cpos, FCt_Cneg );
		double f_ch2     = n_CH2pos / n_CH2neg;
		double f_ch2_err = F_Error( n_CH2pos, n_CH2neg, FCt_CH2pos, FCt_CH2neg );

		double f_avg     = (f_c + f_ch2) / 2.0;
		double f_avg_err = 0.5*sqrt( f_c_err*f_c_err + f_ch2_err*f_ch2_err );

		// Carbon ratio data
		foutC << qmin <<"	"<< qmax <<"	"<< xmin <<"	"<< xmax <<"	"<< f_c <<"	"<< f_c_err <<"	";
		foutC << thisBinSol.getAvgQ2("All") <<"	"<< thisBinSol.getAvgX("All") << endl;
		// CH2 ratio data
		foutCH2 << qmin <<"	"<< qmax <<"	"<< xmin <<"	"<< xmax <<"	"<< f_ch2 <<"	"<< f_ch2_err <<"	";
		foutCH2 << thisBinSol.getAvgQ2("All") <<"	"<< thisBinSol.getAvgX("All") << endl;
		// Average ratio data
		foutAvg << qmin <<"	"<< qmax <<"	"<< xmin <<"	"<< xmax <<"	"<< f_avg <<"	"<< f_avg_err <<"	";
		foutAvg << thisBinSol.getAvgQ2("All") <<"	"<< thisBinSol.getAvgX("All") << endl;

	    }

	} // End loop on X bins
    } // End loop on Q2 bins

    fout.close();
    foutC.close();
    foutCH2.close();
    foutAvg.close();

    cout << "All done! Now make some plots to make sure they look ok. Saving to 'PDF_Plots/' directory...\n";

    int red = 632; int blue = 600; int violet = 880;

    TCanvas* c = new TCanvas("c","c",2200,600); c->Divide(3,1);    

    // Carbon ratio plot
    TGraphErrors* grC = Plot( "../Input_Text_Files/Fall_C_Ratios.txt", blue, "Fall C_{pos}/C_{neg} Ratio; X; f_{C}(x,Q^{2})" );
    c->cd(1);
    grC->Draw("ap");

    // CH2 ratio plot
    TGraphErrors* grCH2 = Plot( "../Input_Text_Files/Fall_CH2_Ratios.txt", red, "Fall CH2_{pos}/CH2_{neg} Ratio; X; f_{CH2}(x,Q^{2})" );
    c->cd(2);
    grCH2->Draw("ap");

    // Carbon ratio plot
    TGraphErrors* grAvg = Plot( "../Input_Text_Files/Solenoid_Scaling_Factors.txt", violet, "Average Solenoid Scaling Ratio; X; f_{avg}(x,Q^{2})" );
    c->cd(3);
    grAvg->Draw("ap");

    c->Print("PDF_Plots/SolScalePlot.pdf","pdf");

}
