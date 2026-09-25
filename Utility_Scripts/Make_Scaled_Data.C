/*************************************************************************************
 *
 * Author: Derek Holmberg
 * Date Created: 3/25/25
 * Last Modified: 9/18/26
 *
 * The purpose of this program is to create scaled versions of the CH2 target data
 * for the summer22 dataset. The purpose of the scaling is to create a pseudo-CD2
 * target set for the summer22 dataset. By multiplying the counts in each x, Q2 bin
 * by the ratio of CD2/CH2 counts in the corresponding bin, the pseudo dataset is
 * created.
 *
 * Why use the scaling instead of just plugging in the CD2 counts from the spring23
 * dataset? Using the scaling method accounts for any systematic differences between
 * the run periods, allowing for an "apples to apples" comparison within the sumemr22
 * run period.
 *
***************************************************************************************/

#include "../Header_Files/Binning_Classes.h"

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

void Make_Scaled_Data(){

    DataSet Scaling_Data;
/*
    string textFilePath;
    cout <<"Enter path to CD2 and CH2 directory: ";
    cin >> textFilePath; cout << endl;
*/
    string textFilePath = "../Latest_Skims/Text_Files/";

    // Track any missed runs
    vector<int> MissedRuns;

    // Used to set the kinematics for the targets
    vector<string> Targets = {};

    //FIXME: Evaluate which CH2, CD2, and ET runs are the most appropriate for the background runs

    // Loop for the spring CH2 data
    for(int Run : Sp23Inb_CH2_Runs){
	string path = textFilePath + "CH2_"+ to_string( Run ) +"_DF_Data.txt";
	if( !Scaling_Data.AddToBins( path, "CH2" ) ) MissedRuns.push_back( Run );
    }

    // Loop for the spring CD2 data
    for(int Run : Sp23Inb_CD2_Runs){
	string path = textFilePath + "CD2_"+ to_string( Run ) +"_DF_Data.txt";
	if( !Scaling_Data.AddToBins( path, "CD2" ) ) MissedRuns.push_back( Run );
    }

    // Get all the ET targets from across the data sets, not just spring.
    // There are no outbending ET runs
    for(int i=16100; i<=17768; i++){
	stringstream s1; s1 << textFilePath <<"/ET_"<< i <<"_DF_Data.txt";
	Scaling_Data.AddToBins( s1.str(), "ET" );
    }

    //Scaling_Data.Print();
   
    // Output file holding the scaling factor for each of the bins
    // NOTE: The N0 and FC0 counts are left unscaled for the time being
    ofstream fout("../Input_Text_Files/SummerAndFall22_CH2_Scaling_Factors.txt");
    //ofstream foutET("../Input_Text_Files/ET_Scaling_Factors.txt");
 
    for(size_t i=0; i<Q2_Bin_Bounds.size()-1; i++){ // Loop over Q2 bins
	double qmin = Q2_Bin_Bounds[i];
	double qmax = Q2_Bin_Bounds[i+1];
	double qmid = (qmin+qmax)/2.0;
	for(size_t j=0; j<X_Bin_Bounds.size()-1; j++){ // Loop over X bins
	    double xmin = X_Bin_Bounds[j];
	    double xmax = X_Bin_Bounds[j+1];
	    double xmid = (xmin+xmax)/2.0;
	    Bin thisBin = Scaling_Data.getThisBin(qmid, xmid);
	    // Get the total counts and total FC charges for all target types.
	    // The Nt and FCt values include only the helicity-latched counts and charges, not the
	    // undefined charges/counts like N0, FC0.
	    double Nt_CH2 = thisBin.getNt("CH2"); double FCt_CH2 = thisBin.getFCt("CH2");
	    double Nt_CD2 = thisBin.getNt("CD2"); double FCt_CD2 = thisBin.getFCt("CD2");
	    double Nt_ET  = thisBin.getNt("ET" ); double FCt_ET  = thisBin.getFCt("ET");
	    if( FCt_CH2 > 0 && FCt_CD2 > 0 && Nt_CH2 > 0 && Nt_CD2 > 0 ){
		// Statistical error on the scaling factor, not sure if this is needed...
		double nCD = Nt_CD2 / FCt_CD2; double nCH = Nt_CH2 / FCt_CH2;
		double f_error = sqrt( (nCD / (nCH*nCH *FCt_CD2)) + ((nCD*nCD) / (nCH*nCH*nCH*FCt_CH2)) );
		cout << "Scale factor for CD2/CH2 bin X = " << xmid <<", Q2 = "<< qmid <<" is: "<< (Nt_CD2/FCt_CD2)/(Nt_CH2/FCt_CH2) <<" +- "<< f_error << endl;
		fout << qmin <<"	"<< qmax <<"	"<< xmin <<"	"<< xmax <<"	"<<  (Nt_CD2/FCt_CD2)/(Nt_CH2/FCt_CH2) <<"	"<< f_error << endl;
/*
		double nET = Nt_ET / FCt_ET;
		double ET_Scale = ET_Scaling_Factors( xmid, q2mid );
		double et_ferror = sqrt(Nt_ET)/FCt_ET;
		cout << "Scale factor for ET bin X = " << xmid <<", Q2 = "<< qmid <<" is: "<< (Nt_CD2/FCt_CD2)/(Nt_CH2/FCt_CH2) <<" +- "<< f_error << endl;
		fout << qmin <<"	"<< qmax <<"	"<< xmin <<"	"<< xmax <<"	"<<  (Nt_CD2/FCt_CD2)/(Nt_CH2/FCt_CH2) <<"	"<< f_error << endl;
*/
	    }
	    else
		fout << qmin <<"	"<< qmax <<"	"<< xmin <<"	"<< xmax <<"	"<< 0 <<"	"<< 0 << endl;
	}
    }
    fout.close();

    // Now go through the process of scaling the Summer22 CH2 files by the new scaling factors
    for( int i=16298; i<=16303; i++){
	stringstream sout; sout << textFilePath << "CD2_"<<i<<"_DF_Data.txt";
	stringstream sscale; sscale << "../Input_Text_Files/SummerAndFall22_CH2_Scaling_Factors.txt";
	stringstream sch2; sch2 << textFilePath << "/CH2_"<<i<<"_DF_Data.txt";
	ofstream foutscale( sout.str().c_str() );
	ifstream finscale( sscale.str().c_str() );
	ifstream finch2( sch2.str().c_str() );
	// Both input files should have the same number of rows, so this should work lol
	string line1, line2;
	getline( finch2, line2 ); // Throw away header row
	foutscale << line2 << endl; // Rewrite the header row into new file
	while( getline( finch2, line2 ) ){
	    getline( finscale, line1 );
	    stringstream sin1(line1);
	    double qmin1, qmax1, xmin1, xmax1, scale;
	    sin1 >> qmin1 >> qmax1 >> xmin1 >> xmax1 >> scale;
	    stringstream sin2(line2);
	    double qmin2, qmax2, xmin2, xmax2, nm, np, fcm, fcp, n0, fc0;
	    sin2 >> qmin2 >> qmax2 >> xmin2 >> xmax2 >> nm >> np >> fcm >> fcp >> n0 >> fc0;
	    
	    //cout << qmin1<<"   "<< qmax1<<"   "<< xmin1<<"   "<< xmax1<<"   "<< endl;
	    //cout << qmin2<<"   "<< qmax2<<"   "<< xmin2<<"   "<< xmax2<<"   "<< nm<<"   "<< np<<"   "<< fcm<<"   "<< fcp << endl;

	    foutscale << qmin2 <<"	"<< qmax2 <<"	"<< xmin2 <<"	"<< xmax2 <<"	"<< scale*nm <<"	"<< scale*np <<"	"<< fcm <<"	"<< fcp;
	    foutscale <<"	"<< n0 <<"	"<< fc0 << endl;
	    if( qmin1 != qmin2 || qmax1 != qmax2 || xmin1 != xmin2 || xmax1 != xmax2 ) cout << "ERROR: Files in run "<<i<<" not the same size.\n";

	}
	foutscale.close(); finscale.close(); finch2.close();
    }

    // Now go through the process of scaling the Fall22 CH2 files by the new scaling factors
    for( int i=17118; i<=17408; i++){
	stringstream sout; sout << "../Input_Text_Files/Scaled_CH2_"<<i<<"_DF_Data.txt";
	stringstream sscale; sscale << "../Input_Text_Files/SummerAndFall22_CH2_Scaling_Factors.txt";
	stringstream sch2; sch2 << textFilePath << "/CH2_"<<i<<"_DF_Data.txt";
	ifstream finch2( sch2.str().c_str() );
	// Both input files should have the same number of rows, so this should work lol
	if( !finch2.fail() ){
	  ofstream foutscale( sout.str().c_str() );
 	  ifstream finscale( sscale.str().c_str() );
	  string line1, line2;
	  getline( finch2, line2 ); // Throw away header row
	  foutscale << line2 << endl; // Rewrite the header row into new file
	  while( getline( finch2, line2 ) ){
	    getline( finscale, line1 );
	    stringstream sin1(line1);
	    double qmin1, qmax1, xmin1, xmax1, scale;
	    sin1 >> qmin1 >> qmax1 >> xmin1 >> xmax1 >> scale;
	    stringstream sin2(line2);
	    double qmin2, qmax2, xmin2, xmax2, nm, np, fcm, fcp, n0, fc0;
	    sin2 >> qmin2 >> qmax2 >> xmin2 >> xmax2 >> nm >> np >> fcm >> fcp >> n0 >> fc0;
	    
	    //cout << qmin1<<"   "<< qmax1<<"   "<< xmin1<<"   "<< xmax1<<"   "<< endl;
	    //cout << qmin2<<"   "<< qmax2<<"   "<< xmin2<<"   "<< xmax2<<"   "<< nm<<"   "<< np<<"   "<< fcm<<"   "<< fcp << endl;

	    foutscale << qmin2 <<"	"<< qmax2 <<"	"<< xmin2 <<"	"<< xmax2 <<"	"<< scale*nm <<"	"<< scale*np <<"	"<< fcm <<"	"<< fcp;
	    foutscale <<"	"<< n0 <<"	"<< fc0 << endl;

	    if( qmin1 != qmin2 || qmax1 != qmax2 || xmin1 != xmin2 || xmax1 != xmax2 ) cout << "ERROR: Files in run "<<i<<" not the same size.\n";

	  }
	  foutscale.close(); finscale.close();
	}
	finch2.close();
    }


    cout << "All done!\n";

}
