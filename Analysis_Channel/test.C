#include "../Header_Files/Binning_Classes.h"
using namespace std;

vector<double> WeightedAverage( vector<double>& Vals, vector<double>& Errs ){

    vector<double> Avg = {0,0};

    // Kill attempt if vectors aren't same size
    int sizeVal = Vals.size();
    int sizeErr = Errs.size();

    if( sizeVal != sizeErr ){
	cout <<"ERROR: Value and error vectors aren't same size. No weighted average computed.\n";
	return Avg;
    }

    // Used for weighted averages
    double num = 0; double den = 0;
    for( int i=0; i<sizeVal; i++ ){
	
	if( Errs[i] != 0 ){
	    num += Vals[i] / pow( Errs[i], 2 );
	    den += 1.0 / pow( Errs[i], 2 );
	}
	else{
	    cout <<"ERROR: Invalid (zero) value for the error. No weighted average computed.\n";
	    return Avg;
	}

    }

    if( den > 0 ){
	Avg[0] = num / den;
	Avg[1] = sqrt( 1.0 / den );
	cout << "Average is "<< Avg[0] <<" +- "<< Avg[1] << endl;
    }
    else{
        cout <<"ERROR: No weighted average computed.\n";
    }

    return Avg;

}


void test(){

    vector<double> Su22 = {0.972588, 1.02842};
    vector<double> Su22Err = {0.0326239, 0.0354384};

    vector<double> Fa22Neg = {0.950429, 0.932961};
    vector<double> Fa22NegErr = {0.0342851, 0.0422452};

    vector<double> Fa22Pos = {0.855632, 0.947196};
    vector<double> Fa22PosErr = {0.0446983, 0.0539796};

    vector<double> Sp23 = {0.918562, 0.876295};
    vector<double> Sp23Err = {0.049903, 0.0527555};

    auto su22 = WeightedAverage( Su22, Su22Err );
    auto fa22neg = WeightedAverage( Fa22Neg, Fa22NegErr );
    auto fa22pos = WeightedAverage( Fa22Pos, Fa22PosErr );
    auto sp23 = WeightedAverage( Sp23, Sp23Err );

    // Check the constants
    cout <<"NH3 constants:\n";
    cout << "A = "<< A << endl;
    cout << "B = "<< B << endl;
    cout << "C = "<< C << endl;
    cout << "a = "<< a << endl;
    cout << "D = "<< D << endl;
    cout << "E = "<< E << endl;
    cout << "F = "<< F << endl;
    cout << "G = "<< G << endl;
    cout << endl;

    cout <<"ND3 constants:\n";
    cout << "pCD = "<< pCD << endl;
    cout << "AA = "<< AA << endl;
    cout << "BB = "<< BB << endl;
    cout << "CC = "<< CC << endl;
    cout << "aa = "<< aa << endl;
    cout << "DD = "<< DD << endl;
    cout << "EE = "<< EE << endl;
    cout << "FF = "<< FF << endl;
    cout << "GG = "<< GG << endl;

}
