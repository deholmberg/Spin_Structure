// This header file contains the information required to make the bins for
// the dilution factors in x, Q2.

#ifndef BINNING_CLASSES_H
#define BINNING_CLASSES_H

//#include "/home/dsizzle/JLab_Stuff/Dilution_Factors/Input_Text_Files/Particle_Masses.h"
#include "/w/hallb-scshelf2102/clas12/holmberg/Spin_Structure/Header_Files/Functions_RGC.h"
#include "/w/hallb-scshelf2102/clas12/holmberg/Spin_Structure/Header_Files/Background_Runs.h"

using namespace std;

const vector<int> Palette = { 632, 800, 416, 600, 840, 880, 900, 616, 820, 432, 920, 860 };

//string THIS_DIR = "/home/dsizzle/JLab_Stuff/Dilution_Factors/Input_Text_Files/";
string THIS_DIR = "/w/hallb-scshelf2102/clas12/holmberg/Spin_Structure/Input_Text_Files/";

// These are parameters used in the calculations of the dilution factors and packing fraction

double dLc= 0.005; // Contraction for carbon of 0.5%
double dPE= 0.0244; // CH2 contraction percentage (written as decimal)
double dL = 0.021; // The percentage of the length for CH2 and CD2 contraction, at 2.1%
double dBath = 0.019; // Percentage expansion for the PTFE (Teflon) bath
double dA = -0.05; // Variation on ammonia density
//double dPE_r = 0.0193; // Values from the paper
//double dPE_l = 0.0193;
double dPE_r = 0.0173 / 0.92;
double dPE_l = 0.0183 / 0.92;

double dC_r  = -0.0016;// + 0.0005;
double dC_l  = -0.0025;// + 0.0005;

double lC = 1.678 * (1.0 - dC_l); // cm (length of carbon target)
double LHe = 5.86 * (1.0 - dBath); // cm (length of LHe bath)
double Lcell = 5.0 *(1.0 - dBath); // cm (length of target cell; assumed same contraction as bath; bath uses PTFE, cell PCTFE)
double lCH = 3.18 * (1.0 - dPE_l); // cm (length of CH2 target)
double lCD = 2.686;// cm (length of CD2 target)
//double pC = 1.7926/(12.0); // mol/cm^3 (carbon target density)
//double pC = 1.7926/(12.0*pow((1.0-dLc),1) ); // mol/cm^3 (carbon target density)
double pC = 1.7926/(12.0*pow((1.0 - dC_r),2) * (1.0 - dC_l) ); // mol/cm^3 (carbon target density)

//double pCH = 0.9425/(14.027); // mol/cm^3 (CH2 target density)
double pCH = 0.9425/(14.027 * pow((1.0-dPE_r),2) * (1.0 - dPE_l) ); // mol/cm^3 (CH2 target density)

double pCD = 1.0979/(16.039); // mol/cm^3 (CD2 target density)
double pA = 0.867/(17.031); // mol/cm^3 (NH3 target density)
//double pA = (1.0 + dA)*0.867/(17.031); // mol/cm^3 (NH3 target density)

double pD = 1.007/(20.048); // mol/cm^3 (ND3 target density)

// These parameters are used in the calculation of the NH3 DF and PF values
double A = -lCH*pCH/(lC*pC);
double B = ( pC*lC*LHe + (pCH-pC)*lC*lCH - lCH*pCH*LHe )/( -lC*pC*LHe );
double C = (pCH-pC)*lCH/(LHe*pC);

double a = lC*lCH*(9*pA*pC - 2*pCH*(pA+3*pC));

double D = -2*lCH*pCH/(9*lC*pC);
double E = (9*lC*pC - 2*lCH*pCH - a/(LHe*pA))/(-9*lC*pC);
double F = -a / (9*lC*pC*LHe*pA);

double G = 6*lCH*pCH / (9*LHe*pA);

// These parameters are used in the calculation of the NH3 DF and PF values
double AA = -lCD*pCD/(lC*pC);
double BB = ( pC*lC*LHe + (pCD-pC)*lC*lCD - lCD*pCD*LHe )/( -lC*pC*LHe );
double CC = (pCD-pC)*lCD/(LHe*pC);

double aa = lC*lCD*(9*pD*pC - 2*pCD*(pD+3*pC));

double DD = -2*lCD*pCD/(9*lC*pC);
double EE = (9*lC*pC - 2*lCD*pCD - aa/(LHe*pD))/(-9*lC*pC);
double FF = -aa / (9*lC*pC*LHe*pD);

double GG = 6*lCD*pCD / (9*LHe*pD);


// Sets the beam energy based on the run number; neglects electron mass in beam energy
void SetBeamEnergy( TLorentzVector& BeamEnergy, int run ){
    if( run >= 16128 && run <= 17065 )      BeamEnergy.SetPxPyPzE( 0, 0, beam_energy1, beam_energy1 );
    else if( run >= 17067 && run <= 17704 ) BeamEnergy.SetPxPyPzE( 0, 0, beam_energy2, beam_energy2 );
    else if( run >= 17720 && run <= 17811 ) BeamEnergy.SetPxPyPzE( 0, 0, beam_energy3, beam_energy3 );
    else cout << "ERROR: Invalid run at "<< run <<". Beam energy not set. Get ready for some nonsense values lol.\n";
}

// Returns the true target polarization for a data set. Based on weighted averages of beam polarizations.
// First entry of returned vector is the Pt, the second is the error.
vector<double> CalculateTruePt( int run, double pbpt, double pbpt_stat_err ){

    // Holds calculated Pt values
    vector<double> TruePt = {0,0};

    // The summer beam polarization
    double pBeamSu22 = 0.8384;
    double pBeamSu22Err = sqrt( pow(0.0086,2) + pow(0.0216,2) ); // stat plus systematic
    // The fall beam polarization
    double pBeamFa22 = 0.8372;
    double pBeamFa22Err = sqrt( pow(0.0045,2) + pow(0.0291,2) );
    // The spring beam polarization
    double pBeamSp23 = 0.8040;
    double pBeamSp23Err = sqrt( pow(0.0061,2) + pow(0.0470,2) );

    if( run >= 16043 && run <= 16772 ){
	TruePt[0] = pbpt / pBeamSu22;
	TruePt[1] = sqrt( pow(pbpt_stat_err,2)/pow(pBeamSu22,2) + pow(pbpt,2)*pow(pBeamSu22Err,2)/pow(pBeamSu22,4) );
    }
    else if( run >= 16843 && run <= 17408 ){
	TruePt[0] = pbpt / pBeamFa22;
	TruePt[1] = sqrt( pow(pbpt_stat_err,2)/pow(pBeamFa22,2) + pow(pbpt,2)*pow(pBeamFa22Err,2)/pow(pBeamFa22,4) );
    }
    else if( run >= 17477 && run <= 17811 ){
	TruePt[0] = pbpt / pBeamSp23;
	TruePt[1] = sqrt( pow(pbpt_stat_err,2)/pow(pBeamSp23,2) + pow(pbpt,2)*pow(pBeamSp23Err,2)/pow(pBeamSp23,4) );
    }

    return TruePt;

}


// Calculates the NH3 packing fraction for each kinematic bin
double NH3PF(double nA, double nCH, double nC, double nET, double nF){

	if( nA < 0 || nCH < 0 || nET < 0 || nF < 0 || nC < 0 ) return 0;

	double numerator = G*(nA - nET);
	double denominator = nCH + D*nC + E*nET + F*nF;

	if( denominator != 0.0 ) return numerator / denominator;
	else{ 
	    //cout << "ERROR: Invalid NH3 packing fraction.\n";
	    //cout << "nA = "<< nA <<" nCH = "<< nCH <<" nC = "<< nC <<" nET = "<< nET <<" nF = "<< nF << endl;
	    return 0.0; // error case
	}
}
// Calculates the error in the NH3 PF for each kinematic bin
double NH3PFError(double nA, double nCH, double nC, double nET, double nF, double dnA, double dnCH, double dnC, double dnET, double dnF){

	if( nA < 0 || nCH < 0 || nET < 0 || nF < 0 || nC < 0 ) return 0;

	double denom = nCH + D*nC + E*nET + F*nF;
	if( denom != 0.0 ){
		double DPA = G / denom;
		double DPET = (-G/denom) + ( G*(nET - nA)*E / pow(denom,2) );
		double DPC = G*(nET - nA)*D / pow(denom,2);
		double DPCH = G*(nET - nA) / pow(denom,2);
		double DPF = G*(nET - nA)*F / pow(denom,2);
		return sqrt( dnA*dnA*DPA*DPA + dnET*dnET*DPET*DPET + dnC*dnC*DPC*DPC + dnCH*dnCH*DPCH*DPCH + dnF*dnF*DPF*DPF );
	}
	else return 0.0; // error case

}

// Calculates the ND3 packing fraction for each kinematic bin
double ND3PF(double nD, double nCD, double nC, double nET, double nF){

	if( nD < 0 || nCD < 0 || nET < 0 || nF < 0 || nC < 0 ) return 0;

	double numerator = GG*(nD - nET);
	double denominator = nCD + DD*nC + EE*nET + FF*nF;

	if( denominator != 0.0 ) return numerator / denominator;
	else{ 
	    //cout << "ERROR: Invalid NH3 packing fraction.\n";
	    //cout << "nD = "<< nD <<" nCD = "<< nCD <<" nC = "<< nC <<" nET = "<< nET <<" nF = "<< nF << endl;
	    return 0.0; // error case
	}
}
// Calculates the error in the ND3 PF for each kinematic bin
double ND3PFError(double nD, double nCD, double nC, double nET, double nF, double dnD, double dnCD, double dnC, double dnET, double dnF){

	if( nD < 0 || nCD < 0 || nET < 0 || nF < 0 || nC < 0 ) return 0;

	double denom = nCD + DD*nC + EE*nET + FF*nF;
	if( denom != 0.0 ){
		double DPD = GG / denom;
		double DPET = (-GG/denom) + ( GG*(nET - nD)*EE / pow(denom,2) );
		double DPC = GG*(nET - nD)*DD / pow(denom,2);
		double DPCD = GG*(nET - nD) / pow(denom,2);
		double DPF = GG*(nET - nD)*FF / pow(denom,2);
		return sqrt( dnD*dnD*DPD*DPD + dnET*dnET*DPET*DPET + dnC*dnC*DPC*DPC + dnCD*dnCD*DPCD*DPCD + dnF*dnF*DPF*DPF );
	}
	else return 0.0; // error case

}


// Calculates the dilution factor for NH3 targets using raw statistics only
double NH3DF(double nA, double nCH, double nC, double nET, double nF){

	if( nA < 0 || nCH < 0 || nET < 0 || nF < 0 || nC < 0 ) return 0;

	double numerator = (nA - nET)*(nCH + A*nC + B*nET + C*nF);
	double denominator = nA*(nCH + D*nC + E*nET + F*nF);

	if( denominator != 0.0 ) return numerator / denominator;
	else return 0; // error case
}

double NH3DFError(double nA, double nCH, double nC, double nET, double nF, double dnA, double dnCH, double dnC, double dnET, double dnF){

  if( nA < 0 || nCH < 0 || nET < 0 || nF < 0 || nC < 0 ) return 0;

  double numTerm = nCH + A*nC + B*nET + C*nF; // Not the full numerator, missing multiplication by (nA - nET)
  double denTerm = nCH + D*nC + E*nET + F*nF; // Not the full denominator, missing the multiplication by nA 
 
  if( denTerm != 0 && nA != 0 ){
     double DNA = numTerm / (nA*denTerm) + ( (nET - nA)*numTerm*denTerm / pow( (nA*denTerm), 2) );
     double DNET= -(nCH + A*nC +2*B*nET + C*nF) / (nA*denTerm) + ( (nET - nA)*numTerm*E*nA / pow( (nA*denTerm), 2) );
     double DNC = A*(nA - nET) / (nA*denTerm) + ( (nET - nA)*D*nA*numTerm / pow( (nA*denTerm), 2 ) );
     double DNCH= (nA - nET) / (nA*denTerm) + ( (nET - nA)*numTerm*nA / pow( (nA*denTerm), 2) );
     double DNF = C*(nA - nET) / (nA*denTerm) + ( (nET - nA)*numTerm*F*nA / pow( (nA*denTerm), 2) );
 
     return sqrt( dnA*dnA*DNA*DNA + dnET*dnET*DNET*DNET + dnC*dnC*DNC*DNC + dnCH*dnCH*DNCH*DNCH + dnF*dnF*DNF*DNF );
  }
  else return 0.0; // error case
  
}

// Calculates the dilution factor for ND3 targets using raw statistics only
double ND3DF(double nD, double nCD, double nC, double nET, double nF){

	if( nD <= 0 || nCD <= 0 || nET <= 0 || nF <= 0 || nC <= 0 ) return 0;

	double numerator = (nD - nET)*(nCD + AA*nC + BB*nET + CC*nF);
	double denominator = nD*(nCD + DD*nC + EE*nET + FF*nF);

	if( denominator != 0.0 ) return numerator / denominator;
	else return 0; // error case
}

double ND3DFError(double nD, double nCD, double nC, double nET, double nF, double dnD, double dnCD, double dnC, double dnET, double dnF){

  if( nD <= 0 || nCD <= 0 || nET <= 0 || nF <= 0 || nC <= 0 ) return 0;

  double numTerm = nCD + AA*nC + BB*nET + CC*nF; // Not the full numerator, missing multiplication by (nD - nET)
  double denTerm = nCD + DD*nC + EE*nET + FF*nF; // Not the full denominator, missing the multiplication by nD 
 
  if( denTerm != 0 && nD != 0 ){
     double DNA = numTerm / (nD*denTerm) + ( (nET - nD)*numTerm*denTerm / pow( (nD*denTerm), 2) );
     double DNET= -(nCD + AA*nC +2*BB*nET + CC*nF) / (nD*denTerm) + ( (nET - nD)*numTerm* EE *nD / pow( (nD*denTerm), 2) );
     double DNC = AA*(nD - nET) / (nD*denTerm) + ( (nET - nD)*DD*nD*numTerm / pow( (nD*denTerm), 2 ) );
     double DNCD= (nD - nET) / (nD*denTerm) + ( (nET - nD)*numTerm*nD / pow( (nD*denTerm), 2) );
     double DNF = CC*(nD - nET) / (nD*denTerm) + ( (nET - nD)*numTerm* FF *nD / pow( (nD*denTerm), 2) );
 
     return sqrt( dnD*dnD*DNA*DNA + dnET*dnET*DNET*DNET + dnC*dnC*DNC*DNC + dnCD*dnCD*DNCD*DNCD + dnF*dnF*DNF*DNF );
  }
  else return 0.0; // error case
  
}


// Calculates the dilution factor using a previously calculated value of the packing fraction (PF)
// The "useLinear" variable controls whether or not the user wants to use the old linear dependence
// on PF instead of the nonlinear relationship; not fully implemented yet...
double NH3DF_Thru_PF(double pf, double nA, double nCH, double nC, double nET, double nF, bool useLinear=false){
 
  if( nA <= 0 || nCH <= 0 || nET <= 0 || nF <= 0 || nC <= 0 ) return 0;
 
  double denom = nCH + D*nC + E*nET + F*nF;
  if( pf != 0.0 && denom != 0 && !useLinear ){

     return (nCH + A*nC + B*nET + C*nF) / (denom + nET*G/pf);

  }
  else return 0.0; // error case

}

// Calculates the error in the dilution factor calculated with the PF as input; does NOT include error
// propagation from the PF right now
// nA isn't used, but I left it in for legacy purposes to not break the code lol
double NH3DF_Thru_PF_Error(double pf, double nA, double nCH, double nC, double nET, double nF, double dpf, double dnA, double dnCH, double dnC, double dnET, double dnF, bool useLinear=false){

  if( nA <= 0 || nCH <= 0 || nET <= 0 || nF <= 0 || nC <= 0 ) return 0;

  double den = nCH + D*nC + F*nF;
  double num = nCH + A*nC + B*nET + C*nF;

  if( pf != 0 && den != 0 ){
     // These variables collect terms
     double M = den + nET*(E + G/pf);
     double DFPF = num * G * nET / (pf*pf*M*M);
     double DFCH = (1.0/M) - ( num /(M*M) );
     double DFC  = (A/M) - ( num * D /(M*M) );
     double DFET = (B/M) - ( num * (E + G/pf) / (M*M) );
     double DFF  = (C/M) - ( num * F /(M*M) );

     return sqrt( dnET*dnET*DFET*DFET + dnC*dnC*DFC*DFC + dnCH*dnCH*DFCH*DFCH + dnF*dnF*DFF*DFF + dpf*dpf*DFF*DFF );
  }
  else return 0.0;

}

// Calculates the dilution factor using a previously calculated value of the packing fraction (PF)
// The "useLinear" variable controls whether or not the user wants to use the old linear dependence
// on PF instead of the nonlinear relationship; not fully implemented yet...
double ND3DF_Thru_PF(double pf, double nD, double nCD, double nC, double nET, double nF, bool useLinear=false){
  
  if( nD <= 0 || nCD <= 0 || nET <= 0 || nF <= 0 || nC <= 0 ) return 0;

  double denom = nCD + DD*nC + EE*nET + FF*nF;
  if( pf != 0.0 && denom != 0 && !useLinear ){

     return (nCD + AA*nC + BB*nET + CC*nF) / (denom + nET*GG/pf);

  }
  else return 0.0; // error case

}

// Calculates the error in the dilution factor calculated with the PF as input; does NOT include error
// propagation from the PF right now
// nD isn't used, but I left it in for legacy purposes to not break the code lol
double ND3DF_Thru_PF_Error(double pf, double nD, double nCD, double nC, double nET, double nF, double dpf, double dnD, double dnCD, double dnC, double dnET, double dnF, bool useLinear=false){

  if( nD <= 0 || nCD <= 0 || nET <= 0 || nF <= 0 || nC <= 0 ) return 0;

  double den = nCD + DD*nC + FF*nF;
  double num = nCD + AA*nC + BB*nET + CC*nF;

  if( pf != 0 && den != 0 ){
     // These variables collect terms
     double MM = den + nET*(EE + GG/pf);
     double DFPF = num * GG * nET / (pf*pf*MM*MM);
     double DFCD = (1.0/MM) - ( num /(MM*MM) );
     double DFC  = (AA/MM) - ( num * DD /(MM*MM) );
     double DFET = (BB/MM) - ( num * (EE + GG/pf) / (MM*MM) );
     double DFF  = (CC/MM) - ( num * FF /(MM*MM) );

     return sqrt( dnET*dnET*DFET*DFET + dnC*dnC*DFC*DFC + dnCD*dnCD*DFCD*DFCD + dnF*dnF*DFF*DFF + dpf*dpf*DFF*DFF );
  }
  else return 0.0;

}

// This function scales the number of counts coming from empty, liquid helium target counts to match the helium
// counts that come from the helium in NH3 targets. This is based on taking the ratio of helium counts in NH3 to
// helium counts in empty targets, generated using Darren's model.
// Evaluated at a given bin in Bjorken X and Q2
double ET_Scale_Factor(double x, double q2){

 if( q2 > 0.0 && x > 0.0 ){
  double a = 30.0 / sqrt(q2);
  double b = -0.058 - 0.00004*pow(q2,3.5);
  double c = 0.825 + 0.12 / pow(q2,0.45);

  return exp(-a*(x+b)) + c;
 }
 //cout << "ERROR: Invalid x, Q2 passed into ET scale factor calculation. Returning scale of zero...\n";
 return 0.0;

}


//Run_Number  Target  Beam_Cur_Req  Beam_Energy  Num_Events  HWP_Status  Target_Pol  Solenoid_Scale  Torus_Scale
// This class holds information about each of the runs from RCDB
class Run{
    private:
	int RunNumber = 0; 
	string TargetType;
	string BeamCurrent = "N/A"; // in nC
	double RMSBeamCurrent = 0;
	double BeamEnergy = 0; // in GeV
        int NumberEvents = 0;
	int HWPStatus = -1; // 1=in, 0=out, -1=uninitialized
	double TargetPolarization = 0; // NMR Tpol for NH3 and ND3; zero for other targets
	double OfflineTPol = 0; // From Ishara's analysis
	double OfflineTPolErr = 0;
	double SolenoidScale = 0; // Current direction in solenoid
	double TorusScale = 0; // Current direction in torus
	string Epoch; // Run range 'epoch' that this run is a part of; used for DF/PF calculation
	int EpochNum; // Number of the epoch saved as an integer
    public:
	// Constructor
	Run() = default;
	// Mutators
	void SetRunInfo(int rn, string tt, double bc, double be, int ne, int hwp, double tp, double ss, double ts, string epoch){
	    RunNumber = rn; TargetType = tt; BeamCurrent = bc; BeamEnergy = be; NumberEvents = ne;
	    HWPStatus = hwp; TargetPolarization = tp; SolenoidScale = ss; TorusScale = ts; Epoch = epoch;
	    int length = epoch.length();
	    string epochnum; epochnum += epoch[length-2]; epochnum += epoch[length-1];
	    //cout << epoch <<" "<< epochnum <<" "<< rn <<"\n";
	    EpochNum = stoi(epochnum); //cout << EpochNum << endl;
	}
	void SetOffileTPol(double tpol, double tpolerr){
	    OfflineTPol = tpol; OfflineTPolErr = tpolerr;
	}
	void SetRMSBeamCurrent(double rms){
	    RMSBeamCurrent = rms;
	}
	// Accessors
	int getRunNumber() const{ return RunNumber;}
	string getTargetType() const{ return TargetType;}
	string getBeamCurrent() const{ return BeamCurrent;}
	double getBeamEnergy() const{ return BeamEnergy;}
	int getNumberEvents() const{ return NumberEvents;}
	int getHWPStatus() const{ return HWPStatus;}
	double getTargetPolarization() const{ return TargetPolarization;}
	double getOfflineTPol() const{ return OfflineTPol;}
	double getOfflineTPolErr() const{ return OfflineTPolErr;}
	double getSolenoidScale() const{ return SolenoidScale;}
	double getTorusScale() const{ return TorusScale;}
	string getEpoch() const{ return Epoch;}
	int getEpochNum() const{ return EpochNum;}
	double getRMSBeamCurrent() const{ return RMSBeamCurrent; }
	void Print(){
	  cout <<RunNumber<<" "<<TargetType<<" "<<BeamCurrent<<" "<<BeamEnergy<<" "<<NumberEvents<<" ";
	  cout <<HWPStatus<<" "<<TargetPolarization<<" "<<SolenoidScale<<" "<<TorusScale<<" "<<Epoch<<" "<< RMSBeamCurrent <<endl;
	}
	// Destructor
	~Run() =  default;
};

// This vector holds non-empty target runs where the misc bit should be ignored for the dilution factor analysis; whether
// or not a misc bit should be applied depends on your analysis channel.
const vector<int> IgnoreMiscBit = 
{
  // Summer 2022 Runs
  16178, // Higher beam current 
  16194, // Good foils run
  16243, // double rate from FT
  16309, // Empty target w/ LHe
  16733, 16734, 16736, // RICH sector 1 off
  
  // Fall 2022 Runs
  16872, 16970, 16975, // Empty w/ LHe
  16976, 16978, 16979, // Foils
  17003, // BMT S1 L4 shows abnormal distribution
  17069, // rich4 ROC issue
  17086, // varied beam current
  17088, // Ran briefly at 8nA, then went back to 6nA
  17164, // N/F looks stable, though it ran at lower beam current
  17179, // RICH errors
  17180, 17181, 17182, 17183, // RICH off
  17188, 17189, // RICH sector 1 partially down
  17191, // RICH sector 1 partially recovered during run
  17211, // BST multiplicity issue
  17215, // Ran briefly at 4nA, then went to 8nA
  17334, // Ran briefly at 4nA, then went to 8nA (excluded anyway)
  17336, // Beam current changed between 5nA and 8nA (excluded anyway)
  
  // Spring 2023 Runs
  17482, // Foils
  17535, 17537, 17538, 17540, 17541, // ND3 runs with CND issue, fine for us
  17544, 17545, 17546, 17547, 17548, 17549, 17550, 17551, 17552, 17553, 17554, 17556, 17557,
  17581, 17583, // BMT issue
  17593, 17594, 17595, 17596, 17597, 17598, 17599, 17600, 17601, // Fluctuating N/F
  17602, 17603, 17604, 17605, 17606, 17607, 17608, 17609, 17610, 17611, // Fluctuating N/F
  17696, // S2 rise in N/F
  17763, 17764, 17765, // Foils
  17766, 17767, 17768, // Empty w/ LHe
  17801 // Hole in PCAL 2, but N/F looked ok
  
};

// This class holds the run for a given run period
class RunPeriod{
    private:
	vector<Run> Runs;
    public:
	// Constructor
	RunPeriod(){ SetRunPeriod(); }
	// Mutator
	void SetRunPeriod(){
	    // Read in the information from the RMS beam charge values.
	    /*
	    string rmsPath = THIS_DIR + "RMS_Outputs_RGC.txt";
	    ifstream rmsIn( rmsPath.c_str() );
	    vector<double> rmsRuns, Currents;
	    if( !rmsIn.fail() ){
		string line;
		getline(rmsIn, line); // Throw away header row
		while( getline(rmsIn, line) ){
		    stringstream sin(line); double run, beamcurrent;
		    sin >> run >> beamcurrent;
		    //cout << run <<" "<< beamcurrent << endl;
		    rmsRuns.push_back(run);
		    Currents.push_back( beamcurrent );
		}
	    }
	    else cout <<"Failed to open RMS Beam current values...\n";
	    rmsIn.close();
	    */
	    // Read in the RGC run information
	    string fPath = THIS_DIR + "RGC_Run_Info.txt";
	    ifstream fin( fPath.c_str() );
	    if( !fin.fail() ){
	        string line;
	        //getline(fin,line); // throw away header row
	        while(getline(fin,line)){
	            stringstream sin(line);
		    int rn; string tt; double bc; double be; int ne; int hwp; double tp; double ss; double ts; string epoch;
		    sin >> rn >> tt >> bc >> be >> ne >> hwp >> tp >> ss >> ts >> epoch;
		    Run thisRun;
		    thisRun.SetRunInfo(rn, tt, bc, be, ne, hwp, tp, ss, ts, epoch);
		    // Now find the appropriate beam current to set...
		    /*
		    double rmsBC = 0; // RMS beam current
		    for(int k=0; k<rmsRuns.size(); k++){
			if( rmsRuns[k] == rn ){ rmsBC = Currents[k]; }
		    }
		    //cout << "Set run "<< rn <<" RMS beam current "<< rmsBC << endl;
		    thisRun.SetRMSBeamCurrent( rmsBC );
		    */
		    Runs.push_back(thisRun);
		}
	    }
	    else cout <<"Couldn't find the 'RGC_Run_Info.txt' input file; run period not set.\n";
	    fin.close();
	    //string fPath2 = THIS_DIR + "Offline_NMR_Values.txt";
	    string fPath2 = THIS_DIR + "NMR_Polarizations_Ishara.txt";
	    ifstream fin2( fPath2.c_str() );
	    if( !fin2.fail() ){
		string line; getline( fin2, line ); // throw away header row
		while(getline(fin2,line)){
		    stringstream sin(line);
		    int run, cell; string species; double avgOnline, avgOffline, avgOfflineErr, runDone;
		    sin >> run >> species >> cell >> avgOnline >> avgOffline >> avgOfflineErr >> runDone;
		    // Now slot into the appropriate run
		    for(size_t i=0; i<Runs.size(); i++){
			if(Runs[i].getRunNumber() == run){
			    Runs[i].SetOffileTPol( avgOffline, avgOfflineErr );
			    i = Runs.size();
			}
		    }
		}
	    }
	    else cout <<"Couldn't find the 'NMR_Polarizations_Ishara.txt' input file; offline NMR target polarizations not set.\n";
	    fin2.close();
	}
	// Accessors

	double getSolenoidScale(int runnum) const{
	    for(size_t i=0; i<Runs.size(); i++){
		if( Runs[i].getRunNumber() == runnum ){
		    return Runs[i].getSolenoidScale();
		}
	    }
	    cout <<"Couldn't find SolenoidScale for run "<< runnum <<". Check inputs.\n";
	    return 0;
	}	
	double getTorusScale(int runnum) const{
	    for(size_t i=0; i<Runs.size(); i++){
		if( Runs[i].getRunNumber() == runnum ){
		    return Runs[i].getTorusScale();
		}
	    }
	    cout <<"Couldn't find TorusScale for run "<< runnum <<". Check inputs.\n";
	    return 0;
	}
	double getTargetPolarization(int runnum) const{
	    for(size_t i=0; i<Runs.size(); i++){
		if( Runs[i].getRunNumber() == runnum ){
		    if( Runs[i].getTargetType() == "NH3" || Runs[i].getTargetType() == "ND3" )
		       	return Runs[i].getTargetPolarization();
		    else return 0;
		}
	    }
	    cout <<"Couldn't find run "<< runnum <<". Check inputs.\n";
	    return 0;
	}
	double getOfflineTPol(int runnum) const{
	    for(size_t i=0; i<Runs.size(); i++){
		if( Runs[i].getRunNumber() == runnum ){
		    if( Runs[i].getTargetType() == "NH3" || Runs[i].getTargetType() == "ND3" )
		       	return Runs[i].getOfflineTPol();
		    else return 0;
		}
	    }
	    cout <<"Couldn't find run "<< runnum <<". Check inputs.\n";
	    return 0;
	}
	double getOfflineTPolErr(int runnum) const{
	    for(size_t i=0; i<Runs.size(); i++){
		if( Runs[i].getRunNumber() == runnum ){
		    if( Runs[i].getTargetType() == "NH3" || Runs[i].getTargetType() == "ND3" )
		       	return Runs[i].getOfflineTPolErr();
		    else return 0;
		}
	    }
	    cout <<"Couldn't find run "<< runnum <<". Check inputs.\n";
	    return 0;
	}

	int getEpochNum(int run) const{
	    for(size_t i=0; i<Runs.size(); i++){
		if( Runs[i].getRunNumber() == run ) return Runs[i].getEpochNum();
	    }
	    return -1;
	}
	int getNumberEvents(int run) const{
	    for(size_t i=0; i<Runs.size(); i++){
		if( Runs[i].getRunNumber() == run ) return Runs[i].getNumberEvents();
	    }
	    return -1;
	}

	string getTargetType(int run) const{ 
	    for(size_t i=0; i<Runs.size(); i++){
		if( Runs[i].getRunNumber() == run ) return Runs[i].getTargetType();
	    }
	    cout <<"Couldn't find run "<< run <<". Check inputs.\n";
	    return "N/A";
	}
	string getEpoch(int run) const{ 
	    for(size_t i=0; i<Runs.size(); i++){
		if( Runs[i].getRunNumber() == run ) return Runs[i].getEpoch();
	    }
	    cout <<"Couldn't find the epoch for "<< run <<". Check inputs.\n";
	    return "N/A";
	}
	int getHWPStatus(int run) const{ 
	    for(size_t i=0; i<Runs.size(); i++){
		if( Runs[i].getRunNumber() == run ) return Runs[i].getHWPStatus();
	    }
	    cout <<"Couldn't find run "<< run <<". Check inputs. Returned HWP error status -1.\n";
	    return -1;
	}
	double getRMSBeamCurrent(int run) const{ 
	    for(size_t i=0; i<Runs.size(); i++){
		if( Runs[i].getRunNumber() == run ) return Runs[i].getRMSBeamCurrent();
	    }
	    cout <<"Couldn't find run "<< run <<". Check inputs. Returned beam current = -1.\n";
	    return -1;
	}

	// If the target is empty or foils, then ignore the "misc" error bit. Else, include the bit
	bool applyMiscBit(int run) const{
	/*
	    for(size_t i=0; i<Runs.size(); i++){
		if( Runs[i].getRunNumber() == run && Runs[i].getTargetType() == "Empty" ) return false;
		else if( Runs[i].getRunNumber() == run && Runs[i].getTargetType() == "Foil" ) return false;
	    }
	*/
	    for(int misc : IgnoreMiscBit ) if( misc == run ) return false;
	    return true;
	}
	vector<int> getRunEpoch(string epoch) const{
	    vector<int> RunsInEpoch;
	    for(size_t i=0; i<Runs.size(); i++){
		if( Runs[i].getEpoch() == epoch ) RunsInEpoch.push_back( Runs[i].getRunNumber() );
	    }
	    if( RunsInEpoch.size() == 0 ) cout << "Couldn't find runs for epoch "<<epoch<<endl;
	    return RunsInEpoch;
	}

	// This is used to get all the runs used in the elastic epochs. For example, it selects all positive
	// or negative runs in a given run period
	vector<int> getElasticEpoch(string Period, string Target, int targetPol) const{

	    // By selecting the sign of the target polarization, you can select the subset of runs from the
	    // epoch that correspond to the sign of the target polarization
	    vector<int> ElasticEpoch;

	    if( abs(targetPol) != 1 ){
		cout <<"ERROR: Proper usage of 'getElasticEpoch' is targetPol = +\\- 1. Check inputs.\n";
		return ElasticEpoch;
	    }

	    for( auto run : Runs ){
		string targ = run.getTargetType();
		int thisRun = run.getRunNumber();
	        double thisTPol = run.getTargetPolarization();
		if( Period == "Su22" && thisRun > 16100 && thisRun < 16800 && thisTPol * 10*targetPol > 0 && targ == Target ) ElasticEpoch.push_back( thisRun );
		else if( Period == "Fa22Neg" && thisRun > 16800 && thisRun < 17185 && thisTPol * 10*targetPol > 0 && targ == Target ) ElasticEpoch.push_back( thisRun );
		else if( Period == "Fa22Pos" && thisRun > 17185 && thisRun < 17450 && thisTPol * 10*targetPol > 0 && targ == Target ) ElasticEpoch.push_back( thisRun );
		else if( Period == "Sp23Inb" && thisRun > 17450 && thisRun <=17768 && thisTPol * 10*targetPol > 0 && targ == Target ) ElasticEpoch.push_back( thisRun );
	    }
	    if( ElasticEpoch.size() == 0 ) cout << "Couldn't find runs for run period "<< Period <<" with P_target ~ "<< targetPol << endl;
	    return ElasticEpoch;
	}

	vector<int> getDataSet(string Dataset, string TorPol, string SolPol) const{
	    // Returns the subset of runs consisting of this dataset, torus polarity, solenoid polarity
	    vector<int> RunsInDataSet;
	    for( auto run : Runs ){
		if( Dataset == "Su22" && TorPol == "Neg" && run.getTorusScale() < 0 && SolPol == "Neg" &&
		    run.getSolenoidScale() < 0 && run.getRunNumber() >= 16128 && run.getRunNumber() <= 16772 ){
		    RunsInDataSet.push_back( run.getRunNumber() );
		}
		else if( Dataset == "Fa22" && TorPol == "Neg" && run.getTorusScale() < 0 && SolPol == "Neg" &&
		    run.getSolenoidScale() < 0 && run.getRunNumber() >= 16843 && run.getRunNumber() <= 17408 ){
		    RunsInDataSet.push_back( run.getRunNumber() );
		}
		else if( Dataset == "Fa22" && TorPol == "Neg" && run.getTorusScale() < 0 && SolPol == "Pos" &&
		    run.getSolenoidScale() > 0 && run.getRunNumber() >= 16843 && run.getRunNumber() <= 17408 ){
		    RunsInDataSet.push_back( run.getRunNumber() );
		}
		// To return all the data in the fall data set
		else if( Dataset == "Fa22" && TorPol == "All" && SolPol == "All" && run.getRunNumber() >= 16843 && run.getRunNumber() <= 17408 ){
		    RunsInDataSet.push_back( run.getRunNumber() );
		}
		else if( Dataset == "Sp23" && TorPol == "Neg" && run.getTorusScale() < 0 && SolPol == "Neg" &&
		    run.getSolenoidScale() < 0 && run.getRunNumber() >= 17482 && run.getRunNumber() <= 17811 ){
		    RunsInDataSet.push_back( run.getRunNumber() );
		}
		else if( Dataset == "Sp23" && TorPol == "Pos" && run.getTorusScale() > 0 && SolPol == "Neg" &&
		    run.getSolenoidScale() < 0 && run.getRunNumber() >= 17482 && run.getRunNumber() <= 17811 ){
		    RunsInDataSet.push_back( run.getRunNumber() );
		}
		else if( Dataset == "Sp23" && TorPol == "All" && SolPol == "All" && run.getRunNumber() >= 17482 && run.getRunNumber() <= 17811 ){
		    RunsInDataSet.push_back( run.getRunNumber() );
		}
	    }
	    return RunsInDataSet;
	}
	vector<int> getDataSetForTarget(string Dataset, string Targ, string TorPol, string SolPol) const{
	    // Returns the subset of runs consisting of this dataset, torus polarity, solenoid polarity
	    vector<int> RunsInDataSet;
	    for( auto run : Runs ){
		if( Dataset == "Su22" && TorPol == "Neg" && run.getTorusScale() < 0 && SolPol == "Neg" &&
		    run.getSolenoidScale() < 0 && run.getRunNumber() >= 16134 && run.getRunNumber() <= 16772 && run.getTargetType() == Targ ){
		    RunsInDataSet.push_back( run.getRunNumber() );
		}
		else if( Dataset == "Fa22" && TorPol == "Neg" && run.getTorusScale() < 0 && SolPol == "Neg" &&
		    run.getSolenoidScale() < 0 && run.getRunNumber() >= 16843 && run.getRunNumber() <= 17408 && run.getTargetType() == Targ ){
		    RunsInDataSet.push_back( run.getRunNumber() );
		}
		else if( Dataset == "Fa22" && TorPol == "Neg" && run.getTorusScale() < 0 && SolPol == "Pos" &&
		    run.getSolenoidScale() > 0 && run.getRunNumber() >= 16843 && run.getRunNumber() <= 17408 && run.getTargetType() == Targ ){
		    RunsInDataSet.push_back( run.getRunNumber() );
		}
		else if( Dataset == "Sp23" && TorPol == "Neg" && run.getTorusScale() < 0 && SolPol == "Neg" &&
		    run.getSolenoidScale() < 0 && run.getRunNumber() >= 17482 && run.getRunNumber() <= 17768 && run.getTargetType() == Targ ){
		    RunsInDataSet.push_back( run.getRunNumber() );
		}
		else if( Dataset == "Sp23" && TorPol == "Pos" && run.getTorusScale() > 0 && SolPol == "Neg" &&
		    run.getSolenoidScale() < 0 && run.getRunNumber() >= 17769 && run.getRunNumber() <= 17811 && run.getTargetType() == Targ ){
		    RunsInDataSet.push_back( run.getRunNumber() );
		}
	    }
	    return RunsInDataSet;
	}

	// Returns a vector of runs from a data set of a given polarization
	vector<int> getPolarizationEpoch(string Dataset, string Targ, string polSign) const{
	    // Returns the subset of runs consisting of this dataset, torus polarity, solenoid polarity
	    vector<int> RunsInDataSet;
	    for( auto run : Runs ){
		// Summer 2022 data
		if( Dataset == "Su22" && polSign == "Neg" && run.getTargetPolarization() < 0 && run.getRunNumber() >= 16134 && run.getRunNumber() <= 16772 && run.getTargetType() == Targ ){
		    RunsInDataSet.push_back( run.getRunNumber() );
		}
		else if( Dataset == "Su22" && polSign == "Pos" && run.getTargetPolarization() > 0 && run.getRunNumber() >= 16134 && run.getRunNumber() <= 16772 && run.getTargetType() == Targ ){
		    RunsInDataSet.push_back( run.getRunNumber() );
		}
		// Fall 2022 data (negative solenoid)
		else if( Dataset == "Fa22Neg" && polSign == "Neg" && run.getTargetPolarization() < 0 && run.getRunNumber() >= 16859 && run.getRunNumber() <= 17183 && run.getTargetType() == Targ ){
		    RunsInDataSet.push_back( run.getRunNumber() );
		}
		else if( Dataset == "Fa22Neg" && polSign == "Pos" && run.getTargetPolarization() > 0 && run.getRunNumber() >= 16859 && run.getRunNumber() <= 17183 && run.getTargetType() == Targ ){
		    RunsInDataSet.push_back( run.getRunNumber() );
		}
		// Fall 2022 data (positive solenoid)
		else if( Dataset == "Fa22Pos" && polSign == "Neg" && run.getTargetPolarization() < 0 && run.getRunNumber() >= 17188 && run.getRunNumber() <= 17408 && run.getTargetType() == Targ ){
		    RunsInDataSet.push_back( run.getRunNumber() );
		}
		else if( Dataset == "Fa22Pos" && polSign == "Pos" && run.getTargetPolarization() > 0 && run.getRunNumber() >= 17188 && run.getRunNumber() <= 17408 && run.getTargetType() == Targ ){
		    RunsInDataSet.push_back( run.getRunNumber() );
		}
		// Spring 2023 data (inbending)
		else if( Dataset == "Sp23" && polSign == "Neg" && run.getTargetPolarization() < 0 && run.getRunNumber() >= 17482 && run.getRunNumber() <= 17768 && run.getTargetType() == Targ ){
		    RunsInDataSet.push_back( run.getRunNumber() );
		}
		else if( Dataset == "Sp23" && polSign == "Pos" && run.getTargetPolarization() > 0 && run.getRunNumber() >= 17482 && run.getRunNumber() <= 17768 && run.getTargetType() == Targ ){
		    RunsInDataSet.push_back( run.getRunNumber() );
		}
		// Spring 2023 data (outbending)
		else if( Dataset == "Sp23" && polSign == "Neg" && run.getTargetPolarization() < 0 && run.getRunNumber() >= 17769 && run.getRunNumber() <= 17811 && run.getTargetType() == Targ ){
		    RunsInDataSet.push_back( run.getRunNumber() );
		}
		else if( Dataset == "Sp23" && polSign == "Pos" && run.getTargetPolarization() > 0 && run.getRunNumber() >= 17769 && run.getRunNumber() <= 17811 && run.getTargetType() == Targ ){
		    RunsInDataSet.push_back( run.getRunNumber() );
		}
	    }
	    return RunsInDataSet;
	}

	vector<int> getAllGoodRuns() const{
	    vector<int> AllGoodRuns;
	    for(size_t i=0; i<Runs.size(); i++){
		AllGoodRuns.push_back( Runs[i].getRunNumber() );
	    }
	    if( AllGoodRuns.size() == 0 ) cout << "Couldn't find any good runs\n";
	    return AllGoodRuns;
    
	}

        vector<int> getDataSetGoodRuns(string DataSet) const{
            vector<int> AllGoodRuns;
            for(size_t i=0; i<Runs.size(); i++){
                int run = Runs[i].getRunNumber();
                if( DataSet == "Su22" && run <= 16772 ) AllGoodRuns.push_back( run );
                else if( DataSet == "Fa22" && run > 16800 && run < 17450 ) AllGoodRuns.push_back( run );
                else if( DataSet == "Sp23" && run > 17450 ) AllGoodRuns.push_back( run );
            }
            if( AllGoodRuns.size() == 0 ) cout << "Couldn't find any good runs\n";
            return AllGoodRuns;

        }

	// HWPstatus = 1 for "in", and = 0 for "out"
	vector<int> getRunsForHWPandTPol( string DataSet, string Target, int HWPstatus, int TPolSign ){
	    vector<int> AllGoodRuns;
	    for(size_t i=0; i<Runs.size(); i++){

		int run = Runs[i].getRunNumber();
		double thisTPol = Runs[i].getTargetPolarization();
		string targ = Runs[i].getTargetType();

		int intTPol = 0;
		if( thisTPol < 0 ) intTPol = -1;
		else if( thisTPol > 0 ) intTPol = 1; // intTPol stays zero if it's unpolarized

		if( DataSet == "Su22" && Target == targ && run <= 16772 && HWPstatus == Runs[i].getHWPStatus() && TPolSign == intTPol ) AllGoodRuns.push_back( run );
		else if( DataSet == "Fa22Neg" && Target == targ && run > 16800 && run <= 17183 && HWPstatus == Runs[i].getHWPStatus() && TPolSign == intTPol ) AllGoodRuns.push_back( run );
		else if( DataSet == "Fa22Pos" && Target == targ && run >= 17188 && run < 17450 && HWPstatus == Runs[i].getHWPStatus() && TPolSign == intTPol ) AllGoodRuns.push_back( run );
		else if( DataSet == "Sp23Inb" && Target == targ && run > 17450 && run <= 17768 && HWPstatus == Runs[i].getHWPStatus() && TPolSign == intTPol ) AllGoodRuns.push_back( run );
		else if( DataSet == "Sp23Out" && Target == targ && run >= 17769 && HWPstatus == Runs[i].getHWPStatus() && TPolSign == intTPol ) AllGoodRuns.push_back( run );

	    }
	    if( AllGoodRuns.size() == 0 ) cout <<"ERROR: Couldn't find any "<< Target <<" runs for "<< DataSet <<", HWP = "<< HWPstatus <<", TPolSign = "<< TPolSign << endl;

	    return AllGoodRuns;
	}

	bool isFoilRun(int run) const{
	    for(size_t i=0; i<Runs.size(); i++){
		if( run == Runs[i].getRunNumber() && Runs[i].getTargetType() == "Foil" ) return true; 
	    }
	    return false;
	}
	bool isGoodRun(int run) const{
	    for(size_t i=0; i<Runs.size(); i++){
		if( run == Runs[i].getRunNumber() ) return true;
	    }
	    return false;
	}
	void Print(){
	    for(size_t i=0; i<Runs.size(); i++){
		Runs[i].Print();
	    }
	}
	// Destructor
	~RunPeriod() = default;
};

// Returns the beam polarization and error
vector<double> PT_Vals(int run, RunPeriod& Period){

    vector<double> pt_vals ={0,0};
    double tpolSign = Period.getTargetPolarization( run );

    if( run < 16800 && tpolSign > 0 ){
	pt_vals[0] = 0.71; pt_vals[1] = 0.03;
    }
    else if( run < 16800 && tpolSign < 0 ){
	pt_vals[0] = 0.66; pt_vals[1] = 0.03;
    }

    else if( run > 16800 && run <= 17183 && tpolSign > 0 ){
	pt_vals[0] = 0.72; pt_vals[1] = 0.02;
    }
    else if( run > 16800 && run <= 17183 && tpolSign < 0 ){
	pt_vals[0] = 0.69; pt_vals[1] = 0.02;
    }

    else if( run > 17183 && run <= 17408 && tpolSign < 0 ){
	pt_vals[0] = 0.70; pt_vals[1] = 0.03;
    }
    else if( run > 17183 && run <= 17408 && tpolSign > 0 ){
	pt_vals[0] = 0.67; pt_vals[1] = 0.03;
    }

    else if( run > 17450 && tpolSign > 0 ){
	pt_vals[0] = 0.67; pt_vals[1] = 0.03;
    }
    else if( run > 17450 && tpolSign < 0 ){
	pt_vals[0] = 0.62; pt_vals[1] = 0.03;
    }

    return pt_vals;

}

// This vector contains the bin boundaries for the Q2 bins
//const vector<double> Q2_Bin_Bounds = { /*1.0, 1.3094, 1.5632, 1.8661,*/ 2.2277,
//	2.6594, 3.1747, 3.7899, 4.5243, 5.4009, 6.4475, 7.6969, 9.1884, 10.9689};

vector<double> Q2_Bin_Bounds = { /*1.0, 1.3094, 1.5632, 1.8661,*/ 2.2277,
	2.6594, 3.1747, 3.7899, 4.5243, 5.4009, 6.4475, 7.6969, 9.1884, 10.9689};

const vector<double> X_Bin_Bounds = {0.075, 0.1, 0.125, 0.15, 0.175, 0.2, 0.225, 0.25, 0.275,
	0.3, 0.325, 0.35, 0.375, 0.4, 0.425, 0.45, 0.475,
	0.5, 0.525, 0.55, 0.575, 0.6, 0.625, 0.65, 0.675,
	0.7, 0.725, 0.75, 0.775, 0.8, 0.825, 0.85, 0.875 };

// This object holds information on the various scaling factors (like generating pseudo-CD2 counts) to be
// used in the calculation of dilution factors and packing fractions.


// This struct holds information about the counts of different target types
class Count{

    private:
	string TargetType = "N/A";
	double NP = 0.0; double NM = 0.0; double N0 = 0.0;
	double FC_P = 0.0; double FC_M = 0.0; double FC0 = 0.0;
	double NormNP = 0.0; double NormNM = 0.0; double NormN0 = 0.0;
	double RadLenCorr = 1; // Radiation length correction factor; used only for ET and F targets!!!
	double RadLenCorrErr = 0; // Error on the radiation length correction; used only for ET and F targets!!!
    public:
	// Constructor
	Count() = default;
	Count(string targtype) : TargetType(targtype){}
	// Mutators
	void SetTarget(string targtype){
		TargetType = targtype;
	}
	void AddCounts(double np, double nm, double n0){
		NP += np; NM += nm; N0 += n0;
	}
	void AddFCCharge(double fcp, double fcm, double fc0){
		FC_P += fcp; FC_M += fcm; FC0 += fc0;
	}
	void SetNormCounts(){
	    if( FC_P != 0.0 && FC_M != 0.0 ){
		NormNP = NP/FC_P; NormNM = NM/FC_M;
		if( FC0 != 0.0 ) NormN0 = N0/FC0;
	    }
	    // Don't print anything if it's an empty bin
	    else if(NP > 0.0 && NM > 0.0){
		cout << "ERROR: Invalid FC charge. Unable to set FC values.\n";
		cout << "FC_P = "<< FC_P <<", FC_M = "<< FC_M << endl;
	    }
	}
	// Accessors
	double getNP() const{return NP;}
	double getNM() const{return NM;}
	double getNt() const{return NP+NM;}
	double getN0() const{return N0;}
	double getFC_P() const{return FC_P;}
	double getFC_M() const{return FC_M;}
	double getFCt() const{return FC_P+FC_M;}
	double getFC0() const{return FC0;}
	double getRadLenCorr() const{return RadLenCorr;}
	double getRadLenCorrErr() const{return RadLenCorrErr;}
	double getNormNP() const{
	    if( FC_P > 0.0 ) return NP / FC_P;
	    else return 0.0;
	}
	double getNormNM() const{
	    if( FC_M > 0.0 ) return NM / FC_M;
	    else return 0.0;
	}
	double getNormN0() const{
	    if( FC0 > 0.0 ) return N0 / FC0;
	    else return 0.0;
	}
	//double getNormNt() const{return NormNP+NormNM;}
	double getNormNt() const{
	    if( FC_P > 0.0 && FC_M > 0.0 )
		return (NP+NM)/(FC_P+FC_M);
	    else
		return 0.0;
	}

	void Print() const{
		cout << "Target = "<< TargetType;
		cout << ",	NP = "<< NP <<",	NM = "<< NM;
	       	cout << ",	FC_P = "<< FC_P <<",	FC_M = "<< FC_M;
		cout << ",	NormNP = "<<NormNP<<",	NormNM = "<<NormNM << endl;
	}
	// Destructor
	~Count() = default;
};


// This struct holds the information for each of the x, Q2 bins
class Bin{

    private:
	// Minimum, maximum, and middle value for the bin (mid = (max+min)/2)
	double Q2_Min = 0.0; double Q2_Max = 0.0; double Q2_Mid = 0.0;
	double X_Min  = 0.0; double X_Max  = 0.0; double X_Mid  = 0.0;
	// Holds the statistics-weighted average values for the x, Q2 bins separately for each target type.
	// These are set by the member function "SetAvgXQ2" in the "DataSet" class, which takes a vector of runs 
	// and a string for the target type to set the values.
	double Q2_Avg_NH3 = 0.0; double Q2_Avg_ND3 = 0.0; double Q2_Avg_C = 0.0; double Q2_Avg_F = 0.0;
	double Q2_Avg_CH2 = 0.0; double Q2_Avg_CD2 = 0.0; double Q2_Avg_ET = 0.0;
	double X_Avg_NH3 = 0.0; double X_Avg_ND3 = 0.0; double X_Avg_C = 0.0; double X_Avg_F = 0.0;
	double X_Avg_CH2 = 0.0; double X_Avg_CD2 = 0.0; double X_Avg_ET = 0.0;
	// These are the same as above, but for multiple target types taken together.
	// Primarily, these will be used for determining the average bin location ofthey're calculated separately from them
	double Q2_Avg_All = 0.0; double X_Avg_All = 0.0;

	// Holds the correction ratios used for generating pseudo data.
	double ScaleCD2 = 0.0; double ScaleCD2err = 0.0; // Used to turn CH2 data into "CD2" data
	double ScaleSol = 0.0; double ScaleSolerr = 0.0; // Used to convert neg. solenoid "ET" and "F" data into pos. solenoid data

	// The following ScaleET and ScaleF have been deprecated for now, but I left them in.
	double ScaleET  = 0.0; double ScaleETerr  = 0.0;
	double ScaleF   = 0.0; double ScaleFerr   = 0.0;

	// Holds the calculated dilution factor and its associated error for this bin
	double DF_NH3 = 0.0; double ErrDF_NH3 = 0.0; // This calculates the DF using unique PF_bath for each bin
	double SysErrDF_NH3 = 0.0; // Systematic error
	double DF_FixedPF_NH3 = 0.0;   // This DF value is calculated using a fixed value for the PF_bath across all bins
	double ErrDF_FixedPF_NH3 = 0.0;
	double DF_ND3 = 0.0; double ErrDF_ND3 = 0.0;
	double SysErrDF_ND3 = 0.0;
	double DF_FixedPF_ND3 = 0.0;   // This DF value is calculated using a fixed value for the PF_bath across all bins
	double ErrDF_FixedPF_ND3 = 0.0;
	// Holds the calculated packing fraction and the associated error for this bin
	double PF_bath_NH3 = 0.0; double ErrPF_bath_NH3 = 0.0;
	double PF_cell_NH3 = 0.0; double ErrPF_cell_NH3 = 0.0;
	double PF_bath_ND3 = 0.0; double ErrPF_bath_ND3 = 0.0;
	double PF_cell_ND3 = 0.0; double ErrPF_cell_ND3 = 0.0;
	// Holds information for various pieces of the DF calculation; for debugging
	double Numerator_NH3_DF_Thru_PF = 0.0; double Numerator_NH3_DF_Thru_PF_Error = 0.0;
	double NA_Counts_NH3_DF_Thru_PF = 0.0; double NA_Counts_NH3_DF_Thru_PF_Error = 0.0;
	double Numerator_ND3_DF_Thru_PF = 0.0; double Numerator_ND3_DF_Thru_PF_Error = 0.0;
	double NA_Counts_ND3_DF_Thru_PF = 0.0; double NA_Counts_ND3_DF_Thru_PF_Error = 0.0;
	// These asymmetries are for sign-checking in the input loop, NOT FOR FINAL BINNING!!!
	double AllRawNoFCInput = 0.0; double AllRawNoFCInputErr = 0.0;
	double AllRawInput = 0.0; double AllRawInputErr = 0.0;
	// Holds the non-FC corrected raw double-spin asymmetry for NH3 and ND3
	double AllRawNoFCNH3 = 0.0; double AllRawNoFCND3 = 0.0;
	double AllRawNoFCNH3Err = 0.0; double AllRawNoFCND3Err = 0.0;
	// Holds the FC asymmetry for NH3 and ND3
	double AFC_NH3 = 0.0; double AFC_ND3 = 0.0;
	double AFC_NH3Err = 0.0; double AFC_ND3Err = 0.0;

	// Holds the raw double-spin asymmetry for ND3
	double AllRawNH3 = 0.0; double AllRawND3 = 0.0;
	double AllRawNH3Err = 0.0; double AllRawND3Err = 0.0;
	// Holds the physical double-spin asymmetry
	double AllPhysNH3 = 0.0; double AllPhysND3 = 0.0;
	double AllPhysNH3Err = 0.0; double AllPhysND3Err = 0.0;
	// Holds the theoretical value of the double-spin asymmetry for NH3 in this bin
	double AllTheoryNH3 = 0.0; double AllTheoryND3 = 0.0;
        // Holds various other kinematic factors (A2 is from Sebastian's theoretical values)
        double DepolFactorNH3= 0.0; double EtaNH3 = 0.0; double A2NH3 = 0.0;
        double DepolFactorND3= 0.0; double EtaND3 = 0.0; double A2ND3 = 0.0;
	// Holds the value of A1 from data and theory
	double A1_NH3 = 0.0; double ErrA1_NH3 = 0.0;
	double A1_ND3 = 0.0; double ErrA1_ND3 = 0.0;
	double A1_Theory_NH3 = 0.0; double A1_Theory_ND3 = 0.0;

	// Hold the count types for every kind of target
	Count NH3_Counts = Count("NH3"); // Ammonia
	Count ND3_Counts = Count("ND3"); // Deuterated ammonia
	Count C_Counts   = Count("C");   // Carbon foils
	Count CH2_Counts = Count("CH2"); // CH2
	Count CD2_Counts = Count("CD2"); // CD2
	Count ET_Counts  = Count("ET");  // Empty target, with LHe
	Count F_Counts   = Count("F");   // Empty target, no LHe (Aluminum foils only)
	Count Input_Loop = Count("Input");//USED ONLY IN "Get_DF_By_Sector.C" FOR WRITING TO INPUT TEXT FILES
    public:
	// Constructors
	Bin() = default;

	// Mutators
	void SetQ2Bins(double qmin, double qmax){
		Q2_Min = qmin; Q2_Max = qmax; Q2_Mid = (qmin+qmax)/2.0;
	}
	void SetXBins(double xmin, double xmax){
		X_Min = xmin; X_Max = xmax; X_Mid = (xmin+xmax)/2.0;
	}
	void SetDF_NH3(double df, double dferr){
		DF_NH3 = df; ErrDF_NH3 = dferr;
	}
	void SetSysErrDF(double syserr, string target){
		if( target == "NH3" ) SysErrDF_NH3 = syserr;
		else if( target == "ND3" ) SysErrDF_ND3 = syserr;
		else cout << "ERROR: Unable to set systematic error for target "<<target<<". Check inputs.\n";
	}
	void SetDF_ND3(double df, double dferr){
		DF_ND3 = df; ErrDF_ND3 = dferr;
	}
	void SetPF(double pf, double pferr, string target, string cellType){
		if( target == "NH3" && cellType == "Bath" ){ PF_bath_NH3 = pf; ErrPF_bath_NH3 = pferr; }
		else if( target == "NH3" && cellType == "Cell" ){ PF_cell_NH3 = pf; ErrPF_cell_NH3 = pferr; }
		else if( target == "ND3" && cellType == "Bath" ){ PF_bath_ND3 = pf; ErrPF_bath_ND3 = pferr; }
		else if( target == "ND3" && cellType == "Cell" ){ PF_cell_ND3 = pf; ErrPF_cell_ND3 = pferr; }
		else cout <<"ERROR: Invalid target or cell type: 'NH3', 'ND3' targets; 'Cell', 'Bath' cell types.\n"; 
	}
	void SetDF_FixedPF_NH3(double df, double dferr){
		DF_FixedPF_NH3 = df; ErrDF_FixedPF_NH3 = dferr;
	}
	void SetDF_FixedPF_ND3(double df, double dferr){
		DF_FixedPF_ND3 = df; ErrDF_FixedPF_ND3 = dferr;
	}
	void SetAllRaw(double allraw, string target){
		if( target == "NH3" ) AllRawNH3 = allraw;
		else if( target == "ND3" ) AllRawND3 = allraw;
	}
	void SetAllRawErr(double allrawerr, string target){
		if( target == "NH3" ) AllRawNH3Err = allrawerr;
		else if( target == "ND3" ) AllRawND3Err = allrawerr;
	}
	void SetAllRawNoFC(double allraw, string target){
		if( target == "NH3" ) AllRawNoFCNH3 = allraw;
		else if( target == "ND3" ) AllRawNoFCND3 = allraw;
	}
	void SetAllRawNoFCErr(double allrawerr, string target){
		if( target == "NH3" ) AllRawNoFCNH3Err = allrawerr;
		else if( target == "ND3" ) AllRawNoFCND3Err = allrawerr;
	}
	void SetAllTheory(double allth, string target){
		if( target == "NH3" ) AllTheoryNH3 = allth;
		else if( target == "ND3" ) AllTheoryND3 = allth;
	}
	void SetAllPhys(double allphys, double errallphys, string target){
		if( target == "NH3" ){ AllPhysNH3 = allphys; AllPhysNH3Err = errallphys;}
		else if( target == "ND3" ){ AllPhysND3 = allphys; AllPhysND3Err = errallphys;}
	}
	void SetKinematicFactors(double depol, double eta, double a1, double a2, string target){
		if( target == "NH3" ){ 
		    DepolFactorNH3 = depol; EtaNH3 = eta; A1_Theory_NH3 = a1; A2NH3 = a2;
		}
		else if( target == "ND3" ){ 
		    DepolFactorND3 = depol; EtaND3 = eta; A1_Theory_ND3 = a1; A2ND3 = a2;
		}
	}
	void SetA1(double a1, double err_a1, string target ){
		if( target == "NH3" ){
		    A1_NH3 = a1; ErrA1_NH3 = err_a1;
		}
		else if( target == "ND3" ){
		    A1_ND3 = a1; ErrA1_ND3 = err_a1;
		}
	}
	void SetRatio( string scaleType, double scale, double err ){
		if( scaleType == "CD2" ){ ScaleCD2 = scale; ScaleCD2err = err; }
		else if( scaleType == "Sol" ){ ScaleSol = scale; ScaleSolerr = err; }
		else if( scaleType == "ET"  ){ ScaleET = scale; ScaleETerr = err; }
		else if( scaleType == "F"   ){ ScaleF = scale; ScaleFerr = err; }
		else cout << "WARNING: Unable to set scaling factors for scaleType = "<< scaleType <<". Check inputs.\n";
	}
	void SetAFC( double afc, string target ){
		if( target == "NH3" ) AFC_NH3 = afc;
		else if( target == "ND3" ) AFC_ND3 = afc;
		else cout << "WARNING: Unable to set FC asymmetry for "<< target <<". Check inputs.\n";
	}
	// This function is used to set the average values of X and Q2 for a given target type
	void SetAvgXQ2(string target, double xavg, double q2avg){
		if( target == "NH3" ){ Q2_Avg_NH3 = q2avg; X_Avg_NH3 = xavg; }
		else if( target == "ND3" ){ Q2_Avg_ND3 = q2avg; X_Avg_ND3 = xavg; }
		else if( target == "CH2" ){ Q2_Avg_CH2 = q2avg; X_Avg_CH2 = xavg; }
		else if( target == "CD2" ){ Q2_Avg_CD2 = q2avg; X_Avg_CD2 = xavg; }
		else if( target == "C"   ){ Q2_Avg_C = q2avg; X_Avg_C = xavg; }
		else if( target == "ET"  ){ Q2_Avg_ET = q2avg; X_Avg_ET = xavg; }
		else if( target == "F"   ){ Q2_Avg_F = q2avg; X_Avg_F = xavg; }
		else if( target == "All" ){ Q2_Avg_All = q2avg; X_Avg_All = xavg; }
		else cout <<"ERROR: Couldn't set target for "<<target<<" with average values Q2 = "<< q2avg <<", X = "<< xavg << endl;
	}
	void CalculateA1(double PBPT, double pbpt_err, string Variation = "ByBin"){
		//cout << "For bin "<< Q2_Mid <<" "<< X_Mid <<": DepolNH3 = "<< DepolFactorNH3 <<", EtaNH3 = "<< EtaNH3 <<", A2_NH3 = "<< A2NH3;
	        //cout << ", DF_FixedPF_NH3 = "<< DF_FixedPF_NH3 <<", PbPt = "<< PBPT << endl;
		double pbpt = abs(PBPT);
		if( Variation == "thruPF" && DepolFactorNH3 > 0 && EtaNH3 > 0 && A2NH3 > 0 && pbpt > 0 && DF_FixedPF_NH3 > 0 ){
		    AllPhysNH3 = (AllRawNH3 / (pbpt*DF_FixedPF_NH3*DepolFactorNH3)) - (EtaNH3*A2NH3);
		    AllPhysNH3Err = (AllRawNH3 / (pbpt*DF_FixedPF_NH3*DepolFactorNH3))*sqrt( pbpt_err*pbpt_err/(pbpt*pbpt) + ErrDF_FixedPF_NH3*ErrDF_FixedPF_NH3/(DF_FixedPF_NH3*DF_FixedPF_NH3) );
		}
		else if( Variation == "ByBin" && DepolFactorNH3 > 0 && EtaNH3 > 0 && A2NH3 > 0 && pbpt > 0 && DF_FixedPF_NH3 > 0 ){
		    AllPhysNH3 = (AllRawNH3 / (pbpt*DF_NH3*DepolFactorNH3)) - (EtaNH3*A2NH3);
		    AllPhysNH3Err = (AllRawNH3 / (pbpt*DF_NH3*DepolFactorNH3))*sqrt( pbpt_err*pbpt_err/(pbpt*pbpt) + ErrDF_NH3*ErrDF_NH3/(DF_NH3*DF_NH3) );
		}
		else{
		    cout << "ERROR: Invalid values for bin "<< Q2_Mid <<" "<< X_Mid << endl;
		}
	}
	// NOTE: The correction factors use the bin evaluated at the ET or Foil average values of x, Q2;
	// They don't use the "All" average value across all target types!!!!!!!!!!!!!!!!!!
	void CalculateDF( string Period, bool useScaling, bool usePseudoData ){
		double nA = NH3_Counts.getNt(); double fA = NH3_Counts.getFCt();
		double nD = ND3_Counts.getNt(); double fD = ND3_Counts.getFCt();
		double nCH= CH2_Counts.getNt(); double fCH= CH2_Counts.getFCt();
		double nCD= CD2_Counts.getNt(); double fCD= CD2_Counts.getFCt();
		double nC = C_Counts.getNt();   double fC = C_Counts.getFCt();
		double nET = ET_Counts.getNt(); double fET = ET_Counts.getFCt();
		double nF = F_Counts.getNt();   double fF = F_Counts.getFCt();

		// Calculate NH3 dilution factor for this bin
		if( fA>0 && fCH>0 && fC>0 && fET>0 && fF>0 && nA>0 && nCH>0 && nC>0 && nET>0 /*&& nF>0*/ ){
			// Propagate the errors of the pseudo data, if applicable
			double ETsf = 1.0; // Factor used to make the radiation length correction to the empty target
			double scale_nET = nET; double scale_nF = nF; // Raw counts used for scaling the ET and F counts
			double errET = sqrt(nET)/fET; double errF = sqrt(nF)/fF; // Default values that may scale if there are scaling factors
			if( useScaling ){ ETsf = ET_Scale_Factor( X_Avg_ET, Q2_Avg_ET ); /*cout << "Set the ET scaling factor.\n";*/ }

			if( usePseudoData && Period == "Fa22Pos" ){
			    scale_nET = nET*ScaleSol;
			    scale_nF  = nF*ScaleSol;
			    // Overwrite the ET and F errors to propagate errors correctly...
			    errET = (1.0/fET)*sqrt( ScaleSol*scale_nET + (pow(scale_nET/ScaleSol,2))*pow(ScaleSolerr,2) );
			    errF  = (1.0/fF )*sqrt( ScaleSol*scale_nF  + (pow(scale_nF /ScaleSol,2))*pow(ScaleSolerr,2) );
			    //cout << "Set the solenoid scaling factor for NH3 data.\n";
			}

			double dfnh3 = NH3DF( nA/fA, nCH/fCH, nC/fC, ETsf*scale_nET/fET, scale_nF/fF );
			double errdfnh3 = NH3DFError( nA/fA, nCH/fCH, nC/fC, ETsf*scale_nET/fET, scale_nF/fF, sqrt(nA)/fA, sqrt(nCH)/fCH, sqrt(nC)/fC, ETsf*errET, errF );

			double pfnh3 = NH3PF( nA/fA, nCH/fCH, nC/fC, ETsf*scale_nET/fET, scale_nF/fF );
			double errpfnh3 = NH3PFError( nA/fA, nCH/fCH, nC/fC, ETsf*scale_nET/fET, scale_nF/fF, sqrt(nA)/fA, sqrt(nCH)/fCH, sqrt(nC)/fC, ETsf*errET, errF );

			//cout <<"DF = "<< dfnh3 <<" +- "<< errdfnh3 <<", PF = "<< pfnh3 <<" +- "<< errpfnh3 << endl;
			//cout <<  nA/fA <<"  "<< nCH/fCH <<"  "<< nC/fC <<"  "<< ETsf*scale_nET/fET <<"  "<< scale_nF/fF << endl; 
			//cout <<  nA <<"  "<< nCH <<"  "<< nC <<"  "<< ETsf*scale_nET <<"  "<< scale_nF << endl; 

			if( dfnh3 > 0.0 && errdfnh3 > 0.0 && pfnh3 > 0.0 && errpfnh3 > 0.0 ){
			    DF_NH3 = dfnh3; ErrDF_NH3 = errdfnh3;
			    PF_bath_NH3 = pfnh3; ErrPF_bath_NH3 = errpfnh3;
			    PF_cell_NH3 = (LHe / Lcell)*pfnh3; ErrPF_cell_NH3 = (LHe / Lcell)*errpfnh3;
/*
			    cout << endl;
			    cout <<"Q2 = "<< Q2_Mid <<", X = "<< X_Mid <<"\n";
			    cout <<" --> scaled_nET = "<< ETsf*scale_nET <<" +- "<< errET << endl;
			    cout <<" --> scaled_nF  = "<< scale_nF  <<" +- "<< errF  << endl;
			    cout << endl;
*/
			}
			else{
			    cout <<"WARNING in CalculateDF() for bin X = "<< X_Mid <<", Q2 = "<< Q2_Mid << endl;
			    cout <<"DF = "<< dfnh3 <<" +- "<< errdfnh3 <<", PF = "<< pfnh3 <<" +- "<< errpfnh3 << endl;
			    cout <<  nA/fA <<"  "<< nCH/fCH <<"  "<< nC/fC <<"  "<< ETsf*scale_nET/fET <<"  "<< scale_nF/fF << endl;
			    cout <<"Zero or negative value for NH3 DF calculated. This value was omitted...\n";
			}
		}
		else if(nA > 0.0){
			cout <<"WARNING in CalculateDF() for bin X = "<< X_Mid <<", Q2 = "<< Q2_Mid << endl;
			cout <<"Some FC charges or counts are zero for NH3 calculation; check inputs.\n";
		}
		// Calculate ND3 dilution factor for this bin
		// Using CH2 in place of CD2 for the time being
		if( fD>0 && fCD>0 && fC>0 && fET>0 && fF>0 && nD>0 && nCD>0 && nC>0 && nET>0 && nF>0 ){
			double ETsf = 1.0;
			double scale_nCD = nCD; double scale_nET = nET; double scale_nF = nF;
			double errET = sqrt(nET)/fET; double errF = sqrt(nF)/fF; // Default values that may scale if there are scaling factors
			double errCD = sqrt(nCD)/fCD;
			if( useScaling ){ ETsf = ET_Scale_Factor( X_Avg_ET, Q2_Avg_ET ); /*cout << "Set the ET scaling factor.\n";*/ }
			if( usePseudoData ){
			    if( Period == "Su22" || Period == "Fa22Neg" || Period == "Fa22Pos" ){
				scale_nCD = nCD*ScaleCD2;
				errCD = (1.0/fCD )*sqrt( ScaleCD2 *scale_nCD  + (pow(scale_nCD /ScaleCD2 ,2))*pow(ScaleCD2err ,2) );
				//cout << "Set CD2 count ratio factors.\n";
			    }
			    // This is a separate check from above
			    if( Period == "Fa22Pos" ){
				scale_nET = nET*ScaleSol;
				scale_nF  = nF*ScaleSol;
				// Overwrite the ET and F errors to propagate errors correctly...
				errET = (1.0/fET)*sqrt( ScaleSol*scale_nET + (pow(scale_nET/ScaleSol,2))*pow(ScaleSolerr,2) );
				errF  = (1.0/fF )*sqrt( ScaleSol*scale_nF  + (pow(scale_nF /ScaleSol,2))*pow(ScaleSolerr,2) );
				//cout << "Set the solenoid scaling factor for ND3 data.\n";
			    }
			}

			double dfnd3 = ND3DF( nD/fD, scale_nCD/fCD, nC/fC, ETsf*scale_nET/fET, scale_nF/fF );
			double errdfnd3 = ND3DFError( nD/fD, scale_nCD/fCD, nC/fC, ETsf*scale_nET/fET, scale_nF/fF, sqrt(nD)/fD, errCD, sqrt(nC)/fC, ETsf*errET, errF );

			double pfnd3 = ND3PF( nD/fD, scale_nCD/fCD, nC/fC, ETsf*scale_nET/fET, scale_nF/fF );
			double errpfnd3 = ND3PFError( nD/fD, scale_nCD/fCD, nC/fC, ETsf*scale_nET/fET, scale_nF/fF, sqrt(nD)/fD, errCD, sqrt(nC)/fC, ETsf*errET, errF );

			if( dfnd3 > 0.0 && errdfnd3 > 0.0 && pfnd3 > 0.0 && errpfnd3 > 0.0 ){
			    DF_ND3 = dfnd3; ErrDF_ND3 = errdfnd3;
			    PF_bath_ND3 = pfnd3; ErrPF_bath_ND3 = errpfnd3;
			    PF_cell_ND3 = (LHe / Lcell)*pfnd3; ErrPF_cell_ND3 = (LHe / Lcell)*errpfnd3;
//			    cout << endl;
/*
			    cout <<"Q2 = "<< Q2_Mid <<", X = "<< X_Mid <<"\n";
			    cout <<" --> scaled_nET = "<< ETsf*scale_nET <<" +- "<< errET << endl;
			    cout <<" --> scaled_nF  = "<< scale_nF  <<" +- "<< errF  << endl;
			    cout <<" --> scaled_nCD = "<< scale_nCD <<" +- "<< errCD << endl;
			    cout << endl;
*/
			}
			else{
			    cout <<"WARNING in CalculateDF() for bin X = "<< X_Mid <<", Q2 = "<< Q2_Mid << endl;
			    cout <<"Zero or negative value for ND3 DF calculated. This value was omitted...\n";
			}
		}
		else if(nD > 0.0){
			cout <<"WARNING in CalculateDF() for bin X = "<< X_Mid <<", Q2 = "<< Q2_Mid << endl;
			cout <<"Some FC charges or counts are zero for ND3 calculation; check inputs.\n";
		}
		
	}


	// NOTE: The correction factors use the bin evaluated at the ET or Foil average values of x, Q2;
	// They don't use the "All" average value across all target types!!!!!!!!!!!!!!!!!!
	void CalculatePF( string Period, bool useScaling, bool usePseudoData ){
		NH3_Counts.SetNormCounts(); double nA = NH3_Counts.getNt(); double fA = NH3_Counts.getFCt();
		ND3_Counts.SetNormCounts(); double nD = ND3_Counts.getNt(); double fD = ND3_Counts.getFCt();
		CH2_Counts.SetNormCounts(); double nCH= CH2_Counts.getNt(); double fCH= CH2_Counts.getFCt();
		CD2_Counts.SetNormCounts(); double nCD= CD2_Counts.getNt(); double fCD= CD2_Counts.getFCt();
		C_Counts.SetNormCounts();   double nC = C_Counts.getNt();   double fC = C_Counts.getFCt();
		ET_Counts.SetNormCounts();  double nET = ET_Counts.getNt(); double fET = ET_Counts.getFCt();
		double ETsf = 1.0;
		if( useScaling ) ETsf = ET_Scale_Factor( X_Avg_ET, Q2_Avg_ET ); 

		F_Counts.SetNormCounts();   double nF = F_Counts.getNt();   double fF = F_Counts.getFCt();
		// Calculate NH3 packing fraction for this bin
		if( fA>0 && fCH>0 && fC>0 && fET>0 && fF>0 && nA>0 && nCH>0 && nC>0 && ETsf*nET>0 && nF>0 ){

			double pfnh3 = NH3PF( nA/fA, nCH/fCH, nC/fC, ETsf*nET/fET, nF/fF );
			double errpfnh3 = NH3PFError( nA/fA, nCH/fCH, nC/fC, ETsf*nET/fET, nF/fF, sqrt(nA)/fA, sqrt(nCH)/fCH, sqrt(nC)/fC, ETsf*sqrt(nET)/fET, sqrt(nF)/fF );
			if( pfnh3 > 0.0 && errpfnh3 > 0.0 ){
			    PF_bath_NH3 = pfnh3; ErrPF_bath_NH3 = errpfnh3;
			    PF_cell_NH3 = (LHe / 5.0)*pfnh3; ErrPF_cell_NH3 = (LHe / 5.0)*errpfnh3;
			}
			else{
			    cout <<"ERROR in CalculatePF() for bin X = "<< X_Mid <<", Q2 = "<< Q2_Mid << endl;
			    cout <<"Zero or negative value for NH3 PF calculated. This value was omitted...\n";
			}
		}
		else if(nA > 0.0){
			cout <<"***********************************************************************************\n";
			cout <<"ERROR in CalculatePF(): Some FC charges or counts are zero for NH3 calculation; check inputs.\n";
			cout <<"For Bin Q^2 = "<< Q2_Mid <<" (GeV^2), X = "<< X_Mid << endl;
			cout <<"  -> fA = " << fA << ", nA =	" << nA << endl;
			cout <<"  -> fCH = " << fCH << ", nCH =	" << nCH << endl;
			cout <<"  -> fC = " << fC << ", nC =	" << nC << endl;
			cout <<"  -> fET = " << fET << ", nET =	" << nET << endl;
			cout <<"  -> fF = " << fF << ", nF =	" << nF << endl;
			cout <<"For the Empty Target (ET):\n";
			cout <<"  -> Average Q^2 = "<< Q2_Avg_ET << endl;
			cout <<"  -> Average X = "<< X_Avg_ET << endl;
			cout <<"  -> ET Scale Factor = "<< ET_Scale_Factor( X_Avg_ET, Q2_Avg_ET ) << endl;
			cout <<"***********************************************************************************\n";
		}
		// Calculate ND3 packing fraction for this bin
		if( fD>0 && fCD>0 && fC>0 && fET>0 && fF>0 && nD>0 && nCD>0 && nC>0 && ETsf*nET>0 && nF>0 ){
			double pfnd3 = ND3PF( nD/fD, nCD/fCD, nC/fC, ETsf*nET/fET, nF/fF );
			double errpfnd3 = ND3PFError( nD/fD, nCD/fCD, nC/fC, ETsf*nET/fET, nF/fF, sqrt(nD)/fD, sqrt(nCD)/fCD, sqrt(nC)/fC, ETsf*sqrt(nET)/fET, sqrt(nF)/fF );
			if( pfnd3 > 0.0 && errpfnd3 > 0.0 ){
			    PF_bath_ND3 = pfnd3; ErrPF_bath_ND3 = errpfnd3;
			    PF_cell_ND3 = (LHe / 5.0)*pfnd3; ErrPF_cell_ND3 = (LHe / 5.0)*errpfnd3;
			}
			else{
			    cout <<"ERROR in CalculatePF() for bin X = "<< X_Mid <<", Q2 = "<< Q2_Mid << endl;
			    cout <<"Zero or negative value for ND3 PF calculated. This value was omitted...\n";
			}
		}
		else if(nD > 0.0){
			cout <<"***********************************************************************************\n";
			cout <<"ERROR in CalculatePF(): Some FC charges or counts are zero for ND3 calculation; check inputs.\n";
			cout <<"For Bin Q^2 = "<< Q2_Mid <<" (GeV^2), X = "<< X_Mid << endl;
			cout <<"  -> fD = " << fD << ", nD =	" << nD << endl;
			cout <<"  -> fCD = " << fCD << ", nCD =	" << nCD << endl;
			cout <<"  -> fC = " << fC << ", nC =	" << nC << endl;
			cout <<"  -> fET = " << fET << ", nET =	" << nET << endl;
			cout <<"  -> fF = " << fF << ", nF =	" << nF << endl;
			cout <<"For the Empty Target (ET):\n";
			cout <<"  -> Average Q^2 = "<< Q2_Avg_ET << endl;
			cout <<"  -> Average X = "<< X_Avg_ET << endl;
			cout <<"  -> ET Scale Factor = "<< ET_Scale_Factor( X_Avg_ET, Q2_Avg_ET ) << endl;
			cout <<"***********************************************************************************\n";
		}
	}
	// This function is called in a member function of the DataSet object, using the same name. Don't call this
	// object on its own; it should only be used by the member function in the DataSet class!!!
	void CalculateDFThruPF( double thisPF, double thisPFerr, string targetType, string Period, bool useScaling, bool usePseudoData ){
		NH3_Counts.SetNormCounts(); double nA = NH3_Counts.getNt(); double fA = NH3_Counts.getFCt();
		ND3_Counts.SetNormCounts(); double nD = ND3_Counts.getNt(); double fD = ND3_Counts.getFCt();
		CH2_Counts.SetNormCounts(); double nCH= CH2_Counts.getNt(); double fCH= CH2_Counts.getFCt();
		CD2_Counts.SetNormCounts(); double nCD= CD2_Counts.getNt(); double fCD= CD2_Counts.getFCt();
		C_Counts.SetNormCounts();   double nC = C_Counts.getNt();   double fC = C_Counts.getFCt();
		ET_Counts.SetNormCounts();  double nET = ET_Counts.getNt(); double fET = ET_Counts.getFCt();
		double ETsf = 1.0;
		if( useScaling ) ETsf = ET_Scale_Factor( X_Avg_ET, Q2_Avg_ET );
		 
		F_Counts.SetNormCounts();   double nF = F_Counts.getNt();   double fF = F_Counts.getFCt();        
		// Calculate NH3 packing fraction for this bin
		if( targetType == "NH3" ){
		    if( fA>0 && fCH>0 && fC>0 && fET>0 && fF>0 && thisPF>0 ){
			double dfnh3 = NH3DF_Thru_PF( thisPF, nA/fA, nCH/fCH, nC/fC, ETsf*nET/fET, nF/fF );
			double errdfnh3 = NH3DF_Thru_PF_Error(thisPF, nA/fA, nCH/fCH, nC/fC, ETsf*nET/fET, nF/fF, 0, sqrt(nA)/fA, sqrt(nCH)/fCH, sqrt(nC)/fC, ETsf*sqrt(nET)/fET, sqrt(nF)/fF );
			if( dfnh3 > 0.0 && errdfnh3 > 0.0 ){
			    DF_FixedPF_NH3 = dfnh3; ErrDF_FixedPF_NH3 = errdfnh3;
			}
			else{
			    cout <<"ERROR in CalculateDFThruPF() for bin X = "<< X_Mid <<", Q2 = "<< Q2_Mid << endl;
			    cout <<"Zero or negative value for NH3 DF calculated. This value was omitted...\n";
			}
			//Numerator_NH3_DF_Thru_PF = NH3DF_Thru_PF_Numerator( thisPF, nA/fA, nCH/fCH, nC/fC, nET/fET, nF/fF );
			//Numerator_NH3_DF_Thru_PF_Error = NH3DF_Thru_PF_Numerator_Error(thisPF, nA/fA, nCH/fCH, nC/fC, nET/fET, nF/fF, 0, sqrt(nA)/fA, sqrt(nCH)/fCH, sqrt(nC)/fC, sqrt(nET)/fET, sqrt(nF)/fF );
			//NA_Counts_NH3_DF_Thru_PF = NH3DF_Thru_PF_NA( thisPF, nA/fA, nCH/fCH, nC/fC, nET/fET, nF/fF );
			//NA_Counts_NH3_DF_Thru_PF_Error = NH3DF_Thru_PF_NA_Error(thisPF, nA/fA, nCH/fCH, nC/fC, nET/fET, nF/fF, 0, sqrt(nA)/fA, sqrt(nCH)/fCH, sqrt(nC)/fC, sqrt(nET)/fET, sqrt(nF)/fF );
		    }
		    else{
			cout <<"ERROR in CalculateDFThruPF() for bin X = "<< X_Mid <<", Q2 = "<< Q2_Mid << endl;
			cout <<"Some FC charges are zero for NH3 calculation; check inputs.\n";
		    }
		}
		// Calculate ND3 packing fraction for this bin
		else if( targetType == "ND3" ){
		    if( fD>0 && fCD>0 && fC>0 && fET>0 && fF>0 && thisPF>0 ){
			double dfnd3 = ND3DF_Thru_PF( thisPF, nD/fD, nCD/fCD, nC/fC, ETsf*nET/fET, nF/fF );
			double errdfnd3 = ND3DF_Thru_PF_Error(thisPF, nD/fD, nCD/fCD, nC/fC, ETsf*nET/fET, nF/fF, 0, sqrt(nD)/fD, sqrt(nCD)/fCD, sqrt(nC)/fC, ETsf*sqrt(nET)/fET, sqrt(nF)/fF );
			if( dfnd3 > 0.0 && errdfnd3 > 0.0 ){
			    DF_FixedPF_ND3 = dfnd3; ErrDF_FixedPF_ND3 = errdfnd3;
			}
			else{
			    cout <<"ERROR in CalculateDFThruPF() for bin X = "<< X_Mid <<", Q2 = "<< Q2_Mid << endl;
			    cout <<"Zero or negative value for ND3 DF calculated. This value was omitted...\n";
			}
			//Numerator_ND3_DF_Thru_PF = ND3DF_Thru_PF_Numerator( thisPF, nD/fD, nCD/fCD, nC/fC, nET/fET, nF/fF );
			//Numerator_ND3_DF_Thru_PF_Error = ND3DF_Thru_PF_Numerator_Error(thisPF, nD/fD, nCD/fCD, nC/fC, nET/fET, nF/fF, 0, sqrt(nD)/fD, sqrt(nCD)/fCD, sqrt(nC)/fC, sqrt(nET)/fET, sqrt(nF)/fF );
			//NA_Counts_ND3_DF_Thru_PF = ND3DF_Thru_PF_NA( thisPF, nD/fD, nCD/fCD, nC/fC, nET/fET, nF/fF );
			//NA_Counts_ND3_DF_Thru_PF_Error = ND3DF_Thru_PF_NA_Error(thisPF, nD/fD, nCD/fCD, nC/fC, nET/fET, nF/fF, 0, sqrt(nD)/fD, sqrt(nCD)/fCD, sqrt(nC)/fC, sqrt(nET)/fET, sqrt(nF)/fF );
			//cout << dfnd3 <<" +- "<< errdfnd3 << endl;
		    }
		    else{
			cout <<"ERROR in CalculateDFThruPF() for bin X = "<< X_Mid <<", Q2 = "<< Q2_Mid << endl;
			cout <<"Some FC charges are zero for ND3 calculation; check inputs.\n";
		    }
		}
		else cout <<"ERROR in CalculateDFThruPF(): Invalid target type. Can only calculate DF(PF) for NH3 and ND3 targets.\n";
	}
	// Calculates the raw double-spin asymmetry for this bin
	void CalculateAllRaw(){
	    NH3_Counts.SetNormCounts();
	    ND3_Counts.SetNormCounts();
	    Input_Loop.SetNormCounts();

	    // These are used only for "GetRawCounts.C" input file
	    double normNPInp = Input_Loop.getNormNP(); double normNMImp = Input_Loop.getNormNM();

	    // Calculate raw asymmetries for the NH3 data
	    if( NH3_Counts.getNormNP() > 0.0 && NH3_Counts.getNormNM() > 0.0 ){ 
		double rawNP = NH3_Counts.getNP(); double rawNM = NH3_Counts.getNM();
		double nm = NH3_Counts.getNormNP(); double np = NH3_Counts.getNormNM();
	        double fcm= NH3_Counts.getFC_P(); double fcp = NH3_Counts.getFC_M();

		AllRawNH3 = (nm - np) / (nm + np);
		//AllRawNH3Err = (2.0/pow((nm + np),2))*sqrt( (np*np*nm)/(fcm*fcm) + (nm*nm*np)/(fcp*fcp));
		AllRawNH3Err = (2.0/pow((nm + np),2))*sqrt( (np*np*nm)/(fcm) + (nm*nm*np)/(fcp));
		//AllRawNH3Err = 0.5 * sqrt( (nm + np)/(nm*np) );
		AllRawNoFCNH3 = (rawNP - rawNM)/(rawNP + rawNM);
		AllRawNoFCNH3Err = (2.0/pow((rawNP + rawNM),2))*sqrt( rawNP*rawNP*rawNM + rawNM*rawNM*rawNP );
		//AllRawNoFCNH3Err = 0.5 * sqrt( (nm + np)/(nm*np) );

		AFC_NH3 = (fcm - fcp) / (fcm + fcp);

	    }
	    else{
		//cout << "ERROR: FC-normalized counts might be zero; check inputs. Set AllRawNH3 = 0\n";
		AllRawNH3 = 0.0;
		AllRawNH3Err = 0.0;
	    }

	    // Calculate raw asymmetries for the ND3 data
	    if( ND3_Counts.getNormNP() > 0.0 && ND3_Counts.getNormNM() > 0.0 ){ 
		double rawNP = ND3_Counts.getNP(); double rawNM = ND3_Counts.getNM();
		double nm = ND3_Counts.getNormNP(); double np = ND3_Counts.getNormNM();
	        double fcm= ND3_Counts.getFC_P(); double fcp = ND3_Counts.getFC_M();

		AllRawND3 = (nm - np) / (nm + np);
		//AllRawND3Err = (2.0/pow((nm + np),2))*sqrt( (np*np*nm)/(fcm*fcm) + (nm*nm*np)/(fcp*fcp));
		AllRawND3Err = (2.0/pow((nm + np),2))*sqrt( (np*np*nm)/(fcm) + (nm*nm*np)/(fcp));
		//AllRawND3Err = 0.5 * sqrt( (nm + np)/(nm*np) );
		AllRawNoFCND3 = (rawNP - rawNM)/(rawNP + rawNM);
		AllRawNoFCND3Err = (2.0/pow((rawNP + rawNM),2))*sqrt( rawNP*rawNP*rawNM + rawNM*rawNM*rawNP );
		//AllRawNoFCND3Err = 0.5 * sqrt( (nm + np)/(nm*np) );

		AFC_ND3 = (fcm - fcp) / (fcm + fcp);
	    }
	    else{
		//cout << "ERROR: FC-normalized counts might be zero; check inputs. Set AllRawND3 = 0\n";
		AllRawND3 = 0.0;
		AllRawND3Err = 0.0;
	    }

	    // Calculate raw asymmetries for the Input data
	    if( Input_Loop.getNormNP() > 0.0 && Input_Loop.getNormNM() > 0.0 ){ 
		double rawNP = Input_Loop.getNP(); double rawNM = Input_Loop.getNM();
		double nm = Input_Loop.getNormNP(); double np = Input_Loop.getNormNM();
	        double fcm= Input_Loop.getFC_P(); double fcp = Input_Loop.getFC_M();

		AllRawInput = (nm - np) / (nm + np);
		//AllRawInputErr = (2.0/pow((nm + np),2))*sqrt( (np*np*nm)/(fcm*fcm) + (nm*nm*np)/(fcp*fcp));
		AllRawInputErr = (2.0/pow((nm + np),2))*sqrt( (np*np*nm)/(fcm) + (nm*nm*np)/(fcp));
		//AllRawInputErr = 0.5 * sqrt( (nm + np)/(nm*np) );
		AllRawNoFCInput = (rawNP - rawNM)/(rawNP + rawNM);
		AllRawNoFCInputErr = (2.0/pow((rawNP + rawNM),2))*sqrt( rawNP*rawNP*rawNM + rawNM*rawNM*rawNP );
		//AllRawNoFCInputErr = 0.5 * sqrt( (nm + np)/(nm*np) );
	    }
	    else{
		//cout << "ERROR: FC-normalized counts might be zero; check inputs. Set AllRawInput = 0\n";
		AllRawInput = 0.0;
		AllRawInputErr = 0.0;
	    }
	    /*
	    double normND3NP = ND3_Counts.getNormNP(); double normND3NM = ND3_Counts.getNormNM();
	    double nmND3 = ND3_Counts.getNP(); double npND3 = ND3_Counts.getNM();
	    if( normND3NP > 0.0 && normND3NM > 0.0 ){
		AllRawND3 = (normND3NP - normND3NM) / (normND3NP + normND3NM);
		AllRawND3Err = (2.0/pow((nmND3 + npND3),2))*sqrt(npND3*npND3*nmND3 + nmND3*nmND3*npND3);
		AllRawNoFCND3 = (nmND3 - npND3)/(nmND3 + npND3);
		AllRawNoFCND3Err = (2.0/pow((nmND3 + npND3),2))*sqrt(npND3*npND3*nmND3 + nmND3*nmND3*npND3);
	    }
	    else{
		//cout << "ERROR: FC-normalized counts might be zero; check inputs. Set AllRawND3 = 0\n";
		AllRawND3 = 0.0;
		AllRawND3Err = 0.0;
	    }
	    */
	}
	// Calculates the physical double-spin asymmetry for this bin
	void CalculateAllPhys( double PbPt, double PbPtErr, string Variation = "ByBin" ){
    	    //CalculateAllRaw(); // Make sure the raw asymmetry is set first
	    if( Variation == "thruPF" ){
	      if( AllRawNH3 != 0 && DF_FixedPF_NH3 != 0){
		AllPhysNH3 = AllRawNH3 / (PbPt*DF_FixedPF_NH3);
		AllPhysNH3Err = (1.0/(PbPt*DF_FixedPF_NH3))*sqrt( pow(AllRawNH3Err,2) + pow(AllRawNH3*ErrDF_FixedPF_NH3,2)/pow(DF_FixedPF_NH3,2) + pow(AllRawNH3*PbPtErr,2)/(PbPt*PbPt)  );
	      }
	      if( AllRawND3 != 0 ){
		AllPhysND3 = AllRawND3 / (PbPt*DF_FixedPF_ND3);
	      }
	    }
	    else if( Variation == "ByBin" ){
	      if( AllRawNH3 != 0 && DF_NH3 != 0){
		AllPhysNH3 = AllRawNH3 / (PbPt*DF_NH3);
		AllPhysNH3Err = (1.0/(PbPt*DF_NH3))*sqrt( pow(AllRawNH3Err,2) + pow(AllRawNH3*ErrDF_NH3,2)/pow(DF_NH3,2) + pow(AllRawNH3*PbPtErr,2)/(PbPt*PbPt)  );
	      }
	      if( AllRawND3 != 0 ){
		AllPhysND3 = AllRawND3 / (PbPt*DF_ND3);
	      }
	    }
	}

	// Automatically scales the ET counts
	void AddCounts(double nm, double np, double n0, string targtype, double xbin, double q2bin ){
		//double etScaleFactor = ET_Scale_Factor( xbin, q2bin );
		if(targtype == "NH3") NH3_Counts.AddCounts(nm, np, n0);
		else if(targtype == "ND3") ND3_Counts.AddCounts(nm, np, n0);
		else if(targtype == "CH2") CH2_Counts.AddCounts(nm, np, n0);
		else if(targtype == "CD2") CD2_Counts.AddCounts(nm, np, n0);
		else if(targtype == "C") C_Counts.AddCounts(nm, np, n0);
		//else if(targtype == "ET") ET_Counts.AddCounts(nm*etScaleFactor, np*etScaleFactor, n0*etScaleFactor);
		else if(targtype == "ET") ET_Counts.AddCounts(nm, np, n0);
		else if(targtype == "F") F_Counts.AddCounts(nm, np, n0);
		else if(targtype == "Input") Input_Loop.AddCounts(nm, np, n0);
		else{
		    cout <<"ERROR: Couldn't find targtype type. Added zero counts\n";
		}
	}
	void AddFCCharge(double fcm, double fcp, double fc0, string targtype){
		if(targtype == "NH3") NH3_Counts.AddFCCharge(fcm, fcp, fc0);
		else if(targtype == "ND3") ND3_Counts.AddFCCharge(fcm, fcp, fc0);
		else if(targtype == "CH2") CH2_Counts.AddFCCharge(fcm, fcp, fc0);
		else if(targtype == "CD2") CD2_Counts.AddFCCharge(fcm, fcp, fc0);
		else if(targtype == "C") C_Counts.AddFCCharge(fcm, fcp, fc0);
		else if(targtype == "ET") ET_Counts.AddFCCharge(fcm, fcp, fc0);
		else if(targtype == "F") F_Counts.AddFCCharge(fcm, fcp, fc0);
		else if(targtype == "Input") Input_Loop.AddFCCharge(fcm, fcp, fc0);
		else{
		    cout <<"ERROR: Couldn't find targtype type. Added zero counts\n";
		}
	}
	void NormalizeAllCounts(){
		NH3_Counts.SetNormCounts(); ND3_Counts.SetNormCounts();
		CH2_Counts.SetNormCounts(); CD2_Counts.SetNormCounts();
		C_Counts.SetNormCounts(); ET_Counts.SetNormCounts();
		F_Counts.SetNormCounts(); Input_Loop.SetNormCounts();
	}
	// Accessors
	double getBinQ2() const{ return Q2_Mid; }
	double getBinQ2Min() const{ return Q2_Min; }
	double getBinQ2Max() const{ return Q2_Max; }
	double getBinX() const{ return X_Mid; }
	double getBinXMin() const{ return X_Min; }
	double getBinXMax() const{ return X_Max; }
	double getDF_NH3() const{ return DF_NH3; }
	double getErrDF_NH3() const{ return ErrDF_NH3; }
	double getSysErrDF_NH3() const{ return SysErrDF_NH3; }
	double getDF_FixedPF_NH3() const{ return DF_FixedPF_NH3; }
	double getErrDF_FixedPF_NH3() const{ return ErrDF_FixedPF_NH3; }
	double getDF_FixedPF_ND3() const{ return DF_FixedPF_ND3; }
	double getErrDF_FixedPF_ND3() const{ return ErrDF_FixedPF_ND3; }
	double getDF_ND3() const{ return DF_ND3; }
	double getErrDF_ND3() const{ return ErrDF_ND3; }
	double getSysErrDF_ND3() const{ return SysErrDF_ND3; }
	double getPF_bath_NH3() const{ return PF_bath_NH3; }
	double getErrPF_bath_NH3() const{ return ErrPF_bath_NH3; }
	double getPF_cell_NH3() const{ return PF_cell_NH3; }
	double getErrPF_cell_NH3() const{ return ErrPF_cell_NH3; }
	double getPF_bath_ND3() const{ return PF_bath_ND3; }
	double getErrPF_bath_ND3() const{ return ErrPF_bath_ND3; }
	double getPF_cell_ND3() const{ return PF_cell_ND3; }
	double getErrPF_cell_ND3() const{ return ErrPF_cell_ND3; }
	double getDepolFactorNH3() const{ return DepolFactorNH3; }
	double getEtaNH3() const{ return EtaNH3; }
	double getA2NH3() const{ return A2NH3; }
	double getDepolFactorND3() const{ return DepolFactorND3; }
	double getEtaND3() const{ return EtaND3; }
	double getA2ND3() const{ return A2ND3; }
	double getA1NH3() const{ return A1_NH3;}
	double getErrA1NH3() const{ return ErrA1_NH3;}
	double getA1ND3() const{ return A1_ND3;}
	double getErrA1ND3() const{ return ErrA1_ND3;}
	double getAllPhysNH3() const{ return AllPhysNH3; }
	double getAllPhysNH3Err() const{ return AllPhysNH3Err; }

	// Returns the total error on DF, including the systematics (hasn't been implemented for DF thru PF though)
	double getTotalErrorDF( string target ) const{
		if( target == "NH3" ){
			if( SysErrDF_NH3 == 0 || ErrDF_NH3 == 0 ) cout << "WARNING: Not all DF_NH3 errors set for bin Q2="<< Q2_Avg_All<<", X="<<X_Avg_All<<endl;
			return sqrt( SysErrDF_NH3*SysErrDF_NH3 + ErrDF_NH3*ErrDF_NH3 );
		}
		else if( target == "ND3" ){
			if( SysErrDF_ND3 == 0 || ErrDF_ND3 == 0 ) cout << "WARNING: Not all DF_ND3 errors set for bin Q2="<< Q2_Avg_All<<", X="<<X_Avg_All<<endl;
			return sqrt( SysErrDF_ND3*SysErrDF_ND3 + ErrDF_ND3*ErrDF_ND3 );
		}
		else{
			cout <<"ERROR: Invalid target type "<< target <<", couldn't return valid error. Returning -1\n";
			return -1;
		}
	}

	double getRatio( string scaleType ) const{
		if( scaleType == "CD2" ){ return ScaleCD2; }
		else if( scaleType == "Sol" ){ return ScaleSol;}
		else if( scaleType == "ET"  ){ return ScaleET; }
		else if( scaleType == "F"   ){ return ScaleF; }
		else{ cout << "WARNING: Unable to return scaling factor for scaleType = "<< scaleType <<". Check inputs.\n"; return 0; };
	}
	double getRatioErr( string scaleType ) const{
		if( scaleType == "CD2" ){ return ScaleCD2err; }
		else if( scaleType == "Sol" ){ return ScaleSolerr;}
		else if( scaleType == "ET"  ){ return ScaleETerr; }
		else if( scaleType == "F"   ){ return ScaleFerr; }
		else{ cout << "WARNING: Unable to return scaling factor error for scaleType = "<< scaleType <<". Check inputs.\n"; return 0; };
	}

	double getAvgQ2(string target) const{
	    if( target == "NH3" ) return Q2_Avg_NH3;
	    else if( target == "ND3" ) return Q2_Avg_ND3;
	    else if( target == "CH2" ) return Q2_Avg_CH2;
	    else if( target == "CD2" ) return Q2_Avg_CD2;
	    else if( target == "C"   ) return Q2_Avg_C;
	    else if( target == "ET"  ) return Q2_Avg_ET;
	    else if( target == "F"   ) return Q2_Avg_F;
	    else if( target == "All" ) return Q2_Avg_All;
	    cout << "ERROR: Invalid target "<< target<<". No Q2_Avg value returned for Q2_Mid = "<<Q2_Mid;
	    cout << ", X_Mid = "<< X_Mid << endl;
	    return 0.0; // Error case
	}
	double getAvgX(string target) const{
	    if( target == "NH3" ) return X_Avg_NH3;
	    else if( target == "ND3" ) return X_Avg_ND3;
	    else if( target == "CH2" ) return X_Avg_CH2;
	    else if( target == "CD2" ) return X_Avg_CD2;
	    else if( target == "C"   ) return X_Avg_C;
	    else if( target == "ET"  ) return X_Avg_ET;
	    else if( target == "F"   ) return X_Avg_F;
	    else if( target == "All" ) return X_Avg_All;
	    cout << "ERROR: Invalid target "<< target<<". No Q2_Avg value returned for Q2_Mid = "<<Q2_Mid;
	    cout << ", X_Mid = "<< X_Mid << endl;
	    return 0.0; // Error case
	}

	double getNumeratorDFThruPF(string target) const{
	    if( target == "NH3" ) return Numerator_NH3_DF_Thru_PF;
	    else if( target == "ND3" ) return Numerator_ND3_DF_Thru_PF;
	    else return 0.0;
	}
	double getNumeratorDFThruPFError(string target) const{
	    if( target == "NH3" ) return Numerator_NH3_DF_Thru_PF_Error;
	    else if( target == "ND3" ) return Numerator_ND3_DF_Thru_PF_Error;
	    else return 0.0;
	}
	double getNACountsDFThruPF(string target) const{
	    if( target == "NH3" ) return NA_Counts_NH3_DF_Thru_PF;
	    else if( target == "ND3" ) return NA_Counts_ND3_DF_Thru_PF;
	    else return 0.0;
	}
	double getNACountsDFThruPFError(string target) const{
	    if( target == "NH3" ) return NA_Counts_NH3_DF_Thru_PF_Error;
	    else if( target == "ND3" ) return NA_Counts_ND3_DF_Thru_PF_Error;
	    else return 0.0;
	}
	double getAllTheory(string target) const{
	    if( target == "NH3" ) return AllTheoryNH3;
	    else if( target == "ND3" ) return AllTheoryND3;
	    else return 0.0;
	}
	double getAllRaw(string target) const{
	    if( target == "NH3" ) return AllRawNH3;
	    else if( target == "ND3" ) return AllRawND3;
	    else if( target == "Input" ) return AllRawInput;
	    else return 0.0;
	}
	double getAllRawErr(string target) const{
	    if( target == "NH3" ) return AllRawNH3Err;
	    else if( target == "ND3" ) return AllRawND3Err;
	    else if( target == "Input" ) return AllRawInputErr;
	    else return 0.0;
	}
	double getAFC(string target) const{
	    if( target == "NH3" ) return AFC_NH3;
	    else if( target == "ND3" ) return AFC_ND3;
	    else return 0.0;
	}
	double getAFC_Err(string target) const{
	    if( target == "NH3" ) return AFC_NH3Err;
	    else if( target == "ND3" ) return AFC_ND3Err;
	    else return 0.0;
	}

	double getAllRawNoFC(string target) const{
	    if( target == "NH3" ) return AllRawNoFCNH3;
	    else if( target == "ND3" ) return AllRawNoFCND3;
	    else if( target == "Input" ) return AllRawNoFCInput;
	    else return 0.0;
	}
	double getAllRawNoFCErr(string target) const{
	    if( target == "NH3" ) return AllRawNoFCNH3Err;
	    else if( target == "ND3" ) return AllRawNoFCND3Err;
	    else if( target == "Input" ) return AllRawNoFCInputErr;
	    else return 0.0;
	}

	double getA1Theory(string target) const{
	    if( target == "NH3" ) return A1_Theory_NH3; //DepolFactorNH3*(A1_Theory_NH3 + EtaNH3*A2NH3);
	    else if( target == "ND3" ) return A1_Theory_ND3; //DepolFactorND3*(A1_Theory_ND3 + EtaND3*A2ND3);
	    else return 0.0;
	}
	bool IsInBin(double qmid, double xmid) const{
	    if( qmid < Q2_Max && qmid >= Q2_Min && xmid < X_Max && xmid >= X_Min ) return true;
	    else return false;
	}
	double getNP(string target) const{
		if(target == "NH3") return NH3_Counts.getNP();
		else if(target == "ND3") return ND3_Counts.getNP();
		else if(target == "CH2") return CH2_Counts.getNP();
		else if(target == "CD2") return CD2_Counts.getNP();
		else if(target == "C") return C_Counts.getNP();
		else if(target == "ET") return ET_Counts.getNP();
		else if(target == "F") return F_Counts.getNP();
		else if(target == "Input") return Input_Loop.getNP();
		else{
		    cout <<"ERROR: Couldn't find target type. Return zero counts\n";
		    return 0.0;
		}
	}
	double getNM(string target) const{
		if(target == "NH3") return NH3_Counts.getNM();
		else if(target == "ND3") return ND3_Counts.getNM();
		else if(target == "CH2") return CH2_Counts.getNM();
		else if(target == "CD2") return CD2_Counts.getNM();
		else if(target == "C") return C_Counts.getNM();
		else if(target == "ET") return ET_Counts.getNM();
		else if(target == "F") return F_Counts.getNM();
		else if(target == "Input") return Input_Loop.getNM();
		else{
		    cout <<"ERROR: Couldn't find target type. Return zero counts\n";
		    return 0.0;
		}
	}
	double getN0(string target) const{
		if(target == "NH3") return NH3_Counts.getN0();
		else if(target == "ND3") return ND3_Counts.getN0();
		else if(target == "CH2") return CH2_Counts.getN0();
		else if(target == "CD2") return CD2_Counts.getN0();
		else if(target == "C") return C_Counts.getN0();
		else if(target == "ET") return ET_Counts.getN0();
		else if(target == "F") return F_Counts.getN0();
		else if(target == "Input") return Input_Loop.getN0();
		else{
		    cout <<"ERROR: Couldn't find target type. Return zero counts\n";
		    return 0.0;
		}
	}
	double getNt(string target) const{
		if(target == "NH3") return NH3_Counts.getNt();
		else if(target == "ND3") return ND3_Counts.getNt();
		else if(target == "CH2") return CH2_Counts.getNt();
		else if(target == "CD2") return CD2_Counts.getNt();
		else if(target == "C") return C_Counts.getNt();
		else if(target == "ET") return ET_Counts.getNt();
		else if(target == "F") return F_Counts.getNt();
		else if(target == "Input") return Input_Loop.getNt();
		else{
		    cout <<"ERROR: Couldn't find target type. Return zero counts\n";
		    return 0.0;
		}
	}
	double getNormNP(string target) const{
		if(target == "NH3") return NH3_Counts.getNormNP();
		else if(target == "ND3") return ND3_Counts.getNormNP();
		else if(target == "CH2") return CH2_Counts.getNormNP();
		else if(target == "CD2") return CD2_Counts.getNormNP();
		else if(target == "C") return C_Counts.getNormNP();
		else if(target == "ET") return ET_Counts.getNormNP();
		else if(target == "F") return F_Counts.getNormNP();
		else if(target == "Input") return Input_Loop.getNormNP();
		else{
		    cout <<"ERROR: Couldn't find target type. Return zero counts\n";
		    return 0.0;
		}
	}
	double getNormNM(string target) const{
		if(target == "NH3") return NH3_Counts.getNormNM();
		else if(target == "ND3") return ND3_Counts.getNormNM();
		else if(target == "CH2") return CH2_Counts.getNormNM();
		else if(target == "CD2") return CD2_Counts.getNormNM();
		else if(target == "C") return C_Counts.getNormNM();
		else if(target == "ET") return ET_Counts.getNormNM();
		else if(target == "F") return F_Counts.getNormNM();
		else if(target == "Input") return Input_Loop.getNormNM();
		else{
		    cout <<"ERROR: Couldn't find target type. Return zero counts\n";
		    return 0.0;
		}
	}
	double getNormN0(string target) const{
		if(target == "NH3") return NH3_Counts.getNormN0();
		else if(target == "ND3") return ND3_Counts.getNormN0();
		else if(target == "CH2") return CH2_Counts.getNormN0();
		else if(target == "CD2") return CD2_Counts.getNormN0();
		else if(target == "C") return C_Counts.getNormN0();
		else if(target == "ET") return ET_Counts.getNormN0();
		else if(target == "F") return F_Counts.getNormN0();
		else if(target == "Input") return Input_Loop.getNormN0();
		else{
		    cout <<"ERROR: Couldn't find target type. Return zero counts\n";
		    return 0.0;
		}
	}
	double getNormNt(string target) const{
		if(target == "NH3") return NH3_Counts.getNormNt();
		else if(target == "ND3") return ND3_Counts.getNormNt();
		else if(target == "CH2") return CH2_Counts.getNormNt();
		else if(target == "CD2") return CD2_Counts.getNormNt();
		else if(target == "C") return C_Counts.getNormNt();
		else if(target == "ET") return ET_Counts.getNormNt();
		else if(target == "F") return F_Counts.getNormNt();
		else if(target == "Input") return Input_Loop.getNormNt();
		else{
		    cout <<"ERROR: Couldn't find target type. Return zero counts\n";
		    return 0.0;
		}
	}
	double getAllNt() const{
		double nh3nt = NH3_Counts.getNt(); double nd3nt = ND3_Counts.getNt();
		double ch2nt = CH2_Counts.getNt(); double cd2nt = CD2_Counts.getNt();
		double carnt = C_Counts.getNt(); double empnt = ET_Counts.getNt();
		double foint = F_Counts.getNt();
		return nh3nt + nd3nt + ch2nt + cd2nt + carnt + empnt + foint;
	}
	double getFC_P(string target) const{
		if(target == "NH3") return NH3_Counts.getFC_P();
		else if(target == "ND3") return ND3_Counts.getFC_P();
		else if(target == "CH2") return CH2_Counts.getFC_P();
		else if(target == "CD2") return CD2_Counts.getFC_P();
		else if(target == "C") return C_Counts.getFC_P();
		else if(target == "ET") return ET_Counts.getFC_P();
		else if(target == "F") return F_Counts.getFC_P();
		else if(target == "Input") return Input_Loop.getFC_P();
		else{
		    cout <<"ERROR: Couldn't find target type. Return zero counts\n";
		    return 0.0;
		}
	}
	double getFC_M(string target) const{
		if(target == "NH3") return NH3_Counts.getFC_M();
		else if(target == "ND3") return ND3_Counts.getFC_M();
		else if(target == "CH2") return CH2_Counts.getFC_M();
		else if(target == "CD2") return CD2_Counts.getFC_M();
		else if(target == "C") return C_Counts.getFC_M();
		else if(target == "ET") return ET_Counts.getFC_M();
		else if(target == "F") return F_Counts.getFC_M();
		else if(target == "Input") return Input_Loop.getFC_M();
		else{
		    cout <<"ERROR: Couldn't find target type. Return zero counts\n";
		    return 0.0;
		}
	}
	double getFC0(string target) const{
		if(target == "NH3") return NH3_Counts.getFC0();
		else if(target == "ND3") return ND3_Counts.getFC0();
		else if(target == "CH2") return CH2_Counts.getFC0();
		else if(target == "CD2") return CD2_Counts.getFC0();
		else if(target == "C") return C_Counts.getFC0();
		else if(target == "ET") return ET_Counts.getFC0();
		else if(target == "F") return F_Counts.getFC0();
		else if(target == "Input") return Input_Loop.getFC0();
		else{
		    cout <<"ERROR: Couldn't find target type. Return zero counts\n";
		    return 0.0;
		}
	}
	double getFCt(string target) const{
		if(target == "NH3") return NH3_Counts.getFCt();
		else if(target == "ND3") return ND3_Counts.getFCt();
		else if(target == "CH2") return CH2_Counts.getFCt();
		else if(target == "CD2") return CD2_Counts.getFCt();
		else if(target == "C") return C_Counts.getFCt();
		else if(target == "ET") return ET_Counts.getFCt();
		else if(target == "F") return F_Counts.getFCt();
		else if(target == "Input") return Input_Loop.getFCt();
		else{
		    cout <<"ERROR: Couldn't find target type. Return zero counts\n";
		    return 0.0;
		}
	}

	double getAllFCt() const{
		double nh3fct = NH3_Counts.getFCt(); double nd3fct = ND3_Counts.getFCt();
		double ch2fct = CH2_Counts.getFCt(); double cd2fct = CD2_Counts.getFCt();
		double carfct = C_Counts.getFCt(); double empfct = ET_Counts.getFCt();
		double foifct = F_Counts.getFCt();
		return nh3fct + nd3fct + ch2fct + cd2fct + carfct + empfct + foifct;
	}
	// Used for making count ratio plots
	double getCountRatio(string target1, string target2, string Period){
		double targ1NormCounts = 0; double targ2NormCounts = 0;

		if(target1 == "NH3") targ1NormCounts = NH3_Counts.getNormNt();
		else if(target1 == "ND3") targ1NormCounts = ND3_Counts.getNormNt();
		else if(target1 == "CH2") targ1NormCounts = CH2_Counts.getNormNt();
		else if(target1 == "CD2") targ1NormCounts = CD2_Counts.getNormNt();
		else if(target1 == "C") targ1NormCounts = C_Counts.getNormNt();
		else if(target1 == "ET") targ1NormCounts = ET_Counts.getNormNt();
		else if(target1 == "F") targ1NormCounts = F_Counts.getNormNt();
		if(target2 == "NH3") targ2NormCounts = NH3_Counts.getNormNt();
		else if(target2 == "ND3") targ2NormCounts = ND3_Counts.getNormNt();
		else if(target2 == "CH2") targ2NormCounts = CH2_Counts.getNormNt();
		else if(target2 == "CD2") targ2NormCounts = CD2_Counts.getNormNt();
		else if(target2 == "C") targ2NormCounts = C_Counts.getNormNt();
		else if(target2 == "ET") targ2NormCounts = ET_Counts.getNormNt();
		else if(target2 == "F") targ2NormCounts = F_Counts.getNormNt();

		double scaleFactor1 = 1.0; double scaleFactor2 = 1.0;
		if( target1 == "ET" ) scaleFactor1 = ET_Scale_Factor( X_Avg_ET, Q2_Avg_ET );
		if( target2 == "ET" ) scaleFactor2 = ET_Scale_Factor( X_Avg_ET, Q2_Avg_ET );

		if( target1 == "CD2" && (Period == "Su22" || Period == "Fa22Neg" || Period == "Fa22Pos") ) scaleFactor1 *= ScaleCD2;
		if( target2 == "CD2" && (Period == "Su22" || Period == "Fa22Neg" || Period == "Fa22Pos") ) scaleFactor2 *= ScaleCD2;
	
		if( ( target1 == "ET" || target1 == "F" ) && Period == "Fa22Pos" ) scaleFactor1 *= ScaleSol;
		if( ( target2 == "ET" || target2 == "F" ) && Period == "Fa22Pos" ) scaleFactor2 *= ScaleSol;

		targ1NormCounts *= scaleFactor1; targ2NormCounts *= scaleFactor2;

		if( targ1NormCounts != 0 && targ2NormCounts != 0 ){
			return targ1NormCounts / targ2NormCounts;
		}

		return 0;
	}
	double getCountRatioErr(string target1, string target2, string Period){
		double targ1Counts = 0; double targ2Counts = 0;
		double targ1FCt = 0; double targ2FCt = 0;

		if(target1 == "NH3") targ1Counts = NH3_Counts.getNt();
		else if(target1 == "ND3") targ1Counts = ND3_Counts.getNt();
		else if(target1 == "CH2") targ1Counts = CH2_Counts.getNt();
		else if(target1 == "CD2") targ1Counts = CD2_Counts.getNt();
		else if(target1 == "C") targ1Counts = C_Counts.getNt();
		else if(target1 == "ET") targ1Counts = ET_Counts.getNt();
		else if(target1 == "F") targ1Counts = F_Counts.getNt();
		if(target2 == "NH3") targ2Counts = NH3_Counts.getNt();
		else if(target2 == "ND3") targ2Counts = ND3_Counts.getNt();
		else if(target2 == "CH2") targ2Counts = CH2_Counts.getNt();
		else if(target2 == "CD2") targ2Counts = CD2_Counts.getNt();
		else if(target2 == "C") targ2Counts = C_Counts.getNt();
		else if(target2 == "ET") targ2Counts = ET_Counts.getNt();
		else if(target2 == "F") targ2Counts = F_Counts.getNt();

		if(target1 == "NH3") targ1FCt = NH3_Counts.getFCt();
		else if(target1 == "ND3") targ1FCt = ND3_Counts.getFCt();
		else if(target1 == "CH2") targ1FCt = CH2_Counts.getFCt();
		else if(target1 == "CD2") targ1FCt = CD2_Counts.getFCt();
		else if(target1 == "C") targ1FCt = C_Counts.getFCt();
		else if(target1 == "ET") targ1FCt = ET_Counts.getFCt();
		else if(target1 == "F") targ1FCt = F_Counts.getFCt();
		if(target2 == "NH3") targ2FCt = NH3_Counts.getFCt();
		else if(target2 == "ND3") targ2FCt = ND3_Counts.getFCt();
		else if(target2 == "CH2") targ2FCt = CH2_Counts.getFCt();
		else if(target2 == "CD2") targ2FCt = CD2_Counts.getFCt();
		else if(target2 == "C") targ2FCt = C_Counts.getFCt();
		else if(target2 == "ET") targ2FCt = ET_Counts.getFCt();
		else if(target2 == "F") targ2FCt = F_Counts.getFCt();

		double scaleFactor1 = 1.0; double scaleFactor2 = 1.0;
		double targ1Err = sqrt( targ1Counts )/targ1FCt;
		double targ2Err = sqrt( targ2Counts )/targ2FCt;

		if( target1 == "ET" ) scaleFactor1 = ET_Scale_Factor( X_Avg_ET, Q2_Avg_ET );
		if( target2 == "ET" ) scaleFactor2 = ET_Scale_Factor( X_Avg_ET, Q2_Avg_ET );

		if( target1 == "CD2" && (Period == "Su22" || Period == "Fa22Neg" || Period == "Fa22Pos") ) scaleFactor1 *= ScaleCD2;
		if( target2 == "CD2" && (Period == "Su22" || Period == "Fa22Neg" || Period == "Fa22Pos") ) scaleFactor2 *= ScaleCD2;
	
		if( ( target1 == "ET" || target1 == "F" ) && Period == "Fa22Pos" ) scaleFactor1 *= ScaleSol;
		if( ( target2 == "ET" || target2 == "F" ) && Period == "Fa22Pos" ) scaleFactor2 *= ScaleSol;

		// Now that the scaling factors have been set, set the errors for targ1 and targ2
		if( ( target1 == "ET" || target1 == "F" ) && Period == "Fa22Pos" && targ1FCt != 0 ){
		    targ1Err = (1.0/targ1FCt)*sqrt( ScaleSol*ScaleSol*targ1Counts + (pow(targ1Counts,2))*pow(ScaleSolerr,2) );
		}
		if( ( target2 == "ET" || target2 == "F" ) && Period == "Fa22Pos" && targ2FCt != 0 ){
		    targ2Err = (1.0/targ2FCt)*sqrt( ScaleSol*ScaleSol*targ2Counts + (pow(targ2Counts,2))*pow(ScaleSolerr,2) );
		}

		if( target1 == "CD2" && ( Period == "Fa22Pos" || Period == "Su22" || Period == "Fa22Neg" ) && targ1FCt != 0 ){
		    targ1Err = (1.0/targ1FCt)*sqrt( ScaleCD2*ScaleCD2*targ1Counts + (pow(targ1Counts,2))*pow(ScaleCD2err,2) );
		}
		if( target2 == "CD2" && ( Period == "Fa22Pos" || Period == "Su22" || Period == "Fa22Neg" ) && targ2FCt != 0 ){
		    targ2Err = (1.0/targ2FCt)*sqrt( ScaleCD2*ScaleCD2*targ2Counts + (pow(targ2Counts,2))*pow(ScaleCD2err,2) );
		}

		if( targ1Counts != 0 && targ2Counts != 0 && targ1FCt != 0 && targ2FCt != 0 ){
			//double targ1Err = sqrt( targ1Counts )/targ1FCt;
			//double targ2Err = sqrt( targ2Counts )/targ2FCt;
			double targ1Norm = targ1Counts / targ1FCt; double targ2Norm = targ2Counts / targ2FCt;
			//return sqrt( targ1Err*targ1Err/(targ2Counts*targ2Counts) + targ2Err*targ2Err*targ1Counts*targ1Counts/pow(targ2Counts,4));
			return sqrt( ( pow(scaleFactor1,2)* targ1Norm/(targ2Norm*targ2Norm*targ1FCt)) + ((pow(scaleFactor2,2)*targ1Norm*targ1Norm)/(pow(targ2Norm,3)*targ2FCt)) );
		}
		
		return 0;

	}

	string ReturnBinContents() const{
		string object = to_string(Q2_Min)+"	"+to_string(Q2_Max)+"	";
		object += to_string(X_Min)+"	"+to_string(X_Max)+"	";
		object += to_string(Input_Loop.getNP()) +"	"+to_string(Input_Loop.getNM())+"	";
		object += to_string(Input_Loop.getFC_P())+"	"+to_string(Input_Loop.getFC_M())+"	";
		object += to_string(Input_Loop.getN0())+"	"+to_string(Input_Loop.getFC0());
		return object;
	}
	void Print( bool forcePrint=false ) const{
	    if( (Q2_Avg_All > 0 && X_Avg_All > 0) || forcePrint ){
		//cout << "=================================== Bin: Q^2 = "<<Q2_Mid <<" X = "<< X_Mid <<" ===================================\n";
		cout << "======== Bin: "<<Q2_Min <<" < Q2 < "<< Q2_Max <<", "<< X_Min <<" < X < "<< X_Max <<" ========\n";
		cout << "Statistics-Weighted Bins Per Target Type:\n";
		cout << "NH3:\n";
		cout << "--> Q^2 = "<< Q2_Avg_NH3 << endl;
		cout << "--> X   = "<< X_Avg_NH3 << endl;
		cout << "ND3:\n";
		cout << "--> Q^2 = "<< Q2_Avg_ND3 << endl;
		cout << "--> X   = "<< X_Avg_ND3 << endl;
		cout << "CH2:\n";
		cout << "--> Q^2 = "<< Q2_Avg_CH2 << endl;
		cout << "--> X   = "<< X_Avg_CH2 << endl;
		cout << "CD2:\n";
		cout << "--> Q^2 = "<< Q2_Avg_CD2 << endl;
		cout << "--> X   = "<< X_Avg_CD2 << endl;
		cout << "C:\n";
		cout << "--> Q^2 = "<< Q2_Avg_C << endl;
		cout << "--> X   = "<< X_Avg_C << endl;
		cout << "ET:\n";
		cout << "--> Q^2 = "<< Q2_Avg_ET << endl;
		cout << "--> X   = "<< X_Avg_ET << endl;
		cout << "F:\n";
		cout << "--> Q^2 = "<< Q2_Avg_F << endl;
		cout << "--> X   = "<< X_Avg_F << endl;
		cout << "All:\n";
		cout << "--> Q^2 = "<< Q2_Avg_All << endl;
		cout << "--> X   = "<< X_Avg_All << endl;
		NH3_Counts.Print(); ND3_Counts.Print();
		CH2_Counts.Print(); CD2_Counts.Print();
		C_Counts.Print(); ET_Counts.Print();
		F_Counts.Print();
		cout << "DF_NH3 = "<<DF_NH3 <<" +/- "<<ErrDF_NH3<<" +/- "<<SysErrDF_NH3<< endl;
		cout << "A1_NH3 = "<<A1_NH3 <<" +/- "<<ErrA1_NH3<< endl;
		cout << "DF_ND3 = "<<DF_ND3 <<" +/- "<<ErrDF_ND3<<" +/- "<<SysErrDF_ND3<< endl;
		cout << "A1_ND3 = "<<A1_ND3 <<" +/- "<<ErrA1_ND3<< endl;
		cout << "==================================================================================================================\n";
	    }
	}
	// Destructor
	~Bin() = default;

};

// This function is used in initializing the DataSet class to create a 2x2 matrix of bins
vector<vector<Bin>> SetBinWidths(){
	// In this 2x2 matrix, each column corresponds to one bin in Q2, and
	// each row is a bin in Bjorken x.
	vector<vector<Bin>> AllBins;
	//cout << "Now making the bin widths\n";
	for(size_t i=0; i<Q2_Bin_Bounds.size()-1; i++){ // Loop over Q2 bins
		double qmin = Q2_Bin_Bounds[i];
		double qmax = Q2_Bin_Bounds[i+1];
		vector<Bin> theseQ2Bins;
		//cout << qmin <<" "<< qmax;
		for(size_t j=0; j<X_Bin_Bounds.size()-1; j++){ // Loop over X bins
			Bin thisBin;
			double xmin = X_Bin_Bounds[j];
			double xmax = X_Bin_Bounds[j+1];
			//cout <<" "<< xmin <<" "<< xmax << endl;
			thisBin.SetQ2Bins(qmin, qmax);
			thisBin.SetXBins(xmin, xmax);
			theseQ2Bins.push_back(thisBin);
		}
		AllBins.push_back(theseQ2Bins);
	}
	// Now initialize the theoretical values for the double-spin asymmetry
	string fPath = THIS_DIR + "Asymmetry_Parameterizations.txt";
	//string fPath = THIS_DIR + "Old_Asymmetry_Parameterizations.txt";

	ifstream allin( fPath.c_str() );
	if( !allin.fail() ){
	    string line;
	    getline(allin,line); // Throw away header row
	    while( getline(allin, line) ){
	        stringstream sin(line);
		//double xbin, Q2bin, W2, Q2, xmid, Q2min, Q2max, Df_NH3, Err_Df_NH3, Q2_avg, W, nu, y, Eprime, theta, eps, D, eta, F1, F2, R, A1, A2, g1, g2;
		//sin >> xbin>>Q2bin>>W2>>Q2>>xmid>>Q2min>>Q2max>>Df_NH3>>Err_Df_NH3>>Q2_avg>>W>>nu>>y>>Eprime>>theta>>eps>>D>>eta>>F1>>F2>>R>>A1>>A2>>g1>>g2;
		//double Q2min, Q2max, qmid, Q2avg, Xmin, Xmax, xmid, Xavg,  W2,  F1,  F2,  R,  A1,  A2,  g1,  g2, D, eta;
		//sin >> Q2min >> Q2max >> qmid >> Q2avg >> Xmin >> Xmax >> xmid >> Xavg >>  W2 >>  F1 >>  F2 >>  R >>  A1 >>  A2 >>  g1 >>  g2 >> D >> eta;
		
		//Q2,	W2,	x,	F1,	F2,	R,	A1,	A2,	g1,	g2,
		double Q2, W2, x, F1, F2, R, A1, A2, g1, g2;
		sin >> Q2 >> W2 >> x >> F1 >> F2 >> R >> A1 >> A2 >> g1 >> g2;
		// Calculate the other relevant quantities
	        double Ep = avgBeamEnergy - ( (W2 + Q2 - (nucleon_mass*nucleon_mass)) / (2*nucleon_mass)   ); // Mean scattered energy for this Q2, W2
        	double theta = 2*asin( sqrt( Q2/(4*avgBeamEnergy*Ep) ) ); // Average scattered angle in radians
	        double tau = Q2 / (4*nucleon_mass*nucleon_mass*x*x); // tau kinematic variable
	        double eps = 1.0 / ( 1.0 + 2*(1.0+tau)*pow( tan(theta/2.0),2)  ); // virtual photon polarization epsilon, with theta passed in with radians
	        double depol = ( 1.0 - (eps*Ep/avgBeamEnergy) ) / ( 1.0 + (eps*R) ); // The depolarization factor
        	double eta = (eps * sqrt(Q2) ) / ( avgBeamEnergy - (eps*Ep) ); // The eta kinematic variable 

		//double Q2, W2, x, F1, F2, R, A1, A2, g1, g2, D, eta;
		//sin >> Q2 >> W2 >> x >> F1 >> F2 >> R >> A1 >> A2 >> g1 >> g2 >>D >> eta;

		//double qmid = Q2; //(Q2min+Q2max) / 2.0;
		//double xmid = x;
		//double qmid = (Q2min+Q2max) / 2.0;

		double ALL_Theory = depol*(A1+eta*A2); // Theoretical value of double spin asymmetry for this bin
		//cout << qmid <<" "<< xmid <<" "<< ALL_Theory << endl;
		for(size_t i=0; i<AllBins.size(); i++){ // Q2 bin loop
		    for(size_t j=0; j<AllBins[i].size(); j++){
		        if( AllBins[i][j].IsInBin( Q2, x )  ){
		            AllBins[i][j].SetAllTheory( ALL_Theory, "NH3" );
			    AllBins[i][j].SetKinematicFactors( depol, eta, A1, A2, "NH3" );
			    //AllBins[i][j].SetAllTheory( ALL_Theory_ND3, "ND3" );
			    // Placeholder for any future ND3 incorporation
		        }
		    }
		}
	    }
	}
	else{
	    cout << "Failed to read in A_ll theory values. Check path to 'Asymmetry_Parameterizations.txt' and try again.\n";
	}
	allin.close();

	return AllBins;
}


class DataSet{

    private:
	vector<vector<Bin>> AllBins;
	Bin NullBin; // Used to pass as an error parameter
	double PF_Avg = 0.0; // Average packing fraction for this data set
	double PF_Avg_Err = 0.0; // Statistical error from counts only
	double BT_Pol = 0.0; // Beam target polarization calculated from the max likelihood method using DIS data
	double BT_Pol_Err = 0.0; // Error on the above
    public:
	// Constructors
	DataSet() : AllBins( SetBinWidths() ){
	    // Automatically load in the scaling factors
	    this->SetScalingFactors();
	}

	// Mutators
	// The following function is used to construct an array of bins from input data
	//void AddToBins( string filePath, string targtype ){
	bool AddToBins( string filePath, string targtype ){
	    ifstream fin(filePath.c_str());
	    bool FCchargeIsSet = false; // Tracks whether or not the FC charge has been appended to the data set
	    if( !fin.fail() ){
		string line;
		getline(fin,line); // throw away header row
		//Q2_Min   Q2_Max   X_Min   X_Max   NP_Counts   NM_Counts   FC_P_Charge   FC_M_Charge
		while(getline(fin,line)){
		    stringstream sin(line);
		    double q2_min, q2_max, x_min, x_max, nm_counts, np_counts, fcm_charge, fcp_charge, n0_counts, fc0_charge;
		    sin >> q2_min >> q2_max >> x_min >> x_max >> nm_counts >> np_counts >> fcm_charge >> fcp_charge >> n0_counts >> fc0_charge;
		    double qmid = (q2_min+q2_max)/2.0; double xmid = (x_min+x_max)/2.0;
		    if( !FCchargeIsSet ){ 
			for(size_t i=0; i<AllBins.size(); i++){
			    for(size_t j=0; j<AllBins[i].size(); j++){
			        AllBins[i][j].AddFCCharge(fcm_charge, fcp_charge, fc0_charge, targtype);
			    }
			}
			FCchargeIsSet = true;
		    }
		    for(size_t i=0; i<AllBins.size(); i++){ // Q2 bin loop
		        for(size_t j=0; j<AllBins[i].size(); j++){
			    if( AllBins[i][j].IsInBin( qmid, xmid ) ){
			        AllBins[i][j].AddCounts( nm_counts, np_counts, n0_counts, targtype, xmid, qmid );
			    }
			}
		    }
		}
	    }
	    else{
		//cout << "Couldn't open file path " << filePath <<".\n";
		//cout << "Please check path and try again.\n";
	    }
	    fin.close();
	    return FCchargeIsSet; // If a file was successfully opened, this will be true.
	}

	// This is a special case where elastic data is being read in
	bool AddToBinsElastic( string filePath, string targtype ){
	    ifstream fin(filePath.c_str());
	    bool FCchargeIsSet = false; // Tracks whether or not the FC charge has been appended to the data set
	    if( !fin.fail() ){
		string line;
		getline(fin,line); // throw away header row
		//Q2_Min   Q2_Max   X_Min   X_Max   NP_Counts   NM_Counts   FC_P_Charge   FC_M_Charge
		while(getline(fin,line)){
		    stringstream sin(line);
		    //3.38465  11.0347  0.225091  52  41  37242.6  37201
		    double q2_avg, theta_avg, A_el, Np_Counts, Nm_Counts, FCp, FCm;
		    sin >> q2_avg >> theta_avg >> A_el >> Np_Counts >> Nm_Counts >> FCp >> FCm;
		    // Because the data is elastic, x=1, so an arbitrary xbin stores the data
		    double xBin = (X_Bin_Bounds[0] + X_Bin_Bounds[1]) / 2.0;
		    if( !FCchargeIsSet ){ 
			for(size_t i=0; i<AllBins.size(); i++){
			    for(size_t j=0; j<AllBins[i].size(); j++){
			        AllBins[i][j].AddFCCharge(FCp, FCm, 0, targtype);
			    }
			}
			FCchargeIsSet = true;
		    }
		    for(size_t i=0; i<AllBins.size(); i++){ // Q2 bin loop
		        for(size_t j=0; j<AllBins[i].size(); j++){
			    if( AllBins[i][j].IsInBin( q2_avg, xBin ) ){
			        AllBins[i][j].AddCounts( Np_Counts, Nm_Counts, 0, targtype, xBin, q2_avg );
			    }
			}
		    }
		}
	    }
	    else{
		//cout << "Couldn't open file path " << filePath <<".\n";
		//cout << "Please check path and try again.\n";
	    }
	    fin.close();
	    return FCchargeIsSet; // If a file was successfully opened, this will be true.
	}

	// This function sets all the relevant info for a given run period and a given epoch.
	// If the boolean "useElasticData" is true, then files will be read in from the "Elastic_Text_Files"
	// directory to use Noemie's data.
	// FIXME: Add functionality for reading in the CD2 data and the scaled CD2 data for the ND3 calculations...
	void SetEpochValues( string Period, string Target, vector<int>& Epoch, vector<int>& MissedRuns, bool useElasticData=false ){
		
	    // Select the background runs that will be used based on the run period
	    vector<int> C_Bkg, CH2_Bkg, ET_Bkg, F_Bkg, CD2_Bkg;
	    if( Period == "Su22" ){
		C_Bkg = Su22_C_Runs; CH2_Bkg = Su22_CH2_Runs; ET_Bkg = Su22_ET_Runs; F_Bkg = Su22_F_Runs;
	    }
	    else if( Period == "Fa22Neg" ){
		C_Bkg = Fa22Neg_C_Runs; CH2_Bkg = Fa22Neg_CH2_Runs; ET_Bkg = Fa22Neg_ET_Runs; F_Bkg = Fa22Neg_F_Runs;
	    }
	    else if( Period == "Fa22Pos" ){
		C_Bkg = Fa22Pos_C_Runs; CH2_Bkg = Fa22Pos_CH2_Runs; 
		// Scaling factors aren't applied to these runs, but it should cancel out in the calculation anyway...
		ET_Bkg = Fa22Neg_ET_Runs; 
		F_Bkg = Fa22Neg_F_Runs;
	    }
	    else if( Period == "Sp23Inb" ){
		C_Bkg = Sp23Inb_C_Runs; CH2_Bkg = Sp23Inb_CH2_Runs; ET_Bkg = Sp23Inb_ET_Runs; F_Bkg = Sp23Inb_F_Runs;
		CD2_Bkg = Sp23Inb_CD2_Runs;
	    }
	    else{
		cout <<"ERROR: Invalid run period "<< Period <<". Valid options are: 'Su22' 'Fa22Neg' 'Fa22Pos' 'Sp23Inb'\n";
		cout <<"       No data was read in!\n";
		return; // Kills the function
	    }

	    // Now add the data from the background runs and the given epoch vector
	    // NH3/ND3 data
	    for(int run : Epoch ){
		//string filePath = "/w/hallb-scshelf2102/clas12/holmberg/Spin_Structure/Latest_Skims/Text_Files/"+ Target +"_"+ to_string(run) +"_DF_Data.txt";
		string filePath;
		if( !useElasticData ){
		    filePath = "/w/hallb-scshelf2102/clas12/holmberg/Spin_Structure/Latest_Skims/Text_Files/"+ Target +"_"+ to_string(run) +"_DF_Data.txt";
		    if( !this->AddToBins( filePath, Target ) ) MissedRuns.push_back( run );
		}
		else{
		    filePath = "/w/hallb-scshelf2102/clas12/holmberg/Spin_Structure/Latest_Skims/Elastic_Text_Files/"+ Target +"_"+ to_string(run) +"_Elastic.txt";
		    if( !this->AddToBinsElastic( filePath, Target ) ) MissedRuns.push_back( run );
		}
	    }
    	    // CH2 data
	    for(int run : CH2_Bkg ){
		//string filePath = "/w/hallb-scshelf2102/clas12/holmberg/Spin_Structure/Latest_Skims/Text_Files/CH2_"+ to_string(run) +"_DF_Data.txt";
		//if( !this->AddToBins( filePath, "CH2" ) ) MissedRuns.push_back( run );
		string filePath;
		if( !useElasticData ){
		    filePath = "/w/hallb-scshelf2102/clas12/holmberg/Spin_Structure/Latest_Skims/Text_Files/CH2_"+ to_string(run) +"_DF_Data.txt";
		    if( !this->AddToBins( filePath, "CH2" ) ) MissedRuns.push_back( run );
		}
		else{
		    filePath = "/w/hallb-scshelf2102/clas12/holmberg/Spin_Structure/Latest_Skims/Elastic_Text_Files/CH2_"+ to_string(run) +"_Elastic.txt";
		    if( !this->AddToBinsElastic( filePath, "CH2" ) ) MissedRuns.push_back( run );
		}
	    }
	    // Carbon data
	    for(int run : C_Bkg ){
		//string filePath = "/w/hallb-scshelf2102/clas12/holmberg/Spin_Structure/Latest_Skims/Text_Files/C_"+ to_string(run) +"_DF_Data.txt";
		//if( !this->AddToBins( filePath, "C" ) ) MissedRuns.push_back( run );
		string filePath;
		if( !useElasticData ){
		    filePath = "/w/hallb-scshelf2102/clas12/holmberg/Spin_Structure/Latest_Skims/Text_Files/C_"+ to_string(run) +"_DF_Data.txt";
		    if( !this->AddToBins( filePath, "C" ) ) MissedRuns.push_back( run );
		}
		else{
		    filePath = "/w/hallb-scshelf2102/clas12/holmberg/Spin_Structure/Latest_Skims/Elastic_Text_Files/C_"+ to_string(run) +"_Elastic.txt";
		    if( !this->AddToBinsElastic( filePath, "C" ) ) MissedRuns.push_back( run );
		}
	    }
	    // ET data
	    for(int run : ET_Bkg ){
		//string filePath = "/w/hallb-scshelf2102/clas12/holmberg/Spin_Structure/Latest_Skims/Text_Files/ET_"+ to_string(run) +"_DF_Data.txt";
		//if( !this->AddToBins( filePath, "ET" ) ) MissedRuns.push_back( run );
		string filePath;
		if( !useElasticData ){
		    filePath = "/w/hallb-scshelf2102/clas12/holmberg/Spin_Structure/Latest_Skims/Text_Files/ET_"+ to_string(run) +"_DF_Data.txt";
		    if( !this->AddToBins( filePath, "ET" ) ) MissedRuns.push_back( run );
		}
		else{
		    filePath = "/w/hallb-scshelf2102/clas12/holmberg/Spin_Structure/Latest_Skims/Elastic_Text_Files/ET_"+ to_string(run) +"_Elastic.txt";
		    if( !this->AddToBinsElastic( filePath, "ET" ) ) MissedRuns.push_back( run );
		}
	    }
	    // Foil data
	    for(int run : F_Bkg ){
		//string filePath = "/w/hallb-scshelf2102/clas12/holmberg/Spin_Structure/Latest_Skims/Text_Files/F_"+ to_string(run) +"_DF_Data.txt";
		//if( !this->AddToBins( filePath, "F" ) ) MissedRuns.push_back( run );
		string filePath;
		if( !useElasticData ){
		    filePath = "/w/hallb-scshelf2102/clas12/holmberg/Spin_Structure/Latest_Skims/Text_Files/F_"+ to_string(run) +"_DF_Data.txt";
		    if( !this->AddToBins( filePath, "F" ) ) MissedRuns.push_back( run );
		}
		else{
		    filePath = "/w/hallb-scshelf2102/clas12/holmberg/Spin_Structure/Latest_Skims/Elastic_Text_Files/F_"+ to_string(run) +"_Elastic.txt";
		    if( !this->AddToBinsElastic( filePath, "F" ) ) MissedRuns.push_back( run );
		}
	    }
/*	    // CD2 data
	    for(int run : CD2_Bkg ){
		//string filePath = "/w/hallb-scshelf2102/clas12/holmberg/Spin_Structure/Latest_Skims/Text_Files/"+ Target +"_"+ to_string(run) +"_DF_Data.txt";
		//if( !this->AddToBins( filePath, "CD2" ) ) MissedRuns.push_back( run );
		string filePath;
		if( !useElasticData ){
		    filePath = "/w/hallb-scshelf2102/clas12/holmberg/Spin_Structure/Latest_Skims/Text_Files/"+ Target +"_"+ to_string(run) +"_DF_Data.txt";
		    if( !this->AddToBins( filePath, "CD2" ) ) MissedRuns.push_back( run );
		}
		else{
		    filePath = "/w/hallb-scshelf2102/clas12/holmberg/Spin_Structure/Latest_Skims/Elastic_Text_Files/"+ Target +"_"+ to_string(run) +"_Elastic.txt";
		    if( !this->AddToBinsElastic( filePath, "CD2" ) ) MissedRuns.push_back( run );
		}
	    }*/
	   
	    if( MissedRuns.size() > 0 ){
		cout <<"ERROR: Missed runs in the epoch or background runs for "<< Period <<". Missing runs:\n";
		for(int r : MissedRuns) cout <<"--> "<< r << endl;
		cout << endl;
	    } 
	    else cout << "Finished reading in data for epoch.\n";
	}

        // "Period" is the run period ("Su22", "Fa22Pos", etc.), and the vector is all the NH3 or ND3 runs in the epoch.
        // To calculate the average values across all runs
	void CalculateAvgXQ2( string Period, vector<int>& Epoch, string Target, bool useElasticData=false ){

	    // First, read in for each target type
	    // Number of x and Q2 bins
	    int nQ2Bins = Q2_Bin_Bounds.size()-1;
	    int nXBins  = X_Bin_Bounds.size()-1;

	    // Select the background runs that will be used based on the run period
	    vector<int> C_Bkg, CH2_Bkg, ET_Bkg, F_Bkg, CD2_Bkg;
	    if( Period == "Su22" ){
		C_Bkg = Su22_C_Runs; CH2_Bkg = Su22_CH2_Runs; ET_Bkg = Su22_ET_Runs; F_Bkg = Su22_F_Runs;
	    }
	    else if( Period == "Fa22Neg" ){
		C_Bkg = Fa22Neg_C_Runs; CH2_Bkg = Fa22Neg_CH2_Runs; ET_Bkg = Fa22Neg_ET_Runs; F_Bkg = Fa22Neg_F_Runs;
	    }
	    else if( Period == "Fa22Pos" ){
		C_Bkg = Fa22Pos_C_Runs; CH2_Bkg = Fa22Pos_CH2_Runs; 
		// Scaling factors aren't applied to these runs, but it should cancel out in the calculation anyway...
		ET_Bkg = Fa22Neg_ET_Runs; 
		F_Bkg = Fa22Neg_F_Runs;
	    }
	    else if( Period == "Sp23Inb" ){
		C_Bkg = Sp23Inb_C_Runs; CH2_Bkg = Sp23Inb_CH2_Runs; ET_Bkg = Sp23Inb_ET_Runs; F_Bkg = Sp23Inb_F_Runs;
		CD2_Bkg = Sp23Inb_CD2_Runs;
	    }
	    else{
		cout <<"ERROR: Invalid run period "<< Period <<". Valid options are: 'Su22' 'Fa22Neg' 'Fa22Pos' 'Sp23Inb'\n";
		cout <<"       No average values of X, Q2 were set!\n";
		return; // Kills the function
	    }

	    // This section of the code is called for setting the bins in the case of elastic scattering data.
	    // FIXME: Add CD2 loop to this
	    /*************************************** ELASTIC DATA EVALUATION *****************************************/
	    if( useElasticData ){
		// Number of elastic Q2 bins to consider
		double nQ2bins = Q2_Bin_Bounds.size()-1;
		// Get the total number of counts to take the stat-weighted average of the kinematics from each run separately
		vector<double> Q2_All(nQ2bins, 0); vector<double> Q2_Ammonia(nQ2bins, 0); vector<double> Q2_C  (nQ2bins, 0); vector<double> Q2_CH2(nQ2bins, 0);
		vector<double> Q2_ET(nQ2bins, 0);  vector<double> Q2_F(nQ2bins, 0);       vector<double> Q2_CD2(nQ2bins, 0);

		vector<double> Counts_All(nQ2bins, 0); vector<double> Counts_Ammonia(nQ2bins, 0); vector<double> Counts_C  (nQ2bins, 0); vector<double> Counts_CH2(nQ2bins, 0);
		vector<double> Counts_ET(nQ2bins, 0);  vector<double> Counts_F(nQ2bins, 0);       vector<double> Counts_CD2(nQ2bins, 0);
		// First loop over the NH3/ND3 data
		for(int run : Epoch){
		    // Loop thru each kinematic bin once the file is opened
		    ifstream fin( string("../Latest_Skims/Elastic_Text_Files/"+ Target +"_"+ to_string(run) +"_Elastic.txt") );
		    if( fin.fail() ){ cout <<"ERROR: Couldn't find elastic file for "<< Target <<" run "<< run <<". Aborting Avg x, Q2 calculation...\n"; return; }
		    string line;
		    getline(fin,line); // Throw away header row
		    //Q2_Mean   Theta_Mean   A_el   N+   N-   FC+   FC-
		    int it = 0; // iterator for tracking the vectors
		    while( getline( fin, line ) ){
			double Q2_Mean, Theta_Mean, A_el, Np, Nm, FCp, FCm;
			stringstream sin(line);
			sin >> Q2_Mean >> Theta_Mean >> A_el >> Np >> Nm >> FCp >> FCm;
			if( Np + Nm > 0 ){
			    Q2_All[it] += (Nm + Np)*Q2_Mean; Q2_Ammonia[it] += (Nm + Np)*Q2_Mean;
			    Counts_All[it] += (Nm + Np); Counts_Ammonia[it] += (Nm + Np);
			}
			it++;
		    }
		    fin.close();
		}
		// Loop over CH2 data
		for(int run : CH2_Bkg){
		    // Loop thru each kinematic bin once the file is opened
		    ifstream fin( string("../Latest_Skims/Elastic_Text_Files/CH2_"+ to_string(run) +"_Elastic.txt") );
		    if( fin.fail() ){ cout <<"ERROR: Couldn't find elastic file for CH2 run "<< run <<". Aborting Avg x, Q2 calculation...\n"; return; }
		    string line;
		    getline(fin,line); // Throw away header row
		    //Q2_Mean   Theta_Mean   A_el   N+   N-   FC+   FC-
		    int it = 0; // iterator for tracking the vectors
		    while( getline( fin, line ) ){
			double Q2_Mean, Theta_Mean, A_el, Np, Nm, FCp, FCm;
			stringstream sin(line);
			sin >> Q2_Mean >> Theta_Mean >> A_el >> Np >> Nm >> FCp >> FCm;
			if( Np + Nm > 0 ){
			    Q2_All[it] += (Nm + Np)*Q2_Mean; Q2_CH2[it] += (Nm + Np)*Q2_Mean;
			    Counts_All[it] += (Nm + Np); Counts_CH2[it] += (Nm + Np);
			}
			it++;
		    }
		    fin.close();
		}
		// Loop over C data
		for(int run : C_Bkg){
		    // Loop thru each kinematic bin once the file is opened
		    ifstream fin( string("../Latest_Skims/Elastic_Text_Files/C_"+ to_string(run) +"_Elastic.txt") );
		    if( fin.fail() ){ cout <<"ERROR: Couldn't find elastic file for C run "<< run <<". Aborting Avg x, Q2 calculation...\n"; return; }
		    string line;
		    getline(fin,line); // Throw away header row
		    //Q2_Mean   Theta_Mean   A_el   N+   N-   FC+   FC-
		    int it = 0; // iterator for tracking the vectors
		    while( getline( fin, line ) ){
			double Q2_Mean, Theta_Mean, A_el, Np, Nm, FCp, FCm;
			stringstream sin(line);
			sin >> Q2_Mean >> Theta_Mean >> A_el >> Np >> Nm >> FCp >> FCm;
			if( Np + Nm > 0 ){
			    Q2_All[it] += (Nm + Np)*Q2_Mean; Q2_C[it] += (Nm + Np)*Q2_Mean;
			    Counts_All[it] += (Nm + Np); Counts_C[it] += (Nm + Np);
			}
			it++;
		    }
		    fin.close();
		}
		// Loop over ET data
		for(int run : ET_Bkg){
		    // Loop thru each kinematic bin once the file is opened
		    ifstream fin( string("../Latest_Skims/Elastic_Text_Files/ET_"+ to_string(run) +"_Elastic.txt") );
		    if( fin.fail() ){ cout <<"ERROR: Couldn't find elastic file for ET run "<< run <<". Aborting Avg x, Q2 calculation...\n"; return; }
		    string line;
		    getline(fin,line); // Throw away header row
		    //Q2_Mean   Theta_Mean   A_el   N+   N-   FC+   FC-
		    int it = 0; // iterator for tracking the vectors
		    while( getline( fin, line ) ){
			double Q2_Mean, Theta_Mean, A_el, Np, Nm, FCp, FCm;
			stringstream sin(line);
			sin >> Q2_Mean >> Theta_Mean >> A_el >> Np >> Nm >> FCp >> FCm;
			if( Np + Nm > 0 ){
			    Q2_All[it] += (Nm + Np)*Q2_Mean; Q2_ET[it] += (Nm + Np)*Q2_Mean;
			    Counts_All[it] += (Nm + Np); Counts_ET[it] += (Nm + Np);
			}
			it++;
		    }
		    fin.close();
		}
		// Loop over F data
		for(int run : F_Bkg){
		    // Loop thru each kinematic bin once the file is opened
		    ifstream fin( string("../Latest_Skims/Elastic_Text_Files/F_"+ to_string(run) +"_Elastic.txt") );
		    if( fin.fail() ){ cout <<"ERROR: Couldn't find elastic file for F run "<< run <<". Aborting Avg x, Q2 calculation...\n"; return; }
		    string line;
		    getline(fin,line); // Throw away header row
		    //Q2_Mean   Theta_Mean   A_el   N+   N-   FC+   FC-
		    int it = 0; // iterator for tracking the vectors
		    while( getline( fin, line ) ){
			double Q2_Mean, Theta_Mean, A_el, Np, Nm, FCp, FCm;
			stringstream sin(line);
			sin >> Q2_Mean >> Theta_Mean >> A_el >> Np >> Nm >> FCp >> FCm;
			if( Np + Nm > 0 ){
			    Q2_All[it] += (Nm + Np)*Q2_Mean; Q2_F[it] += (Nm + Np)*Q2_Mean;
			    Counts_All[it] += (Nm + Np); Counts_F[it] += (Nm + Np);
			}
			it++;
		    }
		    fin.close();
		}

		// Now set the average values.
		// For elastics, I just set to the first x bin
		double xBin = (X_Bin_Bounds[0] + X_Bin_Bounds[1]) / 2.0;
		for( int q=0; q<nQ2Bins; q++ ){
		    if( Counts_All[q]     !=0 ){ double avgQ2_All     = Q2_All[q] / Counts_All[q];         this->SetAvgXQ2( "All", xBin, avgQ2_All ); }
		    if( Counts_Ammonia[q] !=0 ){ double avgQ2_Ammonia = Q2_Ammonia[q] / Counts_Ammonia[q]; this->SetAvgXQ2( Target, xBin, avgQ2_Ammonia ); }
		    if( Counts_CH2[q]     !=0 ){ double avgQ2_CH2     = Q2_CH2[q] / Counts_CH2[q];         this->SetAvgXQ2( "CH2", xBin, avgQ2_CH2 ); }
		    if( Counts_C[q]       !=0 ){ double avgQ2_C       = Q2_C[q] / Counts_C[q];             this->SetAvgXQ2( "C", xBin, avgQ2_C ); }
		    if( Counts_ET[q]      !=0 ){ double avgQ2_ET      = Q2_ET[q] / Counts_ET[q];           this->SetAvgXQ2( "ET", xBin, avgQ2_ET ); }
		    if( Counts_F[q]       !=0 ){ double avgQ2_F       = Q2_F[q] / Counts_F[q];             this->SetAvgXQ2( "F", xBin, avgQ2_F ); }
		    //if( Counts_CD2[q]   !=0 ){ double avgQ2_CD2     = Q2_CD2[q] / Counts_CD2[q];         this->SetAvgXQ2( "CD2", xBin, avgQ2_CD2 ); }

		}
		cout << "Finished reading all elastic data for setting average Q2 values.\n";
		return; // Kills the function to ignore the DIS stuff below...
	    }

	    /****************************************** DIS DATA EVALUATION ******************************************/

	    // This stores the average value of Q2 for each of the x, Q2 bins across ALL target types in the files
	    TProfile2D* Q2BinData = new TProfile2D({"Q2BinData","Q2 Bin Data; Q2 [GeV^{2}]; X", nQ2Bins, Q2_Bin_Bounds.data(), nXBins, X_Bin_Bounds.data()});
	    TProfile2D* XBinData  = new TProfile2D({"XBinData","X Bin Data; Q2 [GeV^{2}]; X", nQ2Bins, Q2_Bin_Bounds.data(), nXBins, X_Bin_Bounds.data()});
	    // Now create profiles for each target type
	    // Carbon
	    TProfile2D* Q2BinData_C = new TProfile2D({"Q2BinData_C","Q2 Bin Data; Q2 [GeV^{2}]; X", nQ2Bins, Q2_Bin_Bounds.data(), nXBins, X_Bin_Bounds.data()});
	    TProfile2D* XBinData_C  = new TProfile2D({"XBinData_C","X Bin Data; Q2 [GeV^{2}]; X", nQ2Bins, Q2_Bin_Bounds.data(), nXBins, X_Bin_Bounds.data()});
	    for( int Run : C_Bkg ){
		string filePath = "/w/hallb-scshelf2102/clas12/holmberg/Spin_Structure/Latest_Skims/ROOT_Files/C_"+ to_string(Run) +"_Data.root";
		ifstream test( filePath );
		if( test.good() ){ // Proceed for a valid file path
		    TFile* Input = new TFile( filePath.c_str() );
		    auto q2In = (TProfile2D*)Input->Get("Q2BinData"); Q2BinData->Add( q2In ); Q2BinData_C->Add( q2In );
		    auto xIn  = (TProfile2D*)Input->Get("XBinData");  XBinData->Add( xIn );   XBinData_C->Add( xIn );
		    Input->Close();
		}
		else{ // Kill the operation if one of the runs is missing!!
		    cout <<"ERROR: Unable to locate ROOT file:\n";
		    cout <<"       "<< filePath << endl;
		    cout <<"       No average values of X, Q2 were set!\n";
		    return; // Kills the function
		}
	    }

	    // CH2
	    TProfile2D* Q2BinData_CH2 = new TProfile2D({"Q2BinData_CH2","Q2 Bin Data; Q2 [GeV^{2}]; X", nQ2Bins, Q2_Bin_Bounds.data(), nXBins, X_Bin_Bounds.data()});
	    TProfile2D* XBinData_CH2  = new TProfile2D({"XBinData_CH2","X Bin Data; Q2 [GeV^{2}]; X", nQ2Bins, Q2_Bin_Bounds.data(), nXBins, X_Bin_Bounds.data()});
	    for( int Run : CH2_Bkg ){
		string filePath = "/w/hallb-scshelf2102/clas12/holmberg/Spin_Structure/Latest_Skims/ROOT_Files/CH2_"+ to_string(Run) +"_Data.root";
		ifstream test( filePath );
		if( test.good() ){ // Proceed for a valid file path
		    TFile* Input = new TFile( filePath.c_str() );
		    auto q2In = (TProfile2D*)Input->Get("Q2BinData"); Q2BinData->Add( q2In ); Q2BinData_CH2->Add( q2In );
		    auto xIn  = (TProfile2D*)Input->Get("XBinData");  XBinData->Add( xIn );   XBinData_CH2->Add( xIn );
		    Input->Close();
		}
		else{ // Kill the operation if one of the runs is missing!!
		    cout <<"ERROR: Unable to locate ROOT file:\n";
		    cout <<"       "<< filePath << endl;
		    cout <<"       No average values of X, Q2 were set!\n";
		    return; // Kills the function
		}
	    }

	    // Empty
	    TProfile2D* Q2BinData_ET = new TProfile2D({"Q2BinData_ET","Q2 Bin Data; Q2 [GeV^{2}]; X", nQ2Bins, Q2_Bin_Bounds.data(), nXBins, X_Bin_Bounds.data()});
	    TProfile2D* XBinData_ET  = new TProfile2D({"XBinData_ET","X Bin Data; Q2 [GeV^{2}]; X", nQ2Bins, Q2_Bin_Bounds.data(), nXBins, X_Bin_Bounds.data()});
	    for( int Run : ET_Bkg ){
		string filePath = "/w/hallb-scshelf2102/clas12/holmberg/Spin_Structure/Latest_Skims/ROOT_Files/Empty_"+ to_string(Run) +"_Data.root";
		ifstream test( filePath );
		if( test.good() ){ // Proceed for a valid file path
		    TFile* Input = new TFile( filePath.c_str() );
		    auto q2In = (TProfile2D*)Input->Get("Q2BinData"); Q2BinData->Add( q2In ); Q2BinData_ET->Add( q2In );
		    auto xIn  = (TProfile2D*)Input->Get("XBinData");  XBinData->Add( xIn );   XBinData_ET->Add( xIn );
		    Input->Close();
		}
		else{ // Kill the operation if one of the runs is missing!!
		    cout <<"ERROR: Unable to locate ROOT file:\n";
		    cout <<"       "<< filePath << endl;
		    cout <<"       No average values of X, Q2 were set!\n";
		    return; // Kills the function
		}
	    }

	    // Foil
	    TProfile2D* Q2BinData_F = new TProfile2D({"Q2BinData_F","Q2 Bin Data; Q2 [GeV^{2}]; X", nQ2Bins, Q2_Bin_Bounds.data(), nXBins, X_Bin_Bounds.data()});
	    TProfile2D* XBinData_F  = new TProfile2D({"XBinData_F","X Bin Data; Q2 [GeV^{2}]; X", nQ2Bins, Q2_Bin_Bounds.data(), nXBins, X_Bin_Bounds.data()});
	    for( int Run : F_Bkg ){
		string filePath = "/w/hallb-scshelf2102/clas12/holmberg/Spin_Structure/Latest_Skims/ROOT_Files/Foil_"+ to_string(Run) +"_Data.root";
		ifstream test( filePath );
		if( test.good() ){ // Proceed for a valid file path
		    TFile* Input = new TFile( filePath.c_str() );
		    auto q2In = (TProfile2D*)Input->Get("Q2BinData"); Q2BinData->Add( q2In ); Q2BinData_F->Add( q2In );
		    auto xIn  = (TProfile2D*)Input->Get("XBinData");  XBinData->Add( xIn );   XBinData_F->Add( xIn );
		    Input->Close();
		}
		else{ // Kill the operation if one of the runs is missing!!
		    cout <<"ERROR: Unable to locate ROOT file:\n";
		    cout <<"       "<< filePath << endl;
		    cout <<"       No average values of X, Q2 were set!\n";
		    return; // Kills the function
		}
	    }

	    // CD2
	    TProfile2D* Q2BinData_CD2 = new TProfile2D({"Q2BinData_CD2","Q2 Bin Data; Q2 [GeV^{2}]; X", nQ2Bins, Q2_Bin_Bounds.data(), nXBins, X_Bin_Bounds.data()});
	    TProfile2D* XBinData_CD2  = new TProfile2D({"XBinData_CD2","X Bin Data; Q2 [GeV^{2}]; X", nQ2Bins, Q2_Bin_Bounds.data(), nXBins, X_Bin_Bounds.data()});
	    for( int Run : CD2_Bkg ){
		string filePath = "/w/hallb-scshelf2102/clas12/holmberg/Spin_Structure/Latest_Skims/ROOT_Files/CD2_"+ to_string(Run) +"_Data.root";
		ifstream test( filePath );
		if( test.good() ){ // Proceed for a valid file path
		    TFile* Input = new TFile( filePath.c_str() );
		    auto q2In = (TProfile2D*)Input->Get("Q2BinData"); Q2BinData->Add( q2In ); Q2BinData_CD2->Add( q2In );
		    auto xIn  = (TProfile2D*)Input->Get("XBinData");  XBinData->Add( xIn );   XBinData_CD2->Add( xIn );
		    Input->Close();
		}
		else{ // Kill the operation if one of the runs is missing!!
		    cout <<"ERROR: Unable to locate ROOT file:\n";
		    cout <<"       "<< filePath << endl;
		    cout <<"       No average values of X, Q2 were set!\n";
		    return; // Kills the function
		}
	    }
	    // Ammonia
	    TProfile2D* Q2BinData_Ammonia = new TProfile2D({"Q2BinData_Ammonia","Q2 Bin Data; Q2 [GeV^{2}]; X", nQ2Bins, Q2_Bin_Bounds.data(), nXBins, X_Bin_Bounds.data()});
	    TProfile2D* XBinData_Ammonia  = new TProfile2D({"XBinData_Ammonia","X Bin Data; Q2 [GeV^{2}]; X", nQ2Bins, Q2_Bin_Bounds.data(), nXBins, X_Bin_Bounds.data()});
	    for( int Run : Epoch ){
		string filePath = "/w/hallb-scshelf2102/clas12/holmberg/Spin_Structure/Latest_Skims/ROOT_Files/"+ Target +"_"+ to_string(Run) +"_Data.root";
		ifstream test( filePath );
		if( test.good() ){ // Proceed for a valid file path
		    TFile* Input = new TFile( filePath.c_str() );
		    auto q2In = (TProfile2D*)Input->Get("Q2BinData"); Q2BinData->Add( q2In ); Q2BinData_Ammonia->Add( q2In );
		    auto xIn  = (TProfile2D*)Input->Get("XBinData");  XBinData->Add( xIn );   XBinData_Ammonia->Add( xIn );
		    Input->Close();
		}
		else{ // Kill the operation if one of the runs is missing!!
		    cout <<"ERROR: Unable to locate ROOT file:\n";
		    cout <<"       "<< filePath << endl;
		    cout <<"       No average values of X, Q2 were set!\n";
		    return; // Kills the function
		}
	    }

	    // If everything was read in successfully, then set the average values for each target type:

	    for( int i=0; i<nQ2Bins; i++ ){
		for( int j=0; j<nXBins; j++ ){

		    int binIndex = XBinData->GetBin(i,j);
		    // Get the average x and Q2 values for this bin for each target type:
		    
		    // All target types
		    double All_X  = XBinData->GetBinContent(  binIndex );
		    double All_Q2 = Q2BinData->GetBinContent( binIndex );
		    this->SetAvgXQ2( "All", All_X, All_Q2 );

		    // C target types
		    double C_X  = XBinData_C->GetBinContent(  binIndex );
		    double C_Q2 = Q2BinData_C->GetBinContent( binIndex );
		    this->SetAvgXQ2( "C", C_X, C_Q2 );

		    // CH2 target types
		    double CH2_X  = XBinData_CH2->GetBinContent(  binIndex );
		    double CH2_Q2 = Q2BinData_CH2->GetBinContent( binIndex );
		    this->SetAvgXQ2( "CH2", CH2_X, CH2_Q2 );

		    // ET target types
		    double ET_X  = XBinData_ET->GetBinContent(  binIndex );
		    double ET_Q2 = Q2BinData_ET->GetBinContent( binIndex );
		    this->SetAvgXQ2( "ET", ET_X, ET_Q2 );

		    // F target types
		    double F_X  = XBinData_F->GetBinContent(  binIndex );
		    double F_Q2 = Q2BinData_F->GetBinContent( binIndex );
		    this->SetAvgXQ2( "F", F_X, F_Q2 );

		    // CD2 target types
		    double CD2_X  = XBinData_CD2->GetBinContent(  binIndex );
		    double CD2_Q2 = Q2BinData_CD2->GetBinContent( binIndex );
		    this->SetAvgXQ2( "CD2", CD2_X, CD2_Q2 );

		    // Ammonia
		    double Ammonia_X  = XBinData_Ammonia->GetBinContent(  binIndex );
		    double Ammonia_Q2 = Q2BinData_Ammonia->GetBinContent( binIndex );
		    this->SetAvgXQ2( Target, Ammonia_X, Ammonia_Q2 );

		} // End loop on X bins
	    } // End loop on Q2 bins

	    // Delete the profiles to prevent memory leaks...
	    XBinData->Delete();         Q2BinData->Delete();
	    XBinData_C->Delete();       Q2BinData_C->Delete();
	    XBinData_CH2->Delete();     Q2BinData_CH2->Delete();
	    XBinData_CD2->Delete();     Q2BinData_CD2->Delete();
	    XBinData_ET->Delete();      Q2BinData_ET->Delete();
	    XBinData_F->Delete();       Q2BinData_F->Delete();
	    XBinData_Ammonia->Delete(); Q2BinData_Ammonia->Delete();

	    cout << "Average X, Q2 set successfully!\n";
	}

        // This version of "CalculateAvgXQ2()" calculates the average kinematics for only the target
        // types included in the "Targets" vector. For example, if a user put in "CH2" and "CD2" for the 
        // vector, then the calculation for the "All" category (across all target types) would be
        // calculated for the "CH2" and "CD2" data. This is really only used with the scaling factor calculation
        // right now...
	void CalculateAvgXQ2( string Period, vector<string>& Targets ){

	    if( Targets.size() == 0 ){
		cout <<"ERROR: No targets in the vector. Check inputs.\n";
		cout <<"       No average values of X, Q2 were set!\n";
		return; // Kills the function
	    }
	    // First, read in for each target type
	    // Number of x and Q2 bins
	    int nQ2Bins = Q2_Bin_Bounds.size()-1;
	    int nXBins  = X_Bin_Bounds.size()-1;

	    // Select the background runs that will be used based on the run period
	    vector<int> C_Bkg, CH2_Bkg, ET_Bkg, F_Bkg, CD2_Bkg;
	    if( Period == "Su22" ){
		C_Bkg = Su22_C_Runs; CH2_Bkg = Su22_CH2_Runs; ET_Bkg = Su22_ET_Runs; F_Bkg = Su22_F_Runs;
	    }
	    else if( Period == "Fa22Neg" ){
		C_Bkg = Fa22Neg_C_Runs; CH2_Bkg = Fa22Neg_CH2_Runs; ET_Bkg = Fa22Neg_ET_Runs; F_Bkg = Fa22Neg_F_Runs;
	    }
	    else if( Period == "Fa22Pos" ){
		C_Bkg = Fa22Pos_C_Runs; CH2_Bkg = Fa22Pos_CH2_Runs; 
		// Scaling factors aren't applied to these runs, but it should cancel out in the calculation anyway...
		ET_Bkg = Fa22Neg_ET_Runs; 
		F_Bkg = Fa22Neg_F_Runs;
	    }
	    else if( Period == "Sp23Inb" ){
		C_Bkg = Sp23Inb_C_Runs; CH2_Bkg = Sp23Inb_CH2_Runs; ET_Bkg = Sp23Inb_ET_Runs; F_Bkg = Sp23Inb_F_Runs;
		CD2_Bkg = Sp23Inb_CD2_Runs;
	    }
	    else if( Period == "Solenoid_Scale" ){
		cout <<"Setting solenoid scaling factors...\n";
	    }
	    else{
		cout <<"ERROR: Invalid run period "<< Period <<". Valid options are: 'Su22' 'Fa22Neg' 'Fa22Pos' 'Sp23Inb'\n";
		cout <<"       No average values of X, Q2 were set!\n";
		return; // Kills the function
	    }

	    // This stores the average value of Q2 for each of the x, Q2 bins across ALL target types in the files
	    TProfile2D* Q2BinData = new TProfile2D({"Q2BinData","Q2 Bin Data; Q2 [GeV^{2}]; X", nQ2Bins, Q2_Bin_Bounds.data(), nXBins, X_Bin_Bounds.data()});
	    TProfile2D* XBinData  = new TProfile2D({"XBinData","X Bin Data; Q2 [GeV^{2}]; X", nQ2Bins, Q2_Bin_Bounds.data(), nXBins, X_Bin_Bounds.data()});

	    // Now loop over and create profiles for each of the target types in the vector
	    for(auto Target : Targets ){
	        // Just in case I forget about swapping ET and F with Empty and Foil...
	        if( Target == "ET" ) Target = "Empty";
	        else if( Target == "F" ) Target = "Foil";

	        // Set the runs to use
	        vector<int> Runs;
		// These first checks are for the special case wehre I'm calculating solenoid scaling factors
		if( Period == "Solenoid_Scale" && Target == "C" ){        Target = "C";   Runs = Fa22Pos_C_Runs; }
		else if( Period == "Solenoid_Scale" && Target == "ET" ){  Target = "C";   Runs = Fa22Neg_C_Runs; }
		else if( Period == "Solenoid_Scale" && Target == "CH2" ){ Target = "CH2"; Runs = Fa22Pos_CH2_Runs; }
		else if( Period == "Solenoid_Scale" && Target == "CD2" ){ Target = "CH2"; Runs = Fa22Neg_CH2_Runs; }
		// This is for other calculations
	        else if( Target == "C" ) Runs = C_Bkg;
	        else if( Target == "CH2" ) Runs = CH2_Bkg;
		else if( Target == "CD2" ) Runs = CD2_Bkg;
		else if( Target == "Empty" ) Runs = ET_Bkg;
		else if( Target == "Foil" ) Runs = F_Bkg;
		else{
		    cout <<"ERROR: Invalid target type "<< Target <<". This function is for background runs only. Check inputs.\n";
		    cout <<"       No average values of X, Q2 were set!\n";
		    return; // Kills the function
	        }

	    	for( int Run : Runs ){
		    string filePath = "/w/hallb-scshelf2102/clas12/holmberg/Spin_Structure/Latest_Skims/ROOT_Files/"+ Target +"_"+ to_string(Run) +"_Data.root";
		    ifstream test( filePath );
		    if( test.good() ){ // Proceed for a valid file path
		        TFile* Input = new TFile( filePath.c_str() );
		        auto q2In = (TProfile2D*)Input->Get("Q2BinData"); Q2BinData->Add( q2In );
		        auto xIn  = (TProfile2D*)Input->Get("XBinData");  XBinData->Add( xIn );
			//q2In->Delete();
			//xIn->Delete();
		        Input->Close();
		    }
		    else{ // Kill the operation if one of the runs is missing!!
		        cout <<"ERROR: Unable to locate ROOT file:\n";
		        cout <<"       "<< filePath << endl;
		        cout <<"       No average values of X, Q2 were set!\n";
		        return; // Kills the function
		    }
	    	}
	    } // End of target loop

	    // If everything was read in successfully, then set the average values for each target type:

	    for( int i=0; i<nQ2Bins; i++ ){
		for( int j=0; j<nXBins; j++ ){

		    int binIndex = XBinData->GetBin(i,j);
		    // Get the average x and Q2 values for this bin for each target type:
		    
		    // All target types
		    double All_X  = XBinData->GetBinContent(  binIndex );
		    double All_Q2 = Q2BinData->GetBinContent( binIndex );
		    this->SetAvgXQ2( "All", All_X, All_Q2 );


		} // End loop on X bins
	    } // End loop on Q2 bins

	    // Delete the profiles to prevent memory leaks...
	    XBinData->Delete();         Q2BinData->Delete();

	    cout << "Average X, Q2 set successfully using targets: ";
	    for(auto t : Targets) cout << t <<", ";
	    cout << endl;
	}

	// Helper function used in both "CalculateAvgXQ2()" functions above...
	void SetAvgXQ2(string Target, double xavg, double q2avg ){
	    for( int j=0; j<AllBins.size(); j++ )
	      for( int k=0; k<AllBins[j].size(); k++ )
		if( AllBins[j][k].IsInBin( q2avg, xavg ) ){ AllBins[j][k].SetAvgXQ2(Target, xavg, q2avg); break; }
	}

	// This reads in all the appropriate scaling factors for generating the CD2 and fall positive solenoid ET and F pseudo data
	void SetScalingFactors(){
	  ifstream finCD2("/w/hallb-scshelf2102/clas12/holmberg/Spin_Structure/Input_Text_Files/SuFa22_CH2_Scaling_Factors.txt");
	  ifstream finSol("/w/hallb-scshelf2102/clas12/holmberg/Spin_Structure/Input_Text_Files/Solenoid_Scaling_Factors.txt");
	  if( !finCD2.fail() && !finSol.fail() ){
	    string line;
	    while( getline(finCD2, line ) ){
	      double qmin, qmax, xmin, xmax, sf, sferr, qavg, xavg;
	      stringstream sin(line);
	      sin >> qmin >> qmax >> xmin >> xmax >> sf >> sferr >> qavg >> xavg;
	      for(size_t i=0; i<AllBins.size(); i++){
	 	for(size_t j=0; j<AllBins[i].size(); j++){
		    if( AllBins[i][j].IsInBin( qavg, xavg ) ){ AllBins[i][j].SetRatio( "CD2", sf, sferr ); }
		}
	      }
	    }
	    while( getline(finSol, line ) ){
	      double qmin, qmax, xmin, xmax, sf, sferr, qavg, xavg;
	      stringstream sin(line);
	      sin >> qmin >> qmax >> xmin >> xmax >> sf >> sferr >> qavg >> xavg;
	      for(size_t i=0; i<AllBins.size(); i++){
	 	for(size_t j=0; j<AllBins[i].size(); j++){
		    if( AllBins[i][j].IsInBin( qavg, xavg ) ){ AllBins[i][j].SetRatio( "Sol", sf, sferr ); }
		}
	      }
	    }

	  }
	  else cout <<"ERROR: Couldn't open all pseudo-data scaling factor files. Make sure files are in 'Input_Text_Files/'\n";
	  finCD2.close(); finSol.close();
	}


	void SlotCountsIntoBin(double qval, double xval, double hel_state){ // Used in Get_DF script
	    // Helicity state corresponds to either NP or NM, corrected for HWP and target polarization sign
	    for(size_t i=0; i<AllBins.size(); i++){
	 	for(size_t j=0; j<AllBins[i].size(); j++){
		    if( AllBins[i][j].IsInBin( qval, xval ) ){
			    					 // nm,  np,  n0
		      if( hel_state == 1 ) AllBins[i][j].AddCounts(1.0, 0.0, 0.0, "Input", xval, qval);
		      else if( hel_state == -1 ) AllBins[i][j].AddCounts(0.0, 1.0, 0.0, "Input", xval, qval);
		      else if( hel_state == 0  ) AllBins[i][j].AddCounts(0.0, 0.0, 1.0, "Input", xval, qval);
		    }
		}
	    }
	}
	void SlotChargeIntoBin(double fcm, double fcp, double fc0=0.0){
	    for(size_t i=0; i<AllBins.size(); i++){
	 	for(size_t j=0; j<AllBins[i].size(); j++){
		    AllBins[i][j].AddFCCharge( fcm, fcp, fc0, "Input" );
		}
	    }    
	}
	void NormalizeAllCounts(){
	    for(size_t i=0; i<AllBins.size(); i++){
	 	for(size_t j=0; j<AllBins[i].size(); j++){
		    AllBins[i][j].NormalizeAllCounts();
		}
	    }
	}
	// The "Period" variable refers to "Su22", "Fa22Neg", "Fa22Pos", "Sp23" for the appropriate scaling factors to use, if applicable
	void CalculateDF( string Period, bool useScaling, bool usePseudoData ){ // This also normalizes the counts when called
	    for(size_t i=0; i<AllBins.size(); i++){
	 	for(size_t j=0; j<AllBins[i].size(); j++){
		    AllBins[i][j].CalculateDF( Period, useScaling, usePseudoData );
		}
	    }
	}
	void CalculatePF( string Period, bool useScaling, bool usePseudoData ){ // This also normalizes the counts when called
	    for(size_t i=0; i<AllBins.size(); i++){
	 	for(size_t j=0; j<AllBins[i].size(); j++){
		    AllBins[i][j].CalculatePF( Period, useScaling, usePseudoData );
		}
	    }
	}
	void CalculateDFThruPF( string targetType, string Period, bool useScaling, bool usePseudoData ){ // Calculates the DF using the average of all the measured PF values
	    // Make sure that the PF is already calculated first:
	    for(size_t i=0; i<AllBins.size(); i++){
	 	for(size_t j=0; j<AllBins[i].size(); j++){
		    AllBins[i][j].CalculatePF( Period, useScaling, usePseudoData );
		}
	    }
	    // Now get the weighted average of the PF for the dataset
	    double totalPF = 0; double denominator = 0;
	    for(size_t i=0; i<AllBins.size(); i++){
	 	for(size_t j=0; j<AllBins[i].size(); j++){
		    double thisPF = 0; double thisPFerr = 0;
		    if( targetType == "NH3" ){
		        thisPF = AllBins[i][j].getPF_bath_NH3();
		        thisPFerr = AllBins[i][j].getErrPF_bath_NH3();
		        if( thisPF != 0.0 && thisPFerr != 0.0 ){ 
			    totalPF += thisPF / pow(thisPFerr,2) ;
			    denominator += 1.0 / pow(thisPFerr,2);
		        }
		    }
		    else if( targetType == "ND3" ){
		        thisPF = AllBins[i][j].getPF_bath_ND3();
		        thisPFerr = AllBins[i][j].getErrPF_bath_ND3();
		        if( thisPF != 0.0 && thisPFerr != 0.0 ){ 
			    totalPF += thisPF / pow(thisPFerr,2) ;
			    denominator += 1.0 / pow(thisPFerr,2);
			}
		    }
		}
	    }
	    // Now take the new weighted average PF and write it to each of the bins
	    if( denominator != 0.0 ){
	        double avgPF = totalPF / denominator;
		double avgPFerr = 1.0 / sqrt( denominator );
	        cout << "Average PF for "<< targetType<<" data set is: "<< avgPF <<" "<< avgPFerr << endl;
		PF_Avg = avgPF;
		PF_Avg_Err = avgPFerr;
	        // Now use this value to set the values of DF_Thru_PF with the new value
	        for(size_t i=0; i<AllBins.size(); i++){
	 	    for(size_t j=0; j<AllBins[i].size(); j++){
		        AllBins[i][j].CalculateDFThruPF( avgPF, 0, targetType, Period, useScaling, usePseudoData );
		    }
	        }
	    }
	    else{
		cout <<"ERROR: Zero passed into weighted avg. PF calculation; new PF values are not set.\n";
	    }
	}
	// Set the raw double-spin asymmetry for NH3
	void CalculateAllRaw(){
	    for(size_t i=0; i<AllBins.size(); i++){
		for(size_t j=0; j<AllBins[i].size(); j++){
		    AllBins[i][j].CalculateAllRaw();
		}
	    }
	}
	void SetDFs( double qval, double xval, double df, double dferr ){
	    for(size_t i=0; i<AllBins.size(); i++){
		for(size_t j=0; j<AllBins[i].size(); j++){
		    double Qmid = AllBins[i][j].getBinQ2(); double Xmid = AllBins[i][j].getBinX();
		    //cout << Epoch.getDF_FixdPF_NH3(qmid,xmid) <<" "<< Epoch.getErrDF_FixedPF_NH3(qmid,xmid) << endl;
		    if( AllBins[i][j].IsInBin( qval, xval ) ){
		        AllBins[i][j].SetDF_FixedPF_NH3( df, dferr ); 
		    }
		}
	    }
	}
	void SetSysErrDF(double syserr, string target, double qmid, double xmid){
	    for(size_t i=0; i<AllBins.size(); i++){
	 	for(size_t j=0; j<AllBins[i].size(); j++){
		    if( AllBins[i][j].IsInBin( qmid, xmid ) ){ 
			AllBins[i][j].SetSysErrDF(syserr, target);
		    }
		}
	    }
	}
	void SetA1( double a1, double a1err, double qmid, double xmid ){
	    for(size_t i=0; i<AllBins.size(); i++){
		for(size_t j=0; j<AllBins[i].size(); j++){
		    if( AllBins[i][j].IsInBin( qmid, xmid ) ) AllBins[i][j].SetA1( a1, a1err, "NH3" );
		}
	    }
	}

	// Accessors
	Bin getThisBin(double qmid, double xmid) const{
	    //Bin& thisBin;
	    for(size_t i=0; i<AllBins.size(); i++){
	 	for(size_t j=0; j<AllBins[i].size(); j++){
		    if( AllBins[i][j].IsInBin( qmid, xmid ) ){ 
			return AllBins[i][j];
		    }
		}
	    }
	    cout <<"ERROR: Couldn't find bin in range specified: Q2 = "<<qmid<<", X = "<<xmid<<endl;
	    // In the fail case, return an empty bin
	    return NullBin;
	}

	double getAFC(string target) const{
	    // The FC charges are common to all bins, so the bin they're grabbed from doesn't matter.
	    double FC_P = AllBins[0][0].getFC_P(target);
	    double FC_M = AllBins[0][0].getFC_M(target);
	    if( FC_P > 0 && FC_M > 0 ) return (FC_P - FC_M) / (FC_P + FC_M);
	    else{
		cout <<"ERROR: Some FC charges zero for "<< target <<". Check inputs\n";
		return -1;
	    }
	}

	double getPF_Avg(){ return PF_Avg; }
	double getPF_Avg_Err(){ return PF_Avg_Err; }
	double getBT_Pol() const{ return BT_Pol; }
	double getBT_Pol_Err() const{ return BT_Pol_Err; }

	// Returns the average x, Q2 value for the given target.
	// "All" returns the average for all target type runs together
	double getAvgQ2(string target, double qmid, double xmid) const{
	    for(size_t i=0; i<AllBins.size(); i++){
	 	for(size_t j=0; j<AllBins[i].size(); j++){
		    if( AllBins[i][j].IsInBin( qmid, xmid ) ){ 
			return AllBins[i][j].getAvgQ2( target );
		    }
		}
	    }
	    return 0.0; // Error case
	}
	double getAvgX(string target, double qmid, double xmid) const{
	    for(size_t i=0; i<AllBins.size(); i++){
	 	for(size_t j=0; j<AllBins[i].size(); j++){
		    if( AllBins[i][j].IsInBin( qmid, xmid ) ){ 
			return AllBins[i][j].getAvgX( target );
		    }
		}
	    }
	    return 0.0;
	}

	// Returns the scaling factors for generating pseudo data for the CD2, ET, and F for the appropriate data sets
	double getRatio(string target, double qmid, double xmid) const{
	    for(size_t i=0; i<AllBins.size(); i++){
	 	for(size_t j=0; j<AllBins[i].size(); j++){
		    if( AllBins[i][j].IsInBin( qmid, xmid ) ){ 
			return AllBins[i][j].getRatio( target );
		    }
		}
	    }
	    return 0.0;
	}
	double getRatioErr(string target, double qmid, double xmid) const{
	    for(size_t i=0; i<AllBins.size(); i++){
	 	for(size_t j=0; j<AllBins[i].size(); j++){
		    if( AllBins[i][j].IsInBin( qmid, xmid ) ){ 
			return AllBins[i][j].getRatioErr( target );
		    }
		}
	    }
	    return 0.0;
	}

/*
	double getAvgPF() const{
	    double total = 0.0; double counts = 0.0;
	    for(size_t i=0; i<AllBins.size(); i++){
	 	for(size_t j=0; j<AllBins[i].size(); j++){
		  if( AllBins[i][j].getPF_bath_NH3() != 0.0 ){
		    total += AllBins[i][j].getPF_bath_NH3();
		    counts++;
		  }
		}
	    }
	    if(counts != 0.0 ) return total/counts;
	    else return 0.0;
	}
*/
	double getNumeratorDFThruPF(double qmid, double xmid, string target) const{
	    for(size_t i=0; i<AllBins.size(); i++){
	 	for(size_t j=0; j<AllBins[i].size(); j++){
		    if( AllBins[i][j].IsInBin( qmid, xmid ) ){ 
			return AllBins[i][j].getNumeratorDFThruPF( target );
		    }
		}
	    }
	    cout <<"ERROR: Couldn't find bin in range specified: Q2 = "<<qmid<<", X = "<<xmid<<endl;
	    return 0.0;
	}
	double getNumeratorDFThruPFError(double qmid, double xmid, string target) const{
	    for(size_t i=0; i<AllBins.size(); i++){
	 	for(size_t j=0; j<AllBins[i].size(); j++){
		    if( AllBins[i][j].IsInBin( qmid, xmid ) ){ 
			return AllBins[i][j].getNumeratorDFThruPFError( target );
		    }
		}
	    }
	    cout <<"ERROR: Couldn't find bin in range specified: Q2 = "<<qmid<<", X = "<<xmid<<endl;
	    return 0.0;
	}
	double getNACountsDFThruPF(double qmid, double xmid, string target) const{
	    for(size_t i=0; i<AllBins.size(); i++){
	 	for(size_t j=0; j<AllBins[i].size(); j++){
		    if( AllBins[i][j].IsInBin( qmid, xmid ) ){ 
			return AllBins[i][j].getNACountsDFThruPF( target );
		    }
		}
	    }
	    cout <<"ERROR: Couldn't find bin in range specified: Q2 = "<<qmid<<", X = "<<xmid<<endl;
	    return 0.0;
	}
	double getNACountsDFThruPFError(double qmid, double xmid, string target) const{
	    for(size_t i=0; i<AllBins.size(); i++){
	 	for(size_t j=0; j<AllBins[i].size(); j++){
		    if( AllBins[i][j].IsInBin( qmid, xmid ) ){ 
			return AllBins[i][j].getNACountsDFThruPFError( target );
		    }
		}
	    }
	    cout <<"ERROR: Couldn't find bin in range specified: Q2 = "<<qmid<<", X = "<<xmid<<endl;
	    return 0.0;
	}

	// Accessors for the dilution factors
	double getDF_NH3(double qmid, double xmid) const{
	    for(size_t i=0; i<AllBins.size(); i++){
	 	for(size_t j=0; j<AllBins[i].size(); j++){
		    if( AllBins[i][j].IsInBin( qmid, xmid ) ){ 
			return AllBins[i][j].getDF_NH3();
		    }
		}
	    }
	    cout <<"ERROR: Couldn't find bin in range specified: Q2 = "<<qmid<<", X = "<<xmid<<endl;
	    return 0.0;
	}
	double getErrDF_NH3(double qmid, double xmid) const{
	    for(size_t i=0; i<AllBins.size(); i++){
	 	for(size_t j=0; j<AllBins[i].size(); j++){
		    if( AllBins[i][j].IsInBin( qmid, xmid ) ){ 
			return AllBins[i][j].getErrDF_NH3();
		    }
		}
	    }
	    cout <<"ERROR: Couldn't find bin in range specified: Q2 = "<<qmid<<", X = "<<xmid<<endl;
	    return 0.0;
	}
	double getSysErrDF_NH3(double qmid, double xmid) const{
	    for(size_t i=0; i<AllBins.size(); i++){
	 	for(size_t j=0; j<AllBins[i].size(); j++){
		    if( AllBins[i][j].IsInBin( qmid, xmid ) ){ 
			return AllBins[i][j].getSysErrDF_NH3();
		    }
		}
	    }
	    cout <<"ERROR: Couldn't find bin in range specified: Q2 = "<<qmid<<", X = "<<xmid<<endl;
	    return 0.0;
	}
	double getDF_ND3(double qmid, double xmid) const{
	    for(size_t i=0; i<AllBins.size(); i++){
	 	for(size_t j=0; j<AllBins[i].size(); j++){
		    if( AllBins[i][j].IsInBin( qmid, xmid ) ){ 
			return AllBins[i][j].getDF_ND3();
		    }
		}
	    }
	    cout <<"ERROR: Couldn't find bin in range specified: Q2 = "<<qmid<<", X = "<<xmid<<endl;
	    return 0.0;
	}
	double getErrDF_ND3(double qmid, double xmid) const{
	    for(size_t i=0; i<AllBins.size(); i++){
	 	for(size_t j=0; j<AllBins[i].size(); j++){
		    if( AllBins[i][j].IsInBin( qmid, xmid ) ){ 
			return AllBins[i][j].getErrDF_ND3();
		    }
		}
	    }
	    cout <<"ERROR: Couldn't find bin in range specified: Q2 = "<<qmid<<", X = "<<xmid<<endl;
	    return 0.0;
	}
	double getSysErrDF_ND3(double qmid, double xmid) const{
	    for(size_t i=0; i<AllBins.size(); i++){
	 	for(size_t j=0; j<AllBins[i].size(); j++){
		    if( AllBins[i][j].IsInBin( qmid, xmid ) ){ 
			return AllBins[i][j].getSysErrDF_ND3();
		    }
		}
	    }
	    cout <<"ERROR: Couldn't find bin in range specified: Q2 = "<<qmid<<", X = "<<xmid<<endl;
	    return 0.0;
	}
	double getTotalErrorDF(string target, double qmid, double xmid) const{
	    for(size_t i=0; i<AllBins.size(); i++){
	 	for(size_t j=0; j<AllBins[i].size(); j++){
		    if( AllBins[i][j].IsInBin( qmid, xmid ) ){ 
			return AllBins[i][j].getTotalErrorDF( target );
		    }
		}
	    }
	    cout <<"ERROR: Couldn't find bin in range specified: Q2 = "<<qmid<<", X = "<<xmid<<endl;
	    return 0.0;
	}

	double getDF_FixedPF_NH3(double qmid, double xmid) const{
	    for(size_t i=0; i<AllBins.size(); i++){
	 	for(size_t j=0; j<AllBins[i].size(); j++){
		    if( AllBins[i][j].IsInBin( qmid, xmid ) ){ 
			return AllBins[i][j].getDF_FixedPF_NH3();
		    }
		}
	    }
	    cout <<"ERROR: Couldn't find bin in range specified: Q2 = "<<qmid<<", X = "<<xmid<<endl;
	    return 0.0;
	}
	double getErrDF_FixedPF_NH3(double qmid, double xmid) const{
	    for(size_t i=0; i<AllBins.size(); i++){
	 	for(size_t j=0; j<AllBins[i].size(); j++){
		    if( AllBins[i][j].IsInBin( qmid, xmid ) ){ 
			return AllBins[i][j].getErrDF_FixedPF_NH3();
		    }
		}
	    }
	    cout <<"ERROR: Couldn't find bin in range specified: Q2 = "<<qmid<<", X = "<<xmid<<endl;
	    return 0.0;
	}
	double getDF_FixedPF_ND3(double qmid, double xmid) const{
	    for(size_t i=0; i<AllBins.size(); i++){
	 	for(size_t j=0; j<AllBins[i].size(); j++){
		    if( AllBins[i][j].IsInBin( qmid, xmid ) ){ 
			return AllBins[i][j].getDF_FixedPF_ND3();
		    }
		}
	    }
	    cout <<"ERROR: Couldn't find bin in range specified: Q2 = "<<qmid<<", X = "<<xmid<<endl;
	    return 0.0;
	}
	double getErrDF_FixedPF_ND3(double qmid, double xmid) const{
	    for(size_t i=0; i<AllBins.size(); i++){
	 	for(size_t j=0; j<AllBins[i].size(); j++){
		    if( AllBins[i][j].IsInBin( qmid, xmid ) ){ 
			return AllBins[i][j].getErrDF_FixedPF_ND3();
		    }
		}
	    }
	    cout <<"ERROR: Couldn't find bin in range specified: Q2 = "<<qmid<<", X = "<<xmid<<endl;
	    return 0.0;
	}

	// Accessors for the bath and cell packing fractions for NH3
	double getPF_bath_NH3(double qmid, double xmid) const{
	    for(size_t i=0; i<AllBins.size(); i++){
	 	for(size_t j=0; j<AllBins[i].size(); j++){
		    if( AllBins[i][j].IsInBin( qmid, xmid ) ){ 
			return AllBins[i][j].getPF_bath_NH3();
		    }
		}
	    }
	    cout <<"ERROR: Couldn't find bin in range specified: Q2 = "<<qmid<<", X = "<<xmid<<endl;
	    return 0.0;
	}
	double getErrPF_bath_NH3(double qmid, double xmid) const{
	    for(size_t i=0; i<AllBins.size(); i++){
	 	for(size_t j=0; j<AllBins[i].size(); j++){
		    if( AllBins[i][j].IsInBin( qmid, xmid ) ){ 
			return AllBins[i][j].getErrPF_bath_NH3();
		    }
		}
	    }
	    cout <<"ERROR: Couldn't find bin in range specified: Q2 = "<<qmid<<", X = "<<xmid<<endl;
	    return 0.0;
	}
	double getPF_cell_NH3(double qmid, double xmid) const{
	    for(size_t i=0; i<AllBins.size(); i++){
	 	for(size_t j=0; j<AllBins[i].size(); j++){
		    if( AllBins[i][j].IsInBin( qmid, xmid ) ){ 
			return AllBins[i][j].getPF_cell_NH3();
		    }
		}
	    }
	    cout <<"ERROR: Couldn't find bin in range specified: Q2 = "<<qmid<<", X = "<<xmid<<endl;
	    return 0.0;
	}
	double getErrPF_cell_NH3(double qmid, double xmid) const{
	    for(size_t i=0; i<AllBins.size(); i++){
	 	for(size_t j=0; j<AllBins[i].size(); j++){
		    if( AllBins[i][j].IsInBin( qmid, xmid ) ){ 
			return AllBins[i][j].getErrPF_cell_NH3();
		    }
		}
	    }
	    cout <<"ERROR: Couldn't find bin in range specified: Q2 = "<<qmid<<", X = "<<xmid<<endl;
	    return 0.0;
	}
	// Accessors for the bath and cell packing fractions for ND3
	double getPF_bath_ND3(double qmid, double xmid) const{
	    for(size_t i=0; i<AllBins.size(); i++){
	 	for(size_t j=0; j<AllBins[i].size(); j++){
		    if( AllBins[i][j].IsInBin( qmid, xmid ) ){ 
			return AllBins[i][j].getPF_bath_ND3();
		    }
		}
	    }
	    cout <<"ERROR: Couldn't find bin in range specified: Q2 = "<<qmid<<", X = "<<xmid<<endl;
	    return 0.0;
	}
	double getErrPF_bath_ND3(double qmid, double xmid) const{
	    for(size_t i=0; i<AllBins.size(); i++){
	 	for(size_t j=0; j<AllBins[i].size(); j++){
		    if( AllBins[i][j].IsInBin( qmid, xmid ) ){ 
			return AllBins[i][j].getErrPF_bath_ND3();
		    }
		}
	    }
	    cout <<"ERROR: Couldn't find bin in range specified: Q2 = "<<qmid<<", X = "<<xmid<<endl;
	    return 0.0;
	}
	double getPF_cell_ND3(double qmid, double xmid) const{
	    for(size_t i=0; i<AllBins.size(); i++){
	 	for(size_t j=0; j<AllBins[i].size(); j++){
		    if( AllBins[i][j].IsInBin( qmid, xmid ) ){ 
			return AllBins[i][j].getPF_cell_ND3();
		    }
		}
	    }
	    cout <<"ERROR: Couldn't find bin in range specified: Q2 = "<<qmid<<", X = "<<xmid<<endl;
	    return 0.0;
	}
	double getErrPF_cell_ND3(double qmid, double xmid) const{
	    for(size_t i=0; i<AllBins.size(); i++){
	 	for(size_t j=0; j<AllBins[i].size(); j++){
		    if( AllBins[i][j].IsInBin( qmid, xmid ) ){ 
			return AllBins[i][j].getErrPF_cell_ND3();
		    }
		}
	    }
	    cout <<"ERROR: Couldn't find bin in range specified: Q2 = "<<qmid<<", X = "<<xmid<<endl;
	    return 0.0;
	}
	double getAllTheory(double qmid, double xmid, string target) const{
	    for(size_t i=0; i<AllBins.size(); i++){
	 	for(size_t j=0; j<AllBins[i].size(); j++){
		    if( AllBins[i][j].IsInBin( qmid, xmid ) ){ 
			return AllBins[i][j].getAllTheory(target);
		    }
		}
	    }
	    cout <<"ERROR: Couldn't find bin in range specified: Q2 = "<<qmid<<", X = "<<xmid<<endl;
	    return 0.0;
	}
	double getAllRaw(double qmid, double xmid, string target) const{
	    for(size_t i=0; i<AllBins.size(); i++){
	 	for(size_t j=0; j<AllBins[i].size(); j++){
		    if( AllBins[i][j].IsInBin( qmid, xmid ) ){ 
			return AllBins[i][j].getAllRaw(target);
		    }
		}
	    }
	    cout <<"ERROR: Couldn't find bin in range specified: Q2 = "<<qmid<<", X = "<<xmid<<endl;
	    return 0.0;
	}
	double getAllRawErr(double qmid, double xmid, string target) const{
	    for(size_t i=0; i<AllBins.size(); i++){
	 	for(size_t j=0; j<AllBins[i].size(); j++){
		    if( AllBins[i][j].IsInBin( qmid, xmid ) ){ 
			return AllBins[i][j].getAllRawErr(target);
		    }
		}
	    }
	    cout <<"ERROR: Couldn't find bin in range specified: Q2 = "<<qmid<<", X = "<<xmid<<endl;
	    return 0.0;
	}
	double getAFC(double qmid, double xmid, string target) const{
	    for(size_t i=0; i<AllBins.size(); i++){
	 	for(size_t j=0; j<AllBins[i].size(); j++){
		    if( AllBins[i][j].IsInBin( qmid, xmid ) ){ 
			return AllBins[i][j].getAFC(target);
		    }
		}
	    }
	    cout <<"ERROR: Couldn't find bin in range specified: Q2 = "<<qmid<<", X = "<<xmid<<endl;
	    return 0.0;
	}
	double getAFC_Err(double qmid, double xmid, string target) const{
	    for(size_t i=0; i<AllBins.size(); i++){
	 	for(size_t j=0; j<AllBins[i].size(); j++){
		    if( AllBins[i][j].IsInBin( qmid, xmid ) ){ 
			return AllBins[i][j].getAFC_Err(target);
		    }
		}
	    }
	    cout <<"ERROR: Couldn't find bin in range specified: Q2 = "<<qmid<<", X = "<<xmid<<endl;
	    return 0.0;
	}

	double getAllRawNoFC(double qmid, double xmid, string target) const{
	    for(size_t i=0; i<AllBins.size(); i++){
	 	for(size_t j=0; j<AllBins[i].size(); j++){
		    if( AllBins[i][j].IsInBin( qmid, xmid ) ){ 
			return AllBins[i][j].getAllRawNoFC(target);
		    }
		}
	    }
	    cout <<"ERROR: Couldn't find bin in range specified: Q2 = "<<qmid<<", X = "<<xmid<<endl;
	    return 0.0;
	}
	double getAllRawNoFCErr(double qmid, double xmid, string target) const{
	    for(size_t i=0; i<AllBins.size(); i++){
	 	for(size_t j=0; j<AllBins[i].size(); j++){
		    if( AllBins[i][j].IsInBin( qmid, xmid ) ){ 
			return AllBins[i][j].getAllRawNoFCErr(target);
		    }
		}
	    }
	    cout <<"ERROR: Couldn't find bin in range specified: Q2 = "<<qmid<<", X = "<<xmid<<endl;
	    return 0.0;
	}

        // Getting the average PF for the NH3 targets
	double getAveragePF_bath_NH3() const{
	    double avgPF = 0.0; double counts = 0.0;
	    for(size_t i=0; i<AllBins.size(); i++){
	 	for(size_t j=0; j<AllBins[i].size(); j++){
		    if( AllBins[i][j].getPF_bath_NH3() > 0.0 ){
			avgPF += AllBins[i][j].getPF_bath_NH3();
			counts++;
		    }
		}
	    }
	    if( counts > 0.0 ) return avgPF / counts;
	    else return 0.0;
	}
	double getAveragePF_cell_NH3() const{
	    double avgPF = 0.0; double counts = 0.0;
	    for(size_t i=0; i<AllBins.size(); i++){
	 	for(size_t j=0; j<AllBins[i].size(); j++){
		    if( AllBins[i][j].getPF_cell_NH3() > 0.0 ){
			avgPF += AllBins[i][j].getPF_cell_NH3();
			counts++;
		    }
		}
	    }
	    if( counts > 0.0 ) return avgPF / counts;
	    else return 0.0;
	}
	// Getting the average PF for the ND3 targets
	double getAveragePF_bath_ND3() const{
	    double avgPF = 0.0; double counts = 0.0;
	    for(size_t i=0; i<AllBins.size(); i++){
	 	for(size_t j=0; j<AllBins[i].size(); j++){
		    if( AllBins[i][j].getPF_bath_ND3() > 0.0 ){
			avgPF += AllBins[i][j].getPF_bath_ND3();
			counts++;
		    }
		}
	    }
	    if( counts > 0.0 ) return avgPF / counts;
	    else return 0.0;
	}
	double getAveragePF_cell_ND3() const{
	    double avgPF = 0.0; double counts = 0.0;
	    for(size_t i=0; i<AllBins.size(); i++){
	 	for(size_t j=0; j<AllBins[i].size(); j++){
		    if( AllBins[i][j].getPF_cell_ND3() > 0.0 ){
			avgPF += AllBins[i][j].getPF_cell_ND3();
			counts++;
		    }
		}
	    }
	    if( counts > 0.0 ) return avgPF / counts;
	    else return 0.0;
	}
	double getCountRatio( double qmid, double xmid, string targ1, string targ2, string Period ){
	    for(size_t i=0; i<AllBins.size(); i++){
	 	for(size_t j=0; j<AllBins[i].size(); j++){
		    if( AllBins[i][j].IsInBin(qmid,xmid) ){
			return AllBins[i][j].getCountRatio(targ1, targ2, Period);
		    }
		}
	    }
	    return 0;
	}
	double getCountRatioErr( double qmid, double xmid, string targ1, string targ2, string Period ){
	    for(size_t i=0; i<AllBins.size(); i++){
	 	for(size_t j=0; j<AllBins[i].size(); j++){
		    if( AllBins[i][j].IsInBin(qmid,xmid) ){
			return AllBins[i][j].getCountRatioErr(targ1, targ2, Period);
		    }
		}
	    }
	    return 0;
	}
	double getAllPhys(double qmid, double xmid) const{
	    for(size_t i=0; i<AllBins.size(); i++){
		for(size_t j=0; j<AllBins[i].size(); j++){
		    if( AllBins[i][j].IsInBin( qmid, xmid ) ){
			return AllBins[i][j].getAllPhysNH3();
		    }
		}
	    }
	    return 0.0;
	}
	double getErrAllPhys(double qmid, double xmid) const{
	    for(size_t i=0; i<AllBins.size(); i++){
		for(size_t j=0; j<AllBins[i].size(); j++){
		    if( AllBins[i][j].IsInBin( qmid, xmid ) ){
			return AllBins[i][j].getAllPhysNH3Err();
		    }
		}
	    }
	    return 0.0;
	}

	double getA1(double qmid, double xmid) const{
	    for(size_t i=0; i<AllBins.size(); i++){
		for(size_t j=0; j<AllBins[i].size(); j++){
		    if( AllBins[i][j].IsInBin( qmid, xmid ) ){
			return AllBins[i][j].getA1NH3();
		    }
		}
	    }
	    return 0.0;
	}
	double getErrA1(double qmid, double xmid) const{
	    for(size_t i=0; i<AllBins.size(); i++){
		for(size_t j=0; j<AllBins[i].size(); j++){
		    if( AllBins[i][j].IsInBin( qmid, xmid ) ){
			return AllBins[i][j].getErrA1NH3();
		    }
		}
	    }
	    return 0.0;
	}
	double getA1Theory(double qmid, double xmid, string target) const{
	    for(size_t i=0; i<AllBins.size(); i++){
		for(size_t j=0; j<AllBins[i].size(); j++){
		    if( AllBins[i][j].IsInBin( qmid, xmid ) ){
			return AllBins[i][j].getA1Theory( target );
		    }
		}
	    }
	    return 0.0;
	}
	void SetAllPhys( double allphys, double allphyserr, double qmid, double xmid ){
	    for(size_t i=0; i<AllBins.size(); i++){
		for(size_t j=0; j<AllBins[i].size(); j++){
		    if( AllBins[i][j].IsInBin( qmid, xmid ) ) AllBins[i][j].SetAllPhys( allphys, allphyserr, "NH3" );
		}
	    }
	}
	void SetAllRaw( double allraw, double allrawerr, double qmid, double xmid ){
	    for(size_t i=0; i<AllBins.size(); i++){
		for(size_t j=0; j<AllBins[i].size(); j++){
		    if( AllBins[i][j].IsInBin( qmid, xmid ) ){
			 AllBins[i][j].SetAllRaw( allraw, "NH3" );
			 AllBins[i][j].SetAllRawErr( allrawerr, "NH3" );
		    }
		}
	    }
	}

	// This function sums up all the np and np counts across all bins and calculates
	// a total asymmetry for all bins. This isn't physical, but since it's supposed
	// to be a positive quantity, it can be used to check if the HWP for a given run
	// is correct. Returns a vector, with first entry the asymmetry, the second the
	// statistical error, the third the total number of normalized (NP+NM) counts,
	// the fourth the error on the counts, the fifth the raw asym. without FC normalization,
	// the sixth the FC charge asymmetry...
	vector<double> IntegratedDISAsymmetry(string target, bool normToFC=true) const{
	    double NP = 0; double NM = 0;
	    double FC_P= AllBins[0][0].getFC_P(target); // This value is the same across all bins 
	    double FC_M= AllBins[0][0].getFC_M(target);
	    vector<double> Asyms = {0,0,0,0,0,0};
	    for(size_t i=0; i<AllBins.size(); i++){
		for(size_t j=0; j<AllBins[i].size(); j++){
		    NP += AllBins[i][j].getNP(target);
		    NM += AllBins[i][j].getNM(target);
		}
	    }
	    double normNP = 0; double normNM = 0;
	    if( FC_P != 0 && FC_M != 0 && normToFC ){ normNP = NP/FC_P; normNM = NM/FC_M;}
	    else if( FC_P != 0 && FC_M != 0 ){ normNP = NP; normNM = NM;}

	    if( (normNP + normNM) != 0 && NP != 0 && NM != 0 ){
		Asyms[0] = (normNP - normNM) / (normNP + normNM);
		//Asyms[1] = 0.5 * sqrt( (NP + NM)/(NP*NM) );
		Asyms[1] = 1.0 / sqrt( NP + NM );
		Asyms[2] = normNP + normNM;
		Asyms[3] = sqrt( NP + NM )/( FC_P + FC_M );
		Asyms[4] = (NP - NM)/(NP + NM);
		//Asyms[5] = -1.0*(FC_P - FC_M)/(FC_P + FC_M); // post-hoc sign correction for mistake in SF_Data files
		Asyms[5] = (FC_P - FC_M)/(FC_P + FC_M); 
	    }
	    else cout << "ERROR: Unable to normalize counts in 'IntegratedDISAsymmetry() function.\n";
	    return Asyms;
	}

	// Returns the total number of counts for a specific target type across all helicity-latched counts.
	// For variable "N", specify "NP" for NP, "NM" for NM, and "Nt" for Nt.
	// Use "Unpol" to return all unpolarized counts.
	double ReturnCounts(string target, string N){

	    double totalCounts = 0.0;

	    if( target == "Unpol" ){
		vector<string> Targs = {"CH2", "CD2", "ET", "F", "C"};
		for(auto vec : AllBins){
		  for(Bin bin : vec){
		    for(string targ : Targs){
			if( N == "Nt" )      totalCounts += bin.getNt( targ );
			else if( N == "NP" ) totalCounts += bin.getNP( targ );
			else if( N == "NM" ) totalCounts += bin.getNM( targ );
		    }
		  }
		}
	    }
	    else{
		for(auto vec : AllBins){
	          for(Bin bin : vec){
		    if( N == "Nt" )      totalCounts += bin.getNt( target );
		    else if( N == "NP" ) totalCounts += bin.getNP( target );
		    else if( N == "NM" ) totalCounts += bin.getNM( target );
		  }
		}
	    }

	    return totalCounts;
	}

	// Using the "ReturnCounts()" function above, this function calculates the "r" factor NP / NM for all
	// runs in this DataSet object
	double CalculateR(){
	    
	    double allNP = this->ReturnCounts("Unpol","NP");
	    double allNM = this->ReturnCounts("Unpol","NM");

	    if( allNM != 0.0 ) return allNP / allNM;
	    else{
		cout << "ERROR: All NM values are zero for r = NP/NM (N^+ / N^-) calculation; check inputs.\n";
		return 0.0;
	    }
	}

	// Prints the total counts, FC charge, and normalized counts for the entire data set
	// for every available target type.
	void PrintTotalCounts() const{
	    double NA = 0; double FCA = 0; // NH3 counts and FC charge across all bins
	    double ND = 0; double FCD = 0; // ND3
	    double NCH= 0; double FCCH= 0; // CH2
	    double NCD= 0; double FCCD= 0; // CD2
	    double NC = 0; double FCC = 0; // Carbon
	    double NET= 0; double FCET= 0; // Empty target w/ LHe bath
	    double NF = 0; double FCF = 0; // Empty target w/o LHe (foils only)
	    for(size_t i=0; i<AllBins.size(); i++){
	 	for(size_t j=0; j<AllBins[i].size(); j++){
		    NA += AllBins[i][j].getNt("NH3"); FCA = AllBins[i][j].getFCt("NH3");
		    ND += AllBins[i][j].getNt("ND3"); FCD = AllBins[i][j].getFCt("ND3");
		    NCH+= AllBins[i][j].getNt("CH2"); FCCH= AllBins[i][j].getFCt("CH2");
		    NCD+= AllBins[i][j].getNt("CD2"); FCCD= AllBins[i][j].getFCt("CD2");
		    NC += AllBins[i][j].getNt("C");   FCC = AllBins[i][j].getFCt("C");
		    NET+= AllBins[i][j].getNt("ET");  FCET= AllBins[i][j].getFCt("ET");
		    NF += AllBins[i][j].getNt("F");   FCF = AllBins[i][j].getFCt("F");
		}
	    }
	    cout << "All Counts for this data set:\n";
	    if( FCA != 0 ) cout << "NH3 Counts: NA = "<< NA <<", FCA = "<< FCA <<", NormNA = " << NA/FCA << endl;
	    else cout << "NH3 Counts: NA = "<< NA <<", FCA = "<< FCA << endl;
	    if( FCD != 0 ) cout << "ND3 Counts: ND = "<< ND <<", FCD = "<< FCD <<", NormND = " << ND/FCD << endl;
	    else cout << "ND3 Counts: ND = "<< ND <<", FCD = "<< FCD << endl;
	    if( FCCH != 0 ) cout << "CH2 Counts: NCH = "<< NCH <<", FCCH = "<< FCCH <<", NormNCH = " << NCH/FCCH << endl;
	    else cout << "CH2 Counts: NCH = "<< NCH <<", FCCH = "<< FCCH << endl;
	    if( FCCD != 0 ) cout << "CD2 Counts: NCD = "<< NCD <<", FCCD = "<< FCCD <<", NormNCD = " << NCD/FCCD << endl;
	    else cout << "CD2 Counts: NCD = "<< NCD <<", FCCD = "<< FCCD << endl;
	    if( FCC != 0 ) cout << "Carbon Counts: NC = "<< NC <<", FCC = "<< FCC <<", NormNC = " << NC/FCC << endl;
	    else cout << "Carbon Counts: NC = "<< NC <<", FCC = "<< FCC << endl;
	    if( FCET != 0 ) cout << "ET Counts: NET = "<< NET <<", FCET = "<< FCET <<", NormNET = " << NET/FCET << endl;
	    else cout << "ET Counts: NET = "<< NET <<", FCET = "<< FCET << endl;
	    if( FCF != 0 ) cout << "Foil Counts: NF = "<< NF <<", FCF = "<< FCF <<", NormNF = " << NF/FCF << endl;
	    else cout << "Foil Counts: NF = "<< NF <<", FCF = "<< FCF << endl;

	}
	void Print( bool forcePrint=false ) const{
	    for(size_t i=0; i<AllBins.size(); i++){
	 	for(size_t j=0; j<AllBins[i].size(); j++){
		    AllBins[i][j].Print( forcePrint );
		}
	    }    
	}
	void PrintScalingFactors() const{
	    for(size_t i=0; i<AllBins.size(); i++){
	 	for(size_t j=0; j<AllBins[i].size(); j++){
		    cout <<"***** For Q2 = "<< AllBins[i][j].getBinQ2() <<", X = "<< AllBins[i][j].getBinX() <<" *****\n";
		    cout <<" --> ScaleCD2 = "<< AllBins[i][j].getRatio("CD2") <<" +- "<< AllBins[i][j].getRatioErr("CD2") << endl;
		    cout <<" --> ScaleSol = "<< AllBins[i][j].getRatio("Sol") <<" +- "<< AllBins[i][j].getRatioErr("Sol") << endl;
		    //cout <<" --> ScaleET = "<< AllBins[i][j].getRatio("ET") <<" +- "<< AllBins[i][j].getRatioErr("ET") << endl;
		    //cout <<" --> ScaleF = "<< AllBins[i][j].getRatio("F") <<" +- "<< AllBins[i][j].getRatioErr("F") << endl;
		    cout <<"********************\n";
		}
	    }    
	}

	void WriteToCSV(string outFilePath, string Period, string del=",", bool useScaling=true, bool usePseudoData=true ){ // Creates comma delimited file
	    string outPath;
	    if( del == ",") outPath = outFilePath + ".csv";
	    else outPath = outFilePath + ".txt";

	    ofstream fout(outPath);
	    if( !fout.fail() ){
		// Write the titles
		fout<<"Bin_Q2"<< del <<"Bin_X"<< del;
		fout<<"DF_NH3"<< del <<"Err_DF_NH3"<< del <<"PF_bath_NH3"<< del <<"ErrPF_bath_NH3"<< del <<"PF_cell_NH3"<< del <<"ErrPF_cell_NH3"<< del;
		fout<<"FC_total"<< del <<"N_NH3"<< del <<"N_C"<< del <<"N_CH2"<< del <<"N_ET"<< del <<"N_F"<< del <<"N_CD2"<< del <<
		      "FC_NH3"<< del <<"FC_C"<< del <<"FC_CH2"<< del <<"FC_ET"<< del <<"FC_F"<< del <<"FC_CD2\n";

		for(size_t i=0; i<AllBins.size(); i++){ // Loop over Q2 bins
	 	    for(size_t j=0; j<AllBins[i].size(); j++){ // Loop over x bins
		      if( AllBins[i][j].getDF_NH3() != 0.0 || AllBins[i][j].getDF_ND3() != 0.0 ){

			double qavg = AllBins[i][j].getAvgQ2("All"); double xavg = AllBins[i][j].getAvgX("All");
			double ETsf = 1;
			double ratioF = 1; double ratioET = 1; double ratioCD2 = 1; // Ratios that are used if corrections are necessary.

			if( useScaling ) ETsf = ET_Scale_Factor(  AllBins[i][j].getAvgX("ET"), AllBins[i][j].getAvgQ2("ET") );
			if( usePseudoData ){
			    if( Period == "Su22" || Period == "Fa22Neg" || Period == "Fa22Pos" ){
				ratioCD2 = AllBins[i][j].getRatio("CD2");
				//errCD = (1.0/fCD )*sqrt( ScaleCD2 *scale_nCD  + (pow(scale_nCD /ScaleCD2 ,2))*pow(ScaleCD2err ,2) );
			    }
			    // This is a separate check from above
			    if( Period == "Fa22Pos" ){
				ratioET = AllBins[i][j].getRatio("Sol");
				ratioF  = AllBins[i][j].getRatio("Sol");
				// Overwrite the ET and F errors to propagate errors correctly...
				//errET = (1.0/fET)*sqrt( ScaleET*scale_nET + (pow(scale_nET/ScaleET,2))*pow(ScaleETerr,2) );
				//errF  = (1.0/fF )*sqrt( ScaleF *scale_nF  + (pow(scale_nF /ScaleF ,2))*pow(ScaleFerr ,2) );
			    }
			}

			fout<<AllBins[i][j].getAvgQ2("All")<< del <<AllBins[i][j].getAvgX("All")<< del ;
			fout<<AllBins[i][j].getDF_NH3()<< del <<AllBins[i][j].getErrDF_NH3()<< del ;
			fout<<AllBins[i][j].getPF_bath_NH3()<< del <<AllBins[i][j].getErrPF_bath_NH3()<< del ;
			fout<<AllBins[i][j].getPF_cell_NH3()<< del <<AllBins[i][j].getErrPF_cell_NH3()<< del ;

			fout<<AllBins[i][j].getAllFCt()<< del ;
			fout<<AllBins[i][j].getNt("NH3")<< del <<AllBins[i][j].getNt("C")<< del <<AllBins[i][j].getNt("CH2")<< del <<AllBins[i][j].getNt("ET")*ETsf*ratioET<< del ;
			fout<<AllBins[i][j].getNt("F")*ratioF<< del <<AllBins[i][j].getNt("CD2")*ratioCD2<< del ;
			fout<<AllBins[i][j].getFCt("NH3")<< del <<AllBins[i][j].getFCt("C")<< del <<AllBins[i][j].getFCt("CH2")<< del <<AllBins[i][j].getFCt("ET")<< del ;
			fout<<AllBins[i][j].getFCt("F")<< del <<AllBins[i][j].getFCt("CD2")<<"\n";

		      }
		    }
	        }
	    }
	    fout.close();
	}

	void ReadDFfromTXT( string inFilePath ){ 
	    
	    ifstream fin( inFilePath );
	    if( !fin.fail() ){
		string line;
		getline(fin,line); // Throw away header row
		while( getline( fin, line ) ){
	
		    double Bin_Q2, Bin_X;
		    double DF_NH3, Err_DF_NH3, PF_bath_NH3, ErrPF_bath_NH3, PF_cell_NH3, ErrPF_cell_NH3;
		    double FC_total, N_NH3, N_C, N_CH2, N_ET, N_F, N_CD2;
		    double FC_NH3, FC_C, FC_CH2, FC_ET, FC_F, FC_CD2;

		    stringstream sin(line);
		    sin >> Bin_Q2 >> Bin_X >>
		           DF_NH3 >> Err_DF_NH3 >> PF_bath_NH3 >> ErrPF_bath_NH3 >> PF_cell_NH3 >> ErrPF_cell_NH3 >>
		           FC_total >> N_NH3 >> N_C >> N_CH2 >> N_ET >> N_F >> N_CD2 >>
		           FC_NH3 >> FC_C >> FC_CH2 >> FC_ET >> FC_F >> FC_CD2;

		    for( int i=0; i<AllBins.size(); i++ ){
			for( int j=0; j<AllBins[i].size(); j++ ){
			  if( AllBins[i][j].IsInBin( Bin_Q2, Bin_X ) ){

				//AllBins[i][j].SetAvgXQ2( "All", Bin_Q2, Bin_X );

				AllBins[i][j].SetDF_NH3( DF_NH3, Err_DF_NH3 );
				AllBins[i][j].SetPF( PF_bath_NH3, ErrPF_bath_NH3, "NH3", "Bath");
				AllBins[i][j].SetPF( PF_cell_NH3, ErrPF_cell_NH3, "NH3", "Cell");
				//AllBins[i][j].SetSysErrDF( syserrNH3, "NH3" ); 
				//AllBins[i][j].SetDF_FixedPF_NH3( dfpfNH3, dfpfNH3err );
			  }
			}
		    }
		} // End of while loop
	    } // End of if checking input file
	    fin.close();
	} // End of function


	// "dlmtr" is the delimiter, so it can be spaced with whitespace or a comma
	// Reads in data from a text file
	void ReadFromTXT(string inFilePath, bool useOnlyDFValues ){ 
	    
	    ifstream fin( inFilePath );
	    if( !fin.fail() ){
		string line;
		getline(fin,line); // Throw away header row
		while( getline( fin, line ) ){
	
		    double avgQ2, avgX, dfNH3, dfNH3err, dfpfNH3, dfpfNH3err, pfbathNH3, pfbathNH3err, pfcellNH3, pfcellNH3err;
		    double dfND3, dfND3err, dfpfND3, dfpfND3err, pfbathND3, pfbathND3err, pfcellND3, pfcellND3err;
		    double fctotal, nNH3, nC, nCH2, nET, nF, nCD2, fcNH3, fcC, fcCH2, fcET, fcF, fcCD2;
		    double syserrNH3 = 0; double syserrND3 = 0; // Initialized to zero for backwards compatibility if these aren't in OG file
		    double nmNH3, npNH3, fcmNH3, fcpNH3, nmND3, npND3, fcmND3, fcpND3;
		    stringstream sin(line);
		    sin >> avgQ2 >> avgX >> dfNH3 >> dfNH3err >> dfpfNH3 >> dfpfNH3err >> pfbathNH3 >> pfbathNH3err >> pfcellNH3 >> pfcellNH3err
		    >> dfND3 >> dfND3err >> dfpfND3 >> dfpfND3err >> pfbathND3 >> pfbathND3err >> pfcellND3 >> pfcellND3err
		    >> fctotal >> nNH3 >> nC >> nCH2 >> nET >> nF >> nCD2 >> fcNH3 >> fcC >> fcCH2 >> fcET >> fcF >> fcCD2 >> syserrNH3 >> syserrND3
		    >> nmNH3 >> npNH3 >> fcmNH3 >> fcpNH3
		    >> nmND3 >> npND3 >> fcmND3 >> fcpND3;

		    for( int i=0; i<AllBins.size(); i++ ){
			for( int j=0; j<AllBins[i].size(); j++ ){
			  if( AllBins[i][j].IsInBin( avgQ2, avgX ) ){

				AllBins[i][j].SetAvgXQ2( "All", avgX, avgQ2 );

				AllBins[i][j].SetDF_NH3( dfNH3, dfNH3err );
				AllBins[i][j].SetPF( pfbathNH3, pfbathNH3err, "NH3", "Bath");
				AllBins[i][j].SetPF( pfcellNH3, pfcellNH3err, "NH3", "Cell");
				AllBins[i][j].SetSysErrDF( syserrNH3, "NH3" ); 
				AllBins[i][j].SetDF_FixedPF_NH3( dfpfNH3, dfpfNH3err );

				AllBins[i][j].SetDF_ND3( dfND3, dfND3err );
				AllBins[i][j].SetPF( pfbathND3, pfbathND3err, "ND3", "Bath");
				AllBins[i][j].SetPF( pfcellND3, pfcellND3err, "ND3", "Cell");
				AllBins[i][j].SetSysErrDF( syserrND3, "ND3" );
				AllBins[i][j].SetDF_FixedPF_ND3( dfpfND3, dfpfND3err );

			     if( !useOnlyDFValues ){

				AllBins[i][j].AddCounts( nmNH3, npNH3, 0, "NH3", 0, 0 );
				AllBins[i][j].AddCounts( nC, 0, 0, "C", 0, 0 );
				AllBins[i][j].AddCounts( nCH2, 0, 0, "CH2", 0, 0 );
				AllBins[i][j].AddCounts( nET, 0, 0, "ET", 0, 0 );
				AllBins[i][j].AddCounts( nF, 0, 0, "F", 0, 0 );
				AllBins[i][j].AddCounts( nCD2, 0, 0, "CD2", 0, 0 );

				AllBins[i][j].AddFCCharge( fcmNH3, fcpNH3, 0, "NH3" );
				AllBins[i][j].AddFCCharge( fcC, 0, 0, "C" );
				AllBins[i][j].AddFCCharge( fcCH2, 0, 0, "CH2" );
				AllBins[i][j].AddFCCharge( fcET, 0, 0, "ET" );
				AllBins[i][j].AddFCCharge( fcF, 0, 0, "F" );
				AllBins[i][j].AddFCCharge( fcCD2, 0, 0, "CD2" );
			     }
			  }
			}
		    }
		}
	    }
	    else cout <<"ERROR: Failed to open file "<< inFilePath <<". Check inputs and try again.\n";

	}

	void WriteRawAsymsTXT(string outFilePath){ // Creates an output file for all the raw asymmetries in a tab delimited format
	    ofstream fout(outFilePath);
	    if( !fout.fail() ){
		// Write titles
		fout <<"Bin_Q2 Bin_X, AllrawNH3, ErrAllrawNH3, AllrawND3, ErrAllrawND3\n";
		for(size_t i=0; i<AllBins.size(); i++){
	 	    for(size_t j=0; j<AllBins[i].size(); j++){
		      if( AllBins[i][j].getDF_NH3() != 0.0 || AllBins[i][j].getDF_ND3() != 0.0){	
		      }
		    }
		}
	    }
	}
	void WriteToTXT(string outFilePath){ // Creates tab delimited file
	    ofstream fout(outFilePath);
	    if( !fout.fail() ){
		// Write the titles
		fout<<"Bin_Q2,Bin_X,DF_NH3,Err_DF_NH3,DF_ND3,Err_DF_ND3,PF_bath_NH3,ErrPF_bath_NH3,FC_total,N_NH3,N_C,N_CH2,N_ET,N_F,N_CD2,FC_NH3,FC_C,FC_CH2,FC_ET,FC_F,FC_CD2\n";
		for(size_t i=0; i<AllBins.size(); i++){
	 	    for(size_t j=0; j<AllBins[i].size(); j++){
		      if( AllBins[i][j].getDF_NH3() != 0.0 || AllBins[i][j].getDF_ND3() != 0.0){
			fout<<AllBins[i][j].getBinQ2()<<"	"<<AllBins[i][j].getBinX()<<"	"<<AllBins[i][j].getDF_NH3()<<"	"<<AllBins[i][j].getErrDF_NH3()<<"	";
			fout<<AllBins[i][j].getDF_ND3()<<"	"<<AllBins[i][j].getErrDF_ND3()<<"	"<<AllBins[i][j].getPF_bath_NH3()<<"	"<<AllBins[i][j].getErrPF_bath_NH3();
			fout<<AllBins[i][j].getAllFCt()<<"	";
			fout<<AllBins[i][j].getNt("NH3")<<"	"<<AllBins[i][j].getNt("C")<<"	"<<AllBins[i][j].getNt("CH2")<<"	"<<AllBins[i][j].getNt("ET")<<"	";
			fout<<AllBins[i][j].getNt("F")<<"	"<<AllBins[i][j].getNt("CD2")<<"	";
			fout<<AllBins[i][j].getFCt("NH3")<<"	"<<AllBins[i][j].getFCt("C")<<"	"<<AllBins[i][j].getFCt("CH2")<<"	"<<AllBins[i][j].getFCt("ET")<<"	";
			fout<<AllBins[i][j].getFCt("F")<<"	"<<AllBins[i][j].getFCt("CD2")<<"\n";
		      }
		    }
	        }
	    }
	    fout.close();
	}
	void WriteRawDFInputs(string outFilePath, vector<vector<TH1D*>>& KinBinsX, vector<vector<TH1D*>>& KinBinsQ2, bool append=false ){ // gives user option to append to file; default is false
	    //ofstream fout(outFilePath);
	    ofstream fout;
	    if ( !append ) fout.open(outFilePath);
	    else fout.open( outFilePath, ios::app );

	    if( !fout.fail() ){
		fout <<" Q2_Min   Q2_Max   X_Min   X_Max   NP_Counts   NM_Counts   FC_P_Charge   FC_M_Charge   N0_Counts   FC0_Charge\n";
		for(size_t i=0; i<AllBins.size(); i++){
		    for(size_t j=0; j<AllBins[i].size(); j++){
			fout << AllBins[i][j].ReturnBinContents() <<"	"<< KinBinsX[i][j]->GetMean() <<"	"<< KinBinsQ2[i][j]->GetMean() << endl;
		    }
		}
	    }
	    fout.close();
	}
	// Different version of the same function above, but looking at Q2 bins only
	void WriteRawDFInputs(string outFilePath, vector<TH1D*>& KinBinsQ2, vector<TH1D*>& KinBinsTheta, bool append=false ){ // gives user option to append to file; default is false
	    //ofstream fout(outFilePath);
	    ofstream fout;
	    if ( !append ) fout.open(outFilePath);
	    else fout.open( outFilePath, ios::app );

	    if( !fout.fail() ){
		fout <<" Q2_Min   Q2_Max   X_Min   X_Max   NP_Counts   NM_Counts   FC_P_Charge   FC_M_Charge   N0_Counts   FC0_Charge\n";
		for(size_t i=0; i<AllBins.size(); i++){
		    for(size_t j=0; j<AllBins[i].size(); j++){
			fout << AllBins[i][j].ReturnBinContents() <<"	"<< KinBinsQ2[i]->GetMean() <<"	"<< KinBinsTheta[i]->GetMean() << endl;
		    }
		}
	    }
	    fout.close();
	}

	void CalculateA1( int Run, string Target, RunPeriod& Period, bool useDISPbPt=true ){
	    // Make sure that everything has been properly calculated.
	    // Read in the dilution factors from the input epoch
	    // Once we have the raw asymmetries and target polarizations, calculate the A1 values for each bin
	    vector<double> pt_vals = {0,0};
	    if( useDISPbPt ) pt_vals = { abs(BT_Pol), abs(BT_Pol_Err) };
	    else  pt_vals = PT_Vals( Run, Period );

	    //cout <<"In A1 calculation, abs( PbPt ) = "<< pt_vals[0] <<" +/- "<< pt_vals[1] << endl;
	    for( size_t i=0; i<AllBins.size(); i++ ){ // Q2 bin loop
		for( size_t j=0; j<AllBins[i].size(); j++){ // x bin loop
		    AllBins[i][j].CalculateAllRaw();
		    AllBins[i][j].CalculateAllPhys( pt_vals[0], pt_vals[1] );
		    //cout << pt_vals[0] << endl;
		    Bin thisBin = AllBins[i][j];
		    //cout << thisBin.getDepolFactorNH3() <<" "<< thisBin.getAllPhysNH3() <<" "<< thisBin.getDF_NH3() << endl;
		    if( thisBin.getDepolFactorNH3() != 0.0 && thisBin.getAllPhysNH3() > 0.0 && thisBin.getDF_NH3() && pt_vals[0] != 0 ){
			//AllBins[i][j].CalculateAllPhys( pt_vals[0], pt_vals[1] );
		        double a1 = (thisBin.getAllPhysNH3()/(thisBin.getDepolFactorNH3())) - thisBin.getEtaNH3()*thisBin.getA2NH3();
		        double a1_err = thisBin.getAllPhysNH3Err()/thisBin.getDepolFactorNH3(); // temporary placeholder
			//cout << "A1 = "<< a1 <<" +- "<< a1_err << endl;
		        AllBins[i][j].SetA1( a1, a1_err, "NH3" );
		    }
		}
	    }
	}

	void MaxLikelihoodPbPt(int Run, RunPeriod& Period, bool normToFC=true){

	    double solPol = Period.getSolenoidScale( Run );
	    //double tarPol = Period.getTargetPolarization( Run ); // Only care about this for ammonia
	    // Tells whether or not to flip the HelP and HelN counts and FC charges
	    //bool flipStates = ( tarPol * solPol ) > 0;

	    // Get the target polarization sign
	    double NMR_Tpol = Period.getTargetPolarization( Run );

	    // Account for the sign of the target polarization.
	    // I essentially have to undo the sign corrections from reading in the data
	    /*
	    int corrFactor = 0;
	    if( NMR_Tpol < 0 ) corrFactor = -1;
	    else if( NMR_Tpol > 0 ) corrFactor = 1; // Allow for corrFactor to be zero for error tracing purposes */
	    int corrFactor = 1;
	    if( solPol > 0 ) corrFactor = -1;

	    double numPbPt = 0.0; double denomPbPt = 0.0; // Numerator and denominator terms used for PbPt calculation
	    double numPbPtErr = 0.0; double denomPbPtErr = 0.0; // Numerator and denominator terms used for error in PbPt calculation
							        // denomPbPtErr term needs to be squared at the end of the calculation!!!
	    
	    for(size_t i=0; i<AllBins.size(); i++){ // Q2 bin loop
		//double total = 0.0; double counts = 0.0;
		for(size_t j=0; j<AllBins[i].size(); j++){ // Loop over x bins for a given Q2 bin
		    double thisNP = AllBins[i][j].getNP("NH3"); // Raw counts NOT normalized to FC charge
		    double thisNM = AllBins[i][j].getNM("NH3"); // Same as above
		    double thisFC_P= AllBins[i][j].getFC_P("NH3");// The FC_P charge from the NH3 targets
		    double thisFC_M= AllBins[i][j].getFC_M("NH3");// Same as above for FC_M
		    //double thisNormNP = AllBins[i][j].getNormNP("NH3"); // NP counts normalized to FC_P
		    //double thisNormNM = AllBins[i][j].getNormNM("NH3"); // NM counts normalized to FC_M
		    //double thisDF = AllBins[i][j].getDF_FixedPF_NH3();
		    //double thisDFErr = AllBins[i][j].getErrDF_FixedPF_NH3();
		    double thisDF = AllBins[i][j].getDF_NH3();
		    double thisDFErr = AllBins[i][j].getErrDF_NH3();
		    double AllTheory = AllBins[i][j].getAllTheory("NH3");
		    if( thisDF > 0.0 && AllTheory != 0.0 && thisDFErr > 0.0 && thisFC_P > 0.0 && thisFC_M > 0.0 && thisNP > 0.0 && thisNM > 0.0 && AllBins[i][j].getBinX() < 0.75 ){
		    //if( thisDF > 0.0 && AllTheory != 0.0 && thisDFErr > 0.0 && thisNormNP > 0.0 && thisNormNM > 0.0 ){
			
			// This is just to temporarily get rid of the FC normalization...
			if( !normToFC ){ 
			    //thisFC_P = 1.0;
			    //thisFC_M = 1.0; 
			    numPbPt += corrFactor * thisDF * AllTheory * (1.0*thisNP - 1.0*thisNM);
			    denomPbPt += thisDF*thisDF * AllTheory*AllTheory * (1.0*thisNP + 1.0*thisNM);
			    numPbPtErr += thisDF*thisDF*AllTheory*AllTheory * (1.0*thisNP + 1.0*thisNM);
			    denomPbPtErr += thisDF*thisDF*AllTheory*AllTheory*(1.0*thisNP + 1.0*thisNM);
			}
			else{
			    /*
			    numPbPt += corrFactor * thisDF * AllTheory * ((thisNP/thisFC_P) - (thisNM/thisFC_M));
			    denomPbPt += thisDF*thisDF * AllTheory*AllTheory * ((thisNP/thisFC_P) + (thisNM/thisFC_M));
			    numPbPtErr += thisDF*thisDF*AllTheory*AllTheory * ((thisNP/thisFC_P) + (thisNM/thisFC_M));
			    denomPbPtErr += thisDF*thisDF*AllTheory*AllTheory*(thisNP + thisNM);
			    */
			    numPbPt += corrFactor * (2.0/(1.0 + (thisFC_P/thisFC_M) )) * (thisNP - (thisFC_P/thisFC_M)*thisNM )* thisDF *AllTheory;
   			    denomPbPt += thisDF*thisDF * AllTheory*AllTheory * (1.0*thisNP + 1.0*thisNM);
			    numPbPtErr += thisDF*thisDF*AllTheory*AllTheory * (1.0*thisNP + 1.0*thisNM);
			    denomPbPtErr += thisDF*thisDF*AllTheory*AllTheory*(1.0*thisNP + 1.0*thisNM);
			}
			/*
			numPbPt += corrFactor * thisDF * AllTheory * (thisNormNP - thisNormNM);
			denomPbPt += thisDF*thisDF * AllTheory*AllTheory * (thisNormNP + thisNormNM);
			numPbPtErr += thisDF*thisDF*AllTheory*AllTheory * ((thisNormNP/thisFC_P) + (thisNormNM/thisFC_M));
			denomPbPtErr += thisDF*thisDF*AllTheory*AllTheory*(thisNormNP + thisNormNM);
			*/
		    }
		}
	    } 
	    if( denomPbPt != 0.0 && denomPbPtErr != 0.0 ){
		BT_Pol = numPbPt / denomPbPt;
		//BT_Pol_Err = sqrt( numPbPtErr )/denomPbPtErr; // Pulled out of square root
		BT_Pol_Err = sqrt( 1.0 / denomPbPtErr );

	        //cout << "Beam-target polaization for run: PbPt = "<< BT_Pol <<" +- "<< BT_Pol_Err << endl;
		//cout << "FC charges are: Fp = "<< AllBins[0][0].getFC_P("NH3") <<", Fm = "<< AllBins[0][0].getFC_M("NH3") << endl;
		//BT_Pol_Err = BT_Pol*0.05;//sqrt( errTerm );
	    }
	    else
		cout << "ERROR: Invalid value of PbPt.\n";
	}

	

	// Destructor
	~DataSet() = default;	
};

// This function is used in the Faraday Cup charge asymmetry systematic error study. This function returns a TCanvas
// with three plots: A_ll normalized to the FC charge, A_ll unnormalized (raw counts), and half the difference between
// A_ll,norm and A_ll,raw. These are plotted for every x, Q2 bin.
// FIXME: Don't forget to correct for the sign of the target polarization and solenoid configuration!!!
TCanvas* FC_Asym_Systematic( DataSet& Data, string target, string period ){

    // Declare the TMultiGraph and TLegend that will hold the plot info for the normalized and NoFC asymmetries
    // as well as the difference plot
    TMultiGraph* mgAll = new TMultiGraph();
    TMultiGraph* mgAllNoFC = new TMultiGraph();
    TMultiGraph* mgDelta = new TMultiGraph();
    TLegend* leg = new TLegend(0.18,0.6,0.58,0.9);
    leg->SetNColumns(3);
    leg->SetHeader("Q^{2} Bins [GeV^{2}]","C");

    // The same FC charge is written to each bin for a target type, so it doens't
    // matter which bin it's drawn from
    double AFC = Data.getAFC( 2.3, 0.31, target);
    TF1* AFC_Line = new TF1("AFC_Line",to_string(AFC).c_str(), 0, 1);
    AFC_Line->SetLineColor(kBlack);

    // Loop over all x, Q2 bins to get the asymmetries for both normalized and raw:
    for(size_t i=0; i<Q2_Bin_Bounds.size()-1; i++){

	double qmid = (Q2_Bin_Bounds[i] + Q2_Bin_Bounds[i+1]) / 2.0;

	// Vectors that are used to hold the plots for each bin
	vector<double> xVals, All, AllErr, AllNoFC, AllNoFCErr, Delta;

	for(size_t j=0; j<X_Bin_Bounds.size()-1; j++){

	    double xmid = (X_Bin_Bounds[j] + X_Bin_Bounds[j+1]) / 2.0;	    
	    double Allval     = Data.getAllRaw( qmid, xmid, target );
	    double Allvalerr  = Data.getAllRawErr( qmid, xmid, target );
	    double Allnofc    = Data.getAllRawNoFC( qmid, xmid, target );
	    double Allnofcerr = Data.getAllRawNoFCErr( qmid, xmid, target );

	    //cout << qmid <<" "<< xmid <<" "<< Allval <<" "<< Allnofc << endl;

	    if( Allval != 0 && Allnofc != 0 ){
		xVals.push_back( xmid ); //FIXME: Add in the average x value for this x, Q2 bin instead of the bin midpoint
		All.push_back( Allval ); AllErr.push_back( Allvalerr );
		AllNoFC.push_back( Allnofc ); AllNoFCErr.push_back( Allnofcerr );
		Delta.push_back( (Allnofc - Allval) );
	    }
	}
	// Now make a plot if there's enough data...
	if( All.size() > 0 && AllNoFC.size() > 0 ){
	    // Normalized All plot
	    TGraphErrors* gr = new TGraphErrors( All.size(), xVals.data(), All.data(), nullptr, AllErr.data() );
	    gr->SetMarkerStyle( kFullCircle );
	    gr->SetMarkerColor( Palette[i] );
	    mgAll->Add( gr, "p" );
	    stringstream s; s <<"Q^{2} = "<< setprecision(4) << qmid;
	    leg->AddEntry( gr, s.str().c_str(), "p" );

	    // Raw All plot (no FC charge normalization)
	    TGraphErrors* grNoFC = new TGraphErrors( AllNoFC.size(), xVals.data(), AllNoFC.data(), nullptr, AllNoFCErr.data() );
	    grNoFC->SetMarkerStyle( kFullCircle );
	    grNoFC->SetMarkerColor( Palette[i] );
	    mgAllNoFC->Add( grNoFC, "p" );

	    // Half the difference between the two previous plots
	    TGraphErrors* grDelta = new TGraphErrors( Delta.size(), xVals.data(), Delta.data(), nullptr, nullptr );
	    grDelta->SetMarkerStyle( kFullCircle );
	    grDelta->SetMarkerColor( Palette[i] );
	    mgDelta->Add( grDelta, "p" );

	}

    } // End of Q2 bins loop

    // Set the titles for the TMultiGraphs
    mgAll->SetTitle( string("A_{||,norm} for "+ period +"; X; A_{||,norm}").c_str() );
    mgAllNoFC->SetTitle( string("A_{||,raw} for "+ period +"; X; A_{||,raw}").c_str() );
    mgDelta->SetTitle( string("A_{||,raw}-A_{||,norm} #approx A_{FC} for "+ period +"; X; #approx A_{FC}").c_str() );

    // Now make the plots...
    string canName = "FC_Asym_Systematic_"+period;
    TCanvas* c = new TCanvas( canName.c_str(), canName.c_str(), 2200, 600 );
    c->Divide(3,1);
    c->cd(1);
    gPad->SetLeftMargin(0.18); gPad->SetRightMargin(0.02);
    mgAll->Draw("ap"); leg->Draw("same");

    c->cd(2);
    gPad->SetLeftMargin(0.18); gPad->SetRightMargin(0.02);
    mgAllNoFC->Draw("ap"); leg->Draw("same");

    c->cd(3);
    gPad->SetLeftMargin(0.18); gPad->SetRightMargin(0.02);
    mgDelta->Draw("ap"); AFC_Line->Draw("same");

    return c;
}

// Reads in A1 data from previous experiments; used in the Make_A1_NH3_Plots function
// Data in the following form: Compass
vector<TGraphErrors*> Old_A1_Data(){
 
    vector<int> palette = { 632, 800, 416, 600, 840, 880, 900, 616, 820, 432, 920, 860 };
    int pint = 0;
   
    vector<vector<double>> xVals(7);
    vector<vector<double>> A1data(7);
    vector<vector<double>> A1err(7);
    vector<vector<double>> zeros(7);
    vector<string> Titles = {"Compass","E143","E155","EMC","Hermes","SMC","EG1b"};
    vector<int> Styles = {3,4,8,21,26,30,34};
    string infile = THIS_DIR + "A1p_DIS_Pushpa.txt";
    ifstream fin(infile);
    string line; getline(fin,line); // Throw out header row
    while( getline( fin, line ) ){
	double Expt, x, q2, w2, a1data, a1stat, a1sys, a1err;
	stringstream sin(line);
	sin >> Expt >> x >> q2 >> w2 >> a1data >> a1stat >> a1sys >> a1err;
	if( Expt != 9 && Expt < 11 ){
	    int index; 
	    if( Expt <= 3 ) index = 0;
	    else if( Expt == 10 ) index = 6;
	    else index = Expt - 3;
	    xVals[index].push_back(x);
	    A1data[index].push_back(a1data);
	    A1err[index].push_back(a1err);
	    zeros[index].push_back(0);
	}
    }
    fin.close();
    // Now make the plots for each data set
    vector<TGraphErrors*> OldData;
    TCanvas* DataPlot = new TCanvas("DataPlot","DataPlot",800,600);
    TMultiGraph* mg = new TMultiGraph();
    TLegend* leg = new TLegend(0.1,0.7,0.5,0.9); leg->SetNColumns(2);
    for(size_t i=0; i<A1data.size(); i++){
	TGraphErrors* gr = new TGraphErrors(xVals[i].size(),xVals[i].data(),A1data[i].data(),zeros[i].data(),A1err[i].data());
	gr->SetMarkerColor(kGray);
	gr->SetLineColor(kGray);
	//gr->SetMarkerColor( palette[pint] );
	//gr->SetLineColor( palette[pint] );
	gr->SetMarkerStyle(Styles[i]);
	leg->AddEntry(gr,Titles[i].c_str(),"p");
	mg->Add(gr,"p");
	OldData.push_back(gr);
	pint++;
    }
    mg->SetTitle("Old A_{1,p} Data; X; A_{1,p}");
    mg->GetYaxis()->SetRangeUser(0,1.2); 
    DataPlot->cd(); mg->Draw("ap"); leg->Draw("same");
    DataPlot->Print("Old_A1p_Data.pdf");
    return OldData;
}

// Makes a plot of both All,phys and A1 for the proton
//vector<TCanvas*> Make_A1_NH3_Plots( vector<PbPt>& AllPhysVals, DataSet& Epoch, string sector, string color ){
//vector<TCanvas*> Make_A1_NH3_Plots( vector<PbPt>& AllPhysVals, string sector, string color ){
vector<TCanvas*> Make_A1_NH3_Plots( vector<DataSet>& AllPhysVals, string sector, string color, vector<int>& Runs, string filePath, string target){

//vector<TGraphErrors*> Make_A1_NH3_Plots( vector<PbPt>& AllPhysVals, DataSet& Epoch, string sector, string color ){

    //PbPt AllA1Vals; AllA1Vals.SetTarget("NH3");
    // This DataSet object holds the weighted average values of all A1 values from the given runs
    DataSet AllA1Vals;

    for( size_t i=0; i<Q2_Bin_Bounds.size()-1; i++ ){ // q2 bin loop
        for( size_t j=0; j<X_Bin_Bounds.size()-1; j++ ){ // x bin loop
	    // These "Num" and "Denom" terms are used to calculate the weighted averages of the physical and raw asymmetreis as well
	    // as the NH3 values of A1
	    double thisBinNum = 0.0;
	    double thisBinDenom = 0.0;
	    double thisBinRawNum = 0.0;
	    double thisBinRawDenom = 0.0;
	    double thisBinA1Num = 0.0;
	    double thisBinA1Denom = 0.0;

	    double qmid = (Q2_Bin_Bounds[i]+Q2_Bin_Bounds[i+1])/2.0;
	    double xmid = (X_Bin_Bounds[j]+X_Bin_Bounds[j+1])/2.0;
	    //cout << "Looking at bin "<< qmid <<" "<< xmid <<":\n";
	    // For this particular x, Q2 bin, loop over all DataSet objects (one per run) and take the weighted average of all values
            for( size_t k=0; k<AllPhysVals.size(); k++){
	        // Loop through each bin and get the average AllPhys
		double thisBinAllPhys = AllPhysVals[k].getAllPhys( qmid, xmid );
		double thisBinAllPhysErr = AllPhysVals[k].getErrAllPhys( qmid, xmid );
		double thisBinAllRaw = AllPhysVals[k].getAllRaw( qmid, xmid, "NH3" );
		double thisBinAllRawErr = AllPhysVals[k].getAllRawErr( qmid, xmid, "NH3" );
		double thisBinA1 = AllPhysVals[k].getA1( qmid, xmid );
		double thisBinA1Err = AllPhysVals[k].getErrA1( qmid, xmid );

	        //cout << "Checking AllPhys for bin "<< qmid <<" "<< xmid <<" "<< thisBinAllPhys <<" +- "<< thisBinAllPhysErr << endl;
		// If one of these quantities doesn't have a valid value, toss it out...
		if( thisBinAllPhys > 0 && thisBinAllPhysErr > 0 && thisBinAllRaw > 0 && thisBinAllRawErr > 0 && thisBinA1 > 0 && thisBinA1Err > 0 ){
		    thisBinNum += thisBinAllPhys / pow(thisBinAllPhysErr,2);
		    thisBinDenom += 1.0 / pow(thisBinAllPhysErr,2);
		    thisBinRawNum += thisBinAllRaw / pow(thisBinAllRawErr,2);
		    thisBinRawDenom += 1.0 / pow(thisBinAllRawErr,2);
		    thisBinA1Num += thisBinA1 / pow( thisBinA1Err, 2 );
		    thisBinA1Denom += 1.0 / pow( thisBinA1Err, 2 );
		    if( i==0 && j==0 ){
			//cout << "PbPt for run "<< AllPhysVals[k].getRunNumber() <<": ";
			//cout << AllPhysVals[k].getBT_Pol() <<" +/- "<< AllPhysVals[k].getBT_Pol_Err() << endl;
		    }
		}
	    }
	    // After the loop, take the weighted average, if enough data is available...
	    if( thisBinDenom > 0.0 && thisBinRawDenom > 0.0 && thisBinA1Denom > 0.0 ){
	        double AvgAllPhys = thisBinNum / thisBinDenom;
	        double AvgAllPhysErr = 1.0 / sqrt( thisBinDenom );
	        double AvgAllRaw = thisBinRawNum / thisBinRawDenom;
	        double AvgAllRawErr = 1.0 / sqrt( thisBinRawDenom );
		double AvgA1NH3 = thisBinA1Num / thisBinA1Denom;
		double AvgA1NH3Err = 1.0 / sqrt( thisBinA1Denom );
	        //cout << "Entering AllRaw for bin "<< qmid <<" "<< xmid <<" "<< AvgAllRaw <<" +- "<< AvgAllRawErr << endl;
		// Write the values to the DataSet object
	        AllA1Vals.SetAllRaw( AvgAllRaw, AvgAllRawErr, qmid, xmid );
	        AllA1Vals.SetAllPhys( AvgAllPhys, AvgAllPhysErr, qmid, xmid );
		AllA1Vals.SetA1( AvgA1NH3, AvgA1NH3Err, qmid, xmid ); 

	        //AllA1Vals.CalculateA1();
	    }
	}
    }
    // AllA1Vals.CalculateA1();

    // Now make the plot
    vector<int> palette = { 1, 632, 800, 400, 416, 600, 840, 880, 900, 616, 820, 432, 920, 860 };
    int pint = 0;

    // Set the weighted average values of X, Q2 for each of the bins, which is relevant for making the A1 values appear
    // in the correct spot along the x axis
    // FIXME: Replace this with the "CalculateAvgXQ2()" function
    //AllA1Vals.SetAvgXQ2( Runs, filePath, target );

    // Hold the plot objects
    TMultiGraph* mgA1_NH3       = new TMultiGraph(); // Holds the A1 values for NH3
    TMultiGraph* mgAllPhys_NH3  = new TMultiGraph(); // Holds the physical asymmetries for NH3
    TMultiGraph* mgAllRaw_NH3   = new TMultiGraph(); // Holds the raw asymmetries for NH3
    TMultiGraph* mgA1_NH3_Q2Avg = new TMultiGraph(); // Holds the A1 NH3 values averaged over all Q2 bins ( i.e., plots A1(x) )

    TLegend* leg = new TLegend(0.1,0.7,0.5,0.9);
    leg->SetNColumns(3);
    leg->SetHeader("Q^{2} Bins (GeV^{2})","C");

    // Add this historical data to the plot
    vector<TGraphErrors*> OldA1Plots = Old_A1_Data(); // Historical A1 data
    for(size_t i=0; i<OldA1Plots.size(); i++){ 
	mgA1_NH3->Add(OldA1Plots[i],"p");
	mgA1_NH3_Q2Avg->Add(OldA1Plots[i],"p");
    }

    // Used for plotting the systematic error band
    int nSteps = X_Bin_Bounds.size()-1;

    for(size_t i=0; i<Q2_Bin_Bounds.size()-1; i++){
	double qmid = (Q2_Bin_Bounds[i]+Q2_Bin_Bounds[i+1])/2.0;
	vector<double> xbins, a1nh3vals, zeros, a1nh3errs, allphysnh3vals, allphysnh3errs, allrawnh3vals, allrawnh3errs;
	for(size_t j=0; j<X_Bin_Bounds.size()-1; j++){
	    double xmid = (X_Bin_Bounds[j]+X_Bin_Bounds[j+1])/2.0;
	    double xavg = AllA1Vals.getAvgX( target, qmid, xmid ); // Statistics-weighted avg. X across all the NH3 or ND3 runs used
	    double thisA1NH3 = AllA1Vals.getA1( qmid, xmid );
	    double thisErrA1NH3 = AllA1Vals.getErrA1( qmid, xmid);
	    double thisAllPhysNH3 = AllA1Vals.getAllPhys( qmid, xmid );
	    double thisAllPhysErrNH3 = AllA1Vals.getErrAllPhys( qmid, xmid );
	    double thisAllRawNH3 = AllA1Vals.getAllRaw( qmid, xmid, "NH3" );
	    double thisAllRawErrNH3 = AllA1Vals.getAllRawErr( qmid, xmid, "NH3" );
	    if( thisA1NH3 != 0.0 && thisAllPhysNH3 != 0.0 /*&& thisAllPhysNH3 < 0.5*/ ){
		xbins.push_back(xavg); 
		zeros.push_back(0.0);
		a1nh3vals.push_back( thisA1NH3 );
		a1nh3errs.push_back( thisErrA1NH3 );
		allphysnh3vals.push_back( thisAllPhysNH3 );
		allphysnh3errs.push_back( thisAllPhysErrNH3 );
		allrawnh3vals.push_back( thisAllRawNH3 );
		allrawnh3errs.push_back( thisAllRawErrNH3 );
		//cout << thisA1NH3 <<" +- "<< thisErrA1NH3 << endl;
		//cout << thisAllPhysNH3 <<" +- "<< thisAllPhysErrNH3 << endl;
	    }
	}
	if( xbins.size() > 0 ){
	    TGraphErrors* gr = new TGraphErrors( xbins.size(), xbins.data(), a1nh3vals.data(), zeros.data(), a1nh3errs.data() );
	    gr->SetMarkerColor( palette[pint] );
	    TGraphErrors* gr2= new TGraphErrors( xbins.size(), xbins.data(), allphysnh3vals.data(), zeros.data(), allphysnh3errs.data() );
	    gr2->SetMarkerColor( palette[pint] );
	    TGraphErrors* grRaw= new TGraphErrors( xbins.size(), xbins.data(), allrawnh3vals.data(), zeros.data(), allrawnh3errs.data() );
	    grRaw->SetMarkerColor( palette[pint] );

	    if(pint < 4 ){ 
		gr->SetMarkerStyle(kFullCircle); 
		gr2->SetMarkerStyle(kFullCircle);  
		grRaw->SetMarkerStyle(kFullCircle);
	    }
	    else if(pint >= 4 && pint < 8){ 
		gr->SetMarkerStyle(kFullTriangleUp); 
		gr2->SetMarkerStyle(kFullTriangleUp); 
		grRaw->SetMarkerStyle(kFullTriangleUp);
	    }
	    else if(pint >= 8){ 
		gr->SetMarkerStyle(kFullSquare); 
		gr2->SetMarkerStyle(kFullSquare); 
		grRaw->SetMarkerStyle(kFullSquare);
	    }
	    mgA1_NH3->Add( gr, "p" );
	    mgAllPhys_NH3->Add( gr2, "p" );
	    mgAllRaw_NH3->Add( grRaw, "p" );
	    stringstream legTitle; legTitle << "Q^{2}=" << fixed << setprecision(3) << qmid;
	    leg->AddEntry( gr, legTitle.str().c_str(), "p");
	    pint++;
	}
    }
    // Now make the plot for the fixed Q2. For a fixed value in X, take the average over all Q2 bins for that value of X.
    // This doesn't take into account the statistics in each bin, it's just the average.
    // While this isn't accurate, I report the average X of each bin to be the average of the weighted averages for each
    // kinematic bin... I know that's not correct and makes like no sense but it's easier and just for show lol.
    vector<double> avgA1vals, Xvals, Zeros, avgA1valsErr, avgA1theory, sysErrorBand;
    for(size_t i=0; i<X_Bin_Bounds.size()-1; i++){
	double xmid = (X_Bin_Bounds[i]+X_Bin_Bounds[i+1])/2.0;
	double thisXBinA1num = 0.0; double thisXBinA1denom = 0.0;
	double sumA1Theory = 0.0; double counts = 0.0;
	double xAvgSum = 0.0;
	for(size_t j=0; j<Q2_Bin_Bounds.size()-1; j++){
	    double qmid = (Q2_Bin_Bounds[j]+Q2_Bin_Bounds[j+1])/2.0;
	    double xAvg = AllA1Vals.getAvgX( target, qmid, xmid );
	    double thisA1NH3 = AllA1Vals.getA1(qmid, xmid);
	    double thisErrA1NH3 = AllA1Vals.getErrA1( qmid, xmid );
	    double thisA1TheoryNH3 = AllA1Vals.getA1Theory( qmid, xmid, "NH3" );
	    if( thisA1NH3 > 0 && thisErrA1NH3 > 0 && thisErrA1NH3 < 0.5 ){
		thisXBinA1num += thisA1NH3 / pow(thisErrA1NH3,2);
		thisXBinA1denom += 1.0 / pow(thisErrA1NH3,2);
		sumA1Theory += thisA1TheoryNH3; counts++;
		xAvgSum += xAvg;
	    }
	}
	if( thisXBinA1denom != 0 && counts != 0 ){
	    double a1Val = thisXBinA1num / thisXBinA1denom;
	    double truexAvg = xAvgSum / counts;
	    avgA1vals.push_back( a1Val );
	    avgA1valsErr.push_back( 1.0/sqrt(thisXBinA1denom) );
	    //Xvals.push_back( xmid ); 
	    Xvals.push_back( truexAvg );
	    Zeros.push_back(0);
	    avgA1theory.push_back( sumA1Theory / counts );
	    // Add to the systematic error band, which ranges from 11% at low x to 6% at high x
	    sysErrorBand.push_back( ( 0.11-(i*0.05/(1.0*nSteps)) )*a1Val );
	}
    }
    vector<TGraphErrors*> A1Plots;
    TGraphErrors* gr3 = new TGraphErrors(Xvals.size(), Xvals.data(), avgA1vals.data(), Zeros.data(), avgA1valsErr.data() );
    // This graph holds the systematic errors
    TGraphErrors* sysErr = new TGraphErrors(Xvals.size(), Xvals.data(), avgA1vals.data(), Zeros.data(), sysErrorBand.data() );
    sysErr->SetMarkerStyle(1); // make an invisible point
    sysErr->SetFillColor(kOrange);

    if( color == "Blue" ){
        gr3->SetMarkerColor( kBlack ); gr3->SetMarkerStyle(kFullCircle);
	gr3->SetMarkerSize(1);
        TGraphErrors* gr4 = new TGraphErrors(Xvals.size(), Xvals.data(), avgA1theory.data(), Zeros.data(), Zeros.data() );
        gr4->SetLineColor( kRed ); gr4->SetLineWidth(1);
	//mgA1_NH3_Q2Avg->Add( sysErr, "A3" );
        mgA1_NH3_Q2Avg->Add( gr3, "p" );
        mgA1_NH3_Q2Avg->Add( gr4, "l" );
	mgA1_NH3->Add( gr4, "l" );
	A1Plots.push_back( gr3 ); A1Plots.push_back( gr4 );
    }
    else{
        gr3->SetMarkerColor( kBlack ); gr3->SetMarkerStyle(kFullSquare);
        mgA1_NH3_Q2Avg->Add( gr3, "p" );
	A1Plots.push_back( gr3 );
    }
    //mgA1_NH3_Q2Avg->Add( gr3, "p" );
    //mgA1_NH3_Q2Avg->Add( gr4, "l" );


    string title = "A_{1,p}(X,Q^{2}) for "+sector+"; X; A_{1,p}";
    mgA1_NH3->SetTitle(title.c_str());
    mgA1_NH3->GetXaxis()->SetLimits(0,0.8); //mgA1_NH3->GetYaxis()->SetRangeUser(0.0,0.4);
    mgA1_NH3->GetYaxis()->SetRangeUser(0,1.2);
    string title2 = "A_{||,phys} for "+sector+"; X; A_{||,phys}";
    mgAllPhys_NH3->SetTitle(title2.c_str());
    mgAllPhys_NH3->GetXaxis()->SetLimits(0,0.8); //mgAllPhys_NH3->GetYaxis()->SetRangeUser(0.0,0.4);
    mgAllPhys_NH3->GetYaxis()->SetRangeUser(0,1.2);
    string title3 = "A_{1,p}(X) for "+sector+"; X; A_{1,p}";
    mgA1_NH3_Q2Avg->SetTitle(title3.c_str());
    mgA1_NH3_Q2Avg->GetXaxis()->SetLimits(0,0.8);
    mgA1_NH3_Q2Avg->GetYaxis()->SetRangeUser(0,1.2);
    string title4 = "A_{||,raw} for "+sector+"; X; A_{||,raw}";
    mgAllRaw_NH3->SetTitle(title4.c_str());
    mgAllRaw_NH3->GetXaxis()->SetLimits(0,0.8); //mgAllRaw_NH3->GetYaxis()->SetRangeUser(0.0,0.4);
    mgAllRaw_NH3->GetYaxis()->SetRangeUser(0,0.14);

    // Make a line showing the SU(6) prediction
    TF1* SU6 = new TF1("SU6","5.0/9.0",-1,1);
    SU6->SetLineWidth(3); SU6->SetLineColor(kBlack);
    SU6->SetLineStyle(2);

    //return mgA1_NH3;
    string c1title = "c1_"+sector;
    string c2title = "c2_"+sector;
    string c3title = "c3_"+sector;
    string c4title = "c4_"+sector;

    TCanvas* c1 = new TCanvas(c1title.c_str(),c1title.c_str(),800,600);
    mgA1_NH3->Draw("AP"); leg->Draw("same"); SU6->Draw("same");

    TCanvas* c2 = new TCanvas(c2title.c_str(),c2title.c_str(),800,600);
    mgAllPhys_NH3->Draw("AP"); leg->Draw("same");

    TCanvas* c3 = new TCanvas(c3title.c_str(),c3title.c_str(),800,600);
    mgA1_NH3_Q2Avg->Draw("ap"); SU6->Draw("same");

    TCanvas* c4 = new TCanvas(c4title.c_str(),c4title.c_str(),800,600);
    mgAllRaw_NH3->Draw("ap"); leg->Draw("same");

    vector<TCanvas*> Plots = {c1, c2, c3, c4};
    
    return Plots;

    // return A1Plots;

}

// This makes a TGraphErrors object of all the target polarizations for the runs in a given epoch
// Make sure that the DFs have been calculated for the "DF_Data" set!!!
/*
TGraphErrors* Run_Range_PbPt( vector<int>& Epoch, DataSet& DF_Data, string Target, RunPeriod& Period ){

    vector<double> runs, pbpt, pbptErr, zeros; // Used for plotting
    for(size_t i=0; i<Epoch.size(); i++){
	PbPt thisTpol;
	thisTpol.ReadInThisRun( Epoch[i], Target, Period );
	thisTpol.SetDFs( DF_Data );
	thisTpol.CalculateAllRaw();
	thisTpol.CalculatePbPt();
	double bt = thisTpol.getBT_Pol(); double btErr = thisTpol.getBT_Pol_Err();
	if( bt != 0.0 && btErr > 0.0 ){
	    runs.push_back( Epoch[i] ); pbpt.push_back( bt );
	    pbptErr.push_back( btErr ); zeros.push_back(0.0);
	}
    }
    if( runs.size() > 0 ){
	//string name = to_string(runs[i]) + "-" + to_string(runs[ runs.size()-1 ]);
	TGraphErrors* PbPt_Vals = new TGraphErrors( runs.size(), runs.data(), pbpt.data(), zeros.data(), pbptErr.data() );
	return PbPt_Vals;
    }
    else{
	cout <<"ERROR: Couldn't calculate any PbPt for this run range. Returning empty TGraphErrors...\n";
	TGraphErrors* empty = new TGraphErrors();
	return empty;
    }

}
*/
// Plots PbPt over several epochs using the elastic method
TCanvas* All_Epochs_PbPt( vector<int>& Epoch, vector<DataSet>& DF_Data, string Target ){

    vector<double> runs, pbpt, pbptErr, zeros, pt, pterr; // Used for plotting
    // Safety check to make sure that DF_Data and Epoch are the same size (i.e. the runs are lined up)
    if( Epoch.size() != DF_Data.size() ){
	cerr <<"ERROR: Runs aren't lined up with the DataSet objects. Check inputs and try again.\n";
	TCanvas* err;
	return err;
    }

    // Used for plotting the weighted average PbPt for each of the runs
    double Su22pbptPosNum = 0; double Su22pbptPosDen = 0;
    double Su22pbptNegNum = 0; double Su22pbptNegDen = 0;
    double Fa22NegpbptPosNum = 0; double Fa22NegpbptPosDen = 0;
    double Fa22NegpbptNegNum = 0; double Fa22NegpbptNegDen = 0;
    double Fa22PospbptPosNum = 0; double Fa22PospbptPosDen = 0;
    double Fa22PospbptNegNum = 0; double Fa22PospbptNegDen = 0;
    double Sp23pbptPosNum = 0; double Sp23pbptPosDen = 0;
    double Sp23pbptNegNum = 0; double Sp23pbptNegDen = 0;

    for(size_t i=0; i<Epoch.size(); i++){ // Loop over epochs
	double bt = DF_Data[i].getBT_Pol(); double btErr = DF_Data[i].getBT_Pol_Err();
	if( bt != 0.0 && btErr > 0.0 ){
	    runs.push_back( Epoch[i] ); 
	    pbpt.push_back( bt );
	    pbptErr.push_back( btErr ); zeros.push_back(0.0);
	    vector<double> truePt = CalculateTruePt( Epoch[i], bt, btErr );
	    pt.push_back( truePt[0] );
	    pterr.push_back( truePt[1] );

	    // Based on the run number, add the value of PbPt to the sum terms for calculating the weighted averages
	    // Summer data
	    if(      Epoch[i] >= 16137 && Epoch[i] <= 16772 && bt > 0 ){ Su22pbptPosNum += bt/(btErr*btErr); Su22pbptPosDen += 1.0/(btErr*btErr); }
	    else if( Epoch[i] >= 16137 && Epoch[i] <= 16772 && bt < 0 ){ Su22pbptNegNum += bt/(btErr*btErr); Su22pbptNegDen += 1.0/(btErr*btErr); }
	    // Fall negative solenoid data
	    else if( Epoch[i] >= 16859 && Epoch[i] <= 17183 && bt > 0 ){ Fa22NegpbptPosNum += bt/(btErr*btErr); Fa22NegpbptPosDen += 1.0/(btErr*btErr); }
	    else if( Epoch[i] >= 16859 && Epoch[i] <= 17183 && bt < 0 ){ Fa22NegpbptNegNum += bt/(btErr*btErr); Fa22NegpbptNegDen += 1.0/(btErr*btErr); }
	    // Fall positive solenoid data
	    else if( Epoch[i] >= 17188 && Epoch[i] <= 17408 && bt > 0 ){ Fa22PospbptPosNum += bt/(btErr*btErr); Fa22PospbptPosDen += 1.0/(btErr*btErr); }
	    else if( Epoch[i] >= 17188 && Epoch[i] <= 17408 && bt < 0 ){ Fa22PospbptNegNum += bt/(btErr*btErr); Fa22PospbptNegDen += 1.0/(btErr*btErr); }
	    // Spring inbending data
	    else if( Epoch[i] >= 17477 && Epoch[i] <= 17768 && bt > 0 ){ Sp23pbptPosNum += bt/(btErr*btErr); Sp23pbptPosDen += 1.0/(btErr*btErr); }
	    else if( Epoch[i] >= 17477 && Epoch[i] <= 17768 && bt < 0 ){ Sp23pbptNegNum += bt/(btErr*btErr); Sp23pbptNegDen += 1.0/(btErr*btErr); }
	}
    }
   
    if( runs.size() > 0 ){

        // Now calculate the weighted averages for each of the run periods. I didn't add safety checks for zero divisions cuz I'm lazy lol
	// Summer data
        double wAvgSu22_PbPt_Pos    = Su22pbptPosNum / Su22pbptPosDen;
        double wAvgSu22_PbPt_PosErr = sqrt( 1.0 / Su22pbptPosDen );
	TF1* wAvgSu22PosLine = new TF1("wAvgSu22PosLine", "pol0", 16000, 16800); 
	wAvgSu22PosLine->SetParameter(0, wAvgSu22_PbPt_Pos);
	wAvgSu22PosLine->SetLineColor(kBlue);
	double PNF_Su22_Pos    = 0.71 / wAvgSu22_PbPt_Pos;
	double errPNF_Su22_Pos = 0.02 / wAvgSu22_PbPt_Pos;
	cout <<"Average Summer DIS PbPt = "<< wAvgSu22_PbPt_Pos <<" +- "<< wAvgSu22_PbPt_PosErr << endl;
	cout <<"Average Summer PNF = "<< PNF_Su22_Pos <<" +- "<< errPNF_Su22_Pos << endl;

        double wAvgSu22_PbPt_Neg    = Su22pbptNegNum / Su22pbptNegDen;
        double wAvgSu22_PbPt_NegErr = sqrt( 1.0 / Su22pbptNegDen );
	TF1* wAvgSu22NegLine = new TF1("wAvgSu22NegLine", "pol0", 16000, 16800); 
	wAvgSu22NegLine->SetParameter(0, wAvgSu22_PbPt_Neg);
	wAvgSu22NegLine->SetLineColor(kBlue);
	double PNF_Su22_Neg    = -0.66 / wAvgSu22_PbPt_Neg;
	double errPNF_Su22_Neg =  0.03 / wAvgSu22_PbPt_Neg;
	cout <<"Average Summer DIS PbPt = "<< wAvgSu22_PbPt_Neg <<" +- "<< wAvgSu22_PbPt_NegErr << endl;
	cout <<"Average Summer PNF = "<< PNF_Su22_Neg <<" +- "<< errPNF_Su22_Neg << endl;
	cout << endl;
	
	// Fall negative solenoid data
        double wAvgFa22Neg_PbPt_Pos    = Fa22NegpbptPosNum / Fa22NegpbptPosDen;
        double wAvgFa22Neg_PbPt_PosErr = sqrt( 1.0 / Fa22NegpbptPosDen );
	TF1* wAvgFa22NegPosLine = new TF1("wAvgFa22NegPosLine", "pol0", 16850, 17183); 
	wAvgFa22NegPosLine->SetParameter(0, wAvgFa22Neg_PbPt_Pos);
	wAvgFa22NegPosLine->SetLineColor(kBlue);
	double PNF_Fa22Neg_Pos    = 0.71 / wAvgFa22Neg_PbPt_Pos;
	double errPNF_Fa22Neg_Pos = 0.02 / wAvgFa22Neg_PbPt_Pos;
	cout <<"Average Fall (Neg.) DIS PbPt = "<< wAvgFa22Neg_PbPt_Pos <<" +- "<< wAvgFa22Neg_PbPt_PosErr << endl;
	cout <<"Average Fall (Neg.) PNF = "<< PNF_Fa22Neg_Pos <<" +- "<< errPNF_Fa22Neg_Pos << endl;

        double wAvgFa22Neg_PbPt_Neg    = Fa22NegpbptNegNum / Fa22NegpbptNegDen;
        double wAvgFa22Neg_PbPt_NegErr = sqrt( 1.0 / Fa22NegpbptNegDen );
	TF1* wAvgFa22NegNegLine = new TF1("wAvgFa22NegNegLine", "pol0", 16850, 17183); 
	wAvgFa22NegNegLine->SetParameter(0, wAvgFa22Neg_PbPt_Neg);
	wAvgFa22NegNegLine->SetLineColor(kBlue);
	double PNF_Fa22Neg_Neg    = -0.66 / wAvgFa22Neg_PbPt_Neg;
	double errPNF_Fa22Neg_Neg =  0.03 / wAvgFa22Neg_PbPt_Neg;
	cout <<"Average Fall (Neg.) DIS PbPt = "<< wAvgFa22Neg_PbPt_Neg <<" +- "<< wAvgFa22Neg_PbPt_NegErr << endl;
	cout <<"Average Fall (Neg.) PNF = "<< PNF_Fa22Neg_Neg <<" +- "<< errPNF_Fa22Neg_Neg << endl;
	cout << endl;

	// Fall positive solenoid data
        double wAvgFa22Pos_PbPt_Pos    = Fa22PospbptPosNum / Fa22PospbptPosDen;
        double wAvgFa22Pos_PbPt_PosErr = sqrt( 1.0 / Fa22PospbptPosDen );
	TF1* wAvgFa22PosPosLine = new TF1("wAvgFa22PosPosLine", "pol0", 17185, 17450); 
	wAvgFa22PosPosLine->SetParameter(0, wAvgFa22Pos_PbPt_Pos);
	wAvgFa22PosPosLine->SetLineColor(kBlue);
	double PNF_Fa22Pos_Pos    = 0.71 / wAvgFa22Pos_PbPt_Pos;
	double errPNF_Fa22Pos_Pos = 0.02 / wAvgFa22Pos_PbPt_Pos;
	cout <<"Average Fall (Pos.) DIS PbPt = "<< wAvgFa22Pos_PbPt_Pos <<" +- "<< wAvgFa22Pos_PbPt_PosErr << endl;
	cout <<"Average Fall (Pos.) PNF = "<< PNF_Fa22Pos_Pos <<" +- "<< errPNF_Fa22Pos_Pos << endl;

        double wAvgFa22Pos_PbPt_Neg    = Fa22PospbptNegNum / Fa22PospbptNegDen;
        double wAvgFa22Pos_PbPt_NegErr = sqrt( 1.0 / Fa22PospbptNegDen );
	TF1* wAvgFa22PosNegLine = new TF1("wAvgFa22PosNegLine", "pol0", 17185, 17450); 
	wAvgFa22PosNegLine->SetParameter(0, wAvgFa22Pos_PbPt_Neg);
	wAvgFa22PosNegLine->SetLineColor(kBlue);
	double PNF_Fa22Pos_Neg    = -0.66 / wAvgFa22Pos_PbPt_Neg;
	double errPNF_Fa22Pos_Neg =  0.03 / wAvgFa22Pos_PbPt_Neg;
	cout <<"Average Fall (Pos.) DIS PbPt = "<< wAvgFa22Pos_PbPt_Neg <<" +- "<< wAvgFa22Pos_PbPt_NegErr << endl;
	cout <<"Average Fall (Pos.) PNF = "<< PNF_Fa22Pos_Neg <<" +- "<< errPNF_Fa22Pos_Neg << endl;
	cout << endl;

	// Spring data
        double wAvgSp23_PbPt_Pos    = Sp23pbptPosNum / Sp23pbptPosDen;
        double wAvgSp23_PbPt_PosErr = sqrt( 1.0 / Sp23pbptPosDen );
	TF1* wAvgSp23PosLine = new TF1("wAvgSp23PosLine", "pol0", 17500, 18000); 
	wAvgSp23PosLine->SetParameter(0, wAvgSp23_PbPt_Pos);
	wAvgSp23PosLine->SetLineColor(kBlue);
	double PNF_Sp23_Pos    = 0.71 / wAvgSp23_PbPt_Pos;
	double errPNF_Sp23_Pos = 0.02 / wAvgSp23_PbPt_Pos;
	cout <<"Average Spring DIS PbPt = "<< wAvgSp23_PbPt_Pos <<" +- "<< wAvgSp23_PbPt_PosErr << endl;
	cout <<"Average Spring PNF = "<< PNF_Sp23_Pos <<" +- "<< errPNF_Sp23_Pos << endl;

        double wAvgSp23_PbPt_Neg    = Sp23pbptNegNum / Sp23pbptNegDen;
        double wAvgSp23_PbPt_NegErr = sqrt( 1.0 / Sp23pbptNegDen );
	TF1* wAvgSp23NegLine = new TF1("wAvgSp23NegLine", "pol0", 17500, 18000); 
	wAvgSp23NegLine->SetParameter(0, wAvgSp23_PbPt_Neg);
	wAvgSp23NegLine->SetLineColor(kBlue);
	double PNF_Sp23_Neg    = -0.66 / wAvgSp23_PbPt_Neg;
	double errPNF_Sp23_Neg =  0.03 / wAvgSp23_PbPt_Neg;
	cout <<"Average Spring DIS PbPt = "<< wAvgSp23_PbPt_Neg <<" +- "<< wAvgSp23_PbPt_NegErr << endl;
	cout <<"Average Spring PNF = "<< PNF_Sp23_Neg <<" +- "<< errPNF_Sp23_Neg << endl;
	cout << endl;

	// Now that the PNF values have been calculated, create scaled values of the DIS PbPt values using this PNF
	vector<double> RunsScaled, PbPtScaled, ZerosScaled, ErrPbPtScaled;
	for(size_t i=0; i<Epoch.size(); i++){ // Loop over epochs
	    double bt = DF_Data[i].getBT_Pol(); double btErr = DF_Data[i].getBT_Pol_Err();
	    // Summer data
	    if(      Epoch[i] >= 16137 && Epoch[i] <= 16772 && bt > 0 ){ bt = bt * PNF_Su22_Pos; btErr = sqrt( pow(PNF_Su22_Pos,2)*btErr*btErr + bt*bt*pow(errPNF_Su22_Pos,2) ); }
	    else if( Epoch[i] >= 16137 && Epoch[i] <= 16772 && bt < 0 ){ bt = bt * PNF_Su22_Neg; btErr = sqrt( pow(PNF_Su22_Neg,2)*btErr*btErr + bt*bt*pow(errPNF_Su22_Neg,2) ); }
	    // Fall negative solenoid data
	    else if( Epoch[i] >= 16859 && Epoch[i] <= 17183 && bt > 0 ){ bt = bt * PNF_Fa22Neg_Pos; btErr = sqrt( pow(PNF_Fa22Neg_Pos,2)*btErr*btErr + bt*bt*pow(errPNF_Fa22Neg_Pos,2) ); }
	    else if( Epoch[i] >= 16859 && Epoch[i] <= 17183 && bt < 0 ){ bt = bt * PNF_Fa22Neg_Neg; btErr = sqrt( pow(PNF_Fa22Neg_Neg,2)*btErr*btErr + bt*bt*pow(errPNF_Fa22Neg_Neg,2) ); }
	    // Fall positive solenoid data
	    else if( Epoch[i] >= 17188 && Epoch[i] <= 17408 && bt > 0 ){ bt = bt * PNF_Fa22Pos_Pos; btErr = sqrt( pow(PNF_Fa22Pos_Pos,2)*btErr*btErr + bt*bt*pow(errPNF_Fa22Pos_Pos,2) ); }
	    else if( Epoch[i] >= 17188 && Epoch[i] <= 17408 && bt < 0 ){ bt = bt * PNF_Fa22Pos_Neg; btErr = sqrt( pow(PNF_Fa22Pos_Neg,2)*btErr*btErr + bt*bt*pow(errPNF_Fa22Pos_Neg,2) ); }
	    // Spring inbending data
	    else if( Epoch[i] >= 17477 && Epoch[i] <= 17768 && bt > 0 ){ bt = bt * PNF_Sp23_Pos; btErr = sqrt( pow(PNF_Sp23_Pos,2)*btErr*btErr + bt*bt*pow(errPNF_Sp23_Pos,2) ); }
	    else if( Epoch[i] >= 17477 && Epoch[i] <= 17768 && bt < 0 ){ bt = bt * PNF_Sp23_Neg; btErr = sqrt( pow(PNF_Sp23_Neg,2)*btErr*btErr + bt*bt*pow(errPNF_Sp23_Neg,2) ); }
	    
	    RunsScaled.push_back( Epoch[i] ); PbPtScaled.push_back( bt ); ErrPbPtScaled.push_back( btErr ); ZerosScaled.push_back( 0 );
	}

	// Now set the error band
	/*
	TGraphErrors* grSu22PosErr = new TGraphErrors();
	grSu22PosErr->SetPoint(0, 16000, wAvgSu22_PbPt_Pos ); grSu22PosErr->SetPointError(0, 0, wAvgSu22_PbPt_PosErr);
	grSu22PosErr->SetPoint(1, 16800, wAvgSu22_PbPt_Pos ); grSu22PosErr->SetPointError(1, 0, wAvgSu22_PbPt_PosErr);
	grSu22PosErr->SetFillColorAlpha(kBlue,0.35);
	*/
        //double wAvgSu22_PbPt_Neg = Su22pbptNegNum / Su22pbptNegDen;

	//string name = to_string(runs[i]) + "-" + to_string(runs[ runs.size()-1 ]);
	TGraphErrors* PbPt_Vals = new TGraphErrors( runs.size(), runs.data(), pbpt.data(), zeros.data(), pbptErr.data() );
	PbPt_Vals->SetMarkerColor(kBlue); PbPt_Vals->SetMarkerStyle(kFullCircle);

	TGraphErrors* TruePt_Vals = new TGraphErrors( runs.size(), runs.data(), pt.data(), zeros.data(), pterr.data() );
	TruePt_Vals->SetMarkerColor(kViolet); TruePt_Vals->SetMarkerStyle(kFullSquare);

	TGraphErrors* grScaledPbPt = new TGraphErrors( RunsScaled.size(), RunsScaled.data(), PbPtScaled.data(), ZerosScaled.data(), ErrPbPtScaled.data() );
	grScaledPbPt->SetMarkerColor(kViolet); grScaledPbPt->SetMarkerStyle(kFullSquare);

	string canName = Target+"_PbPtVals";
	TCanvas* c = new TCanvas(canName.c_str(), canName.c_str(), 1600,1200);

	TLegend* leg = new TLegend(0.7,0.4,0.9,0.6);
	TMultiGraph* mg = new TMultiGraph();
	string mgTitle = Target+" Beam-Target Polarizations; Run Number; PbPt";
	mg->SetTitle( mgTitle.c_str() );
	mg->Add( PbPt_Vals, "p");
	mg->Add( grScaledPbPt, "p");
	mg->GetYaxis()->SetRangeUser(-1,1.2);
	leg->AddEntry( PbPt_Vals, "DIS PbPt","p");
	leg->AddEntry( grScaledPbPt, "Scaled DIS PbPt","p");
	leg->AddEntry( wAvgSu22PosLine, "Avg. DIS PbPt", "lp" );

	//mg->Add( grSu22PosErr, "e3 al" );
	/*
	mg->Add( TruePt_Vals, "p" );
	leg->AddEntry( TruePt_Vals, "(DIS PbPt)/Pb_{avg}", "p");
	*/
	// Now draw the lines;
	
	TF1* Su22PosLine = new TF1("Su22PosLine", "0.71",16000,16800);
	leg->AddEntry( Su22PosLine, "Elastic PbPt", "lp");

	TF1* Su22NegLine = new TF1("Su22NegLine","-0.66",16000,16800);
	TF1* Fa22NSolPosLine = new TF1("Fa22NSolPosLine", "0.72",16850,17183);
	TF1* Fa22NSolNegLine = new TF1("Fa22NSolNegLine","-0.69",16850,17183);
	TF1* Fa22PSolPosLine = new TF1("Fa22PSolPosLine", "0.67",17185,17450);
	TF1* Fa22PSolNegLine = new TF1("Fa22PSolNegLine","-0.70",17185,17450);
	TF1* Sp23PosLine = new TF1("Sp23PosLine", "0.67",17500,18000);
	TF1* Sp23NegLine = new TF1("Sp23NegLine","-0.62",17500,18000);
	TF1* BS_Line_Pos = new TF1("BS_Line_Pos","1",16000,18000); BS_Line_Pos->SetLineColor(kBlack);
	TF1* BS_Line_Neg = new TF1("BS_Line_Neg","-1",16000,18000); BS_Line_Neg->SetLineColor(kBlack);

	c->cd();
	mg->Draw("AP");
	Su22PosLine->Draw("same"); Su22NegLine->Draw("same");
	Fa22NSolPosLine->Draw("same"); Fa22NSolNegLine->Draw("same");
	Fa22PSolPosLine->Draw("same"); Fa22PSolNegLine->Draw("same");
	Sp23PosLine->Draw("same"); Sp23NegLine->Draw("same");
	BS_Line_Pos->Draw("same"); BS_Line_Neg->Draw("same");
	wAvgSu22PosLine->Draw("same"); wAvgSu22NegLine->Draw("same");
	wAvgFa22NegPosLine->Draw("same"); wAvgFa22NegNegLine->Draw("same");
	wAvgFa22PosPosLine->Draw("same"); wAvgFa22PosNegLine->Draw("same");
	wAvgSp23PosLine->Draw("same"); wAvgSp23NegLine->Draw("same");

	leg->Draw("same");

	return c;	
	//mg->Add( NMR_PbPt_Vals, "p");

    }
    else{
	cout <<"ERROR: Couldn't calculate any PbPt for this run range. Returning empty TCanvas...\n";
	//TMultiGraph* empty = new TMultiGraph();
	TCanvas* empty;
	return empty;
    }

}

// Using the maximum likelihood method, this function calculates the physical asymmetry across a set of
// runs (given in the "Epoch" vector) and the corresponding DataSet objects in "DF_Data".
TCanvas* All_Epochs_All_Phys( vector<int>& Epoch, vector<DataSet>& DF_Data, string Target, string pathToData, string period="All" ){

    vector<int> palette = { 1, 632, 800, 400, 416, 600, 840, 880, 900, 616, 820, 432, 920, 860 };

    // Safety check to make sure that DF_Data and Epoch are the same size (i.e. the runs are lined up)
    if( Epoch.size() != DF_Data.size() ){
	cerr <<"ERROR: Runs aren't lined up with the DataSet objects. Check inputs and try again.\n";
	TCanvas* err;
	return err;
    }

    // This is a DataSet object that will hold all of the output A_meas values
    DataSet All_Data;
    //FIXME: Replace this with "CalculateAvgXQ2()" function
    //All_Data.SetAvgXQ2( Epoch, pathToData, "NH3" ); // Set the average values of x, Q2
						    // using all available runs

    string cName = "c_"+period;
    TCanvas* c = new TCanvas(cName.c_str(),cName.c_str(),800,600);
    TLegend* leg = new TLegend(0.1,0.7,0.5,0.9); leg->SetNColumns(3);
    leg->SetHeader("Q^{2} Bin Centers [GeV^{2}]","C");
    TMultiGraph* mg = new TMultiGraph();
    TMultiGraph* mgA1 = new TMultiGraph(); // Plot the A1 values
	
    /*
    // Plots the A1 values averaged over Q2 bins
    TCanvas* cAvg = new TCanvas("cAvg","cAvg",800,600);
    TLegend* legAvg = new TLegend(0.1,0.7,0.5,0.9);
    legAvg->SetHeader("Q^{2} Bin Centers [GeV^{2}]","C");
    TMultiGraph* mgAvg = new TMultiGraph(); 
*/
    // Output file with all the calculated A1 values
    string outputName ="Calculated_A1_Values_"+period+".txt";
    ofstream fout( outputName.c_str() );
    fout <<"x  Q2  A1  A1_error  A1_theory\n";

    // Calculate a max. likelihood value of All,meas for each bin in x, Q2
    for(size_t i=0; i<Q2_Bin_Bounds.size()-1; i++){
	double qmid = (Q2_Bin_Bounds[i]+Q2_Bin_Bounds[i+1])/2.0;
	vector<double> xMeans, allPhys, zeros, allPhysErr, a1Data, a1DataErr, a1TheoryNH3;

	// Used for plotting the expected systematic error band
	double nSteps = X_Bin_Bounds.size()-1;

	// Used to calculate the weighted average of the model for this Q2 bin
	//double modelAvgA1 = 0; double modelCounts = 0;

	for(size_t j=0; j<X_Bin_Bounds.size()-1; j++){
	    double xmid = (X_Bin_Bounds[j]+X_Bin_Bounds[j+1])/2.0;

	    Bin thisAvgAllData = All_Data.getThisBin( qmid, xmid );
	    // Now loop over all the runs
	    double aMeasNum = 0.0; // numerator term in the ML eq.
	    double aMeasDen = 0.0; // denominator term in the ML eq.
	    //cout << "Looking at bin Q2 = "<< qmid <<", X = "<< xmid << endl;
	    for(int k=0; k<Epoch.size(); k++){
		// Get the raw counts un-normalized to FC charge
		Bin thisBin = DF_Data[k].getThisBin( qmid, xmid );
		double NM = thisBin.getNP( Target );
		double NP = thisBin.getNM( Target );
		double FC_P= thisBin.getFC_P( Target );
		double FC_M= thisBin.getFC_M( Target );
		double DF = DF_Data[k].getDF_NH3( qmid, xmid );
		double P  = DF_Data[k].getBT_Pol();
		if( NM > 0 && NP > 0 && DF > 0 && P > 0 && FC_P > 0 && FC_M > 0 ){
		    //NM = NM / FC_P;
		    //NP = NP / FC_M;
		    aMeasNum += (NM - NP)*DF*P;
		    aMeasDen += (NM + NP)*DF*DF*P*P;
		}
	    }
	    if( aMeasNum > 0 && aMeasDen > 0 ){
		double err = sqrt( 1.0 / aMeasDen );
		double aMeasThisBin = aMeasNum / aMeasDen;
		// Get the mean value of x
		xMeans.push_back( All_Data.getAvgX( Target, qmid, xmid ) );
		allPhys.push_back( aMeasThisBin );
		zeros.push_back( 0 );
		allPhysErr.push_back( err );
		//The A1 data
		if( thisAvgAllData.getDepolFactorNH3() > 0.0001 ){
		    double a1 = aMeasThisBin/thisAvgAllData.getDepolFactorNH3() - thisAvgAllData.getEtaNH3()*thisAvgAllData.getA2NH3();
		    double a1Err = err / thisAvgAllData.getDepolFactorNH3();
		    //cout << "A1 = "<< a1 <<" +- "<< a1Err << endl;
		    a1Data.push_back( a1 );
		    a1DataErr.push_back( a1Err );
		    double a1Theory = thisAvgAllData.getA1Theory("NH3"); 
		    if( a1Theory > 0.0001 ) a1TheoryNH3.push_back( a1Theory );
		    fout << All_Data.getAvgX( Target, qmid, xmid ) <<"  "<< All_Data.getAvgQ2( Target, qmid, xmid ) <<"  "<< a1 <<"  "<< a1Err <<"  "<< a1Theory<< endl;

		    //modelAvgA1 += thisAvgAllData.getA1Theory("NH3"); modelCounts++;
		}
	    }
	}

	// Calculate the model average value
	TGraphErrors* gr = new TGraphErrors( xMeans.size(), xMeans.data(), allPhys.data(), zeros.data(), allPhysErr.data() );
	gr->SetMarkerStyle( kFullCircle );
	gr->SetMarkerColor( palette[i] );
	stringstream legEntry; legEntry << "Q^{2} = "<< setprecision(3) << qmid;
	leg->AddEntry( gr, legEntry.str().c_str(), "p" );
	mg->Add( gr, "p" );

	TGraphErrors* grA1 = new TGraphErrors( xMeans.size(), xMeans.data(), a1Data.data(), zeros.data(), a1DataErr.data() );
	grA1->SetMarkerStyle( kFullCircle );
	grA1->SetMarkerColor( palette[i] );
	mgA1->Add( grA1, "p" );
/*
	TGraphErrors* grA1Theory = new TGraphErrors( xMeans.size(), xMeans.data(), a1TheoryNH3.data(), zeros.data(), zeros.data() );
	grA1Theory->SetMarkerStyle( kFullSquare );
	grA1Theory->SetMarkerColor( kGray );
	mgA1->Add( grA1Theory, "p" );
*/

    } // End of Q2_Bin_Bounds loop

    fout.close();

    c->cd();
    /*
    mg->SetTitle("Plot of All A_{||,meas}; X; A_{||,meas}");
    mg->Draw("ap");
    leg->Draw("same");
    */
    string a1Title = "Plot of "+period+" A_{1}^{p}; X; A_{1}^{p}";
    mgA1->SetTitle( a1Title.c_str() );
    mgA1->GetYaxis()->SetRangeUser(0,1);
    mgA1->Draw("ap");
    leg->Draw("same");

    return c;
}

// This function reads in the output from "Calculated_A1_Values.txt" and makes a plot of the weighted
// xbin values of A1 for the data and the theory
TCanvas* A1_Q2_Averages(string period){

    string cName = "a1Plt_"+period;
    TCanvas* a1Plt = new TCanvas(cName.c_str(),cName.c_str(),800,600);
    TMultiGraph* mg = new TMultiGraph();
    TLegend* leg = new TLegend(0.1,0.7,0.4,0.9);
    //leg->SetHeader("Q^{2} Bins [GeV^{2}]","C");

    vector<string> FileRows;
    string infileName = "Calculated_A1_Values_"+period+".txt";
    ifstream fin( infileName.c_str() );
    string line;
    getline(fin,line); // Throw away header row
    while( getline(fin,line) ) FileRows.push_back( line );
    fin.close();

    vector<double> xVals, a1Vals, zeros, a1Errs, a1Models, sysErrBar, sysErrorBand;

    // Used for the sysErrorBand plot
    double nSteps = X_Bin_Bounds.size()-1;

    for(size_t i=0; i<X_Bin_Bounds.size()-1; i++){
        double xmid = ( X_Bin_Bounds[i] + X_Bin_Bounds[i+1] ) / 2.0;
        double xAvg = 0; 
	double a1AvgNum = 0; 
	double a1AvgDen = 0;
	double a1ModAvg = 0;
	double counts = 0;
        for(size_t j=0; j<Q2_Bin_Bounds.size()-1; j++){

	    double qmid = (Q2_Bin_Bounds[j] + Q2_Bin_Bounds[j+1]) / 2.0;
	    // Loop over the data to find the correct bin
	    for(auto row : FileRows ){
		stringstream sin(row);
		double x, q2, a1, a1err, a1theory;
		sin >> x >> q2 >> a1 >> a1err >> a1theory;
		if( x < X_Bin_Bounds[i+1] && x >= X_Bin_Bounds[i] && q2 < Q2_Bin_Bounds[j+1] && q2 >= Q2_Bin_Bounds[j] ){
		    xAvg += x; 
		    a1ModAvg += a1theory;
		    a1AvgNum += a1 / (a1err*a1err);
		    a1AvgDen += 1.0/ (a1err*a1err);
		    counts++;
		}
	    }
        }
	if( counts > 0 && a1AvgDen > 0 ){
	    xVals.push_back( xAvg / counts );
	    a1Vals.push_back( a1AvgNum / a1AvgDen );
	    a1Errs.push_back( sqrt( 1.0 / a1AvgDen ) );
	    zeros.push_back( 0 );
	    sysErrBar.push_back( 0.1 );
	    a1Models.push_back( a1ModAvg / counts );
	    double sysErr = ( 0.11-(i*0.05/(1.0*nSteps)) )*( a1AvgNum / a1AvgDen ); 
	    sysErrorBand.push_back( sysErr );
	}

    }

    // Now make the plot
    TGraphAsymmErrors* sysErr = new TGraphAsymmErrors(xVals.size(), xVals.data(), sysErrBar.data(), zeros.data(), zeros.data(), zeros.data(), sysErrorBand.data() );
    sysErr->SetMarkerStyle(1); // make an invisible point
    sysErr->SetFillColor(kOrange);
    mg->Add( sysErr, "4" );
    leg->AddEntry( sysErr, "Systematic Error", "f" );

    TGraphErrors* grA1 = new TGraphErrors(xVals.size(), xVals.data(), a1Vals.data(), zeros.data(), a1Errs.data() );
    grA1->SetMarkerStyle(kFullCircle);
    grA1->SetMarkerColor(kBlue);
    mg->Add( grA1, "p");
    leg->AddEntry( grA1, "Q^{2}-Avg. A_{1}^{p}", "p");

    TGraph* grMod = new TGraph( xVals.size(), xVals.data(), a1Models.data() );
    grMod->SetMarkerStyle(1);
    grMod->SetLineColor(kRed);
    grMod->SetLineWidth(2);
    mg->Add( grMod, "l");
    leg->AddEntry( grMod, "Q^{2}-Avg. Model A_{1}^{p}", "l");

    string title = "A_{1}^{p} for "+ period +" RG-C Data; X; A_{1}^{p}";
    mg->SetTitle(title.c_str());
    mg->GetYaxis()->SetRangeUser(0,1);
    mg->GetXaxis()->SetRangeUser(0,1);

    // Make a line showing the SU(6) prediction
    TF1* SU6 = new TF1("SU6","5.0/9.0",-1,1);
    SU6->SetLineWidth(3); SU6->SetLineColor(kBlack);
    SU6->SetLineStyle(2);

    a1Plt->cd();
    mg->Draw("ap");
    leg->Draw("same");
    SU6->Draw("same");

    return a1Plt;

}	
/*
// This is jank, but I loop over the bins in a separate order to average over the Q2 bins.
*/
/*
// Plots PbPt over several epochs
TMultiGraph* All_Epochs_PbPt( vector<vector<int>>& Epoch, vector<DataSet>& DF_Data, string Target, RunPeriod& Period, bool useElastic=false ){

    vector<double> runs, pbpt, pbptErr, zeros, NMRpbpt, NMRpbptErr; // Used for plotting
    for(size_t i=0; i<Epoch.size(); i++){ // Loop over epochs
     for(size_t j=0; j<Epoch[i].size(); j++){ // Loop over the runs within an epoch
	PbPt thisTpol;
	cout << "Starting PbPt for run "<< Epoch[i][j] << endl;
	thisTpol.ReadInThisRun( Epoch[i][j], Target, Period );
	thisTpol.SetDFs( DF_Data[i] ); // This assumes, of course, that the epochs in Epoch and DF_Data line up, so be careful of that...
	thisTpol.CalculateAllRaw();
	if( useElastic ) thisTpol.MaxLikelihoodPbPt();
	else thisTpol.CalculatePbPt();

	double bt = thisTpol.getBT_Pol(); double btErr = thisTpol.getBT_Pol_Err();
	if( bt != 0.0 && btErr > 0.0 ){
	    runs.push_back( Epoch[i][j] ); pbpt.push_back( bt );
	    pbptErr.push_back( btErr ); zeros.push_back(0.0);
	    //NMRpbpt.push_back( 0.8 * Period.getTargetPolarization( Epoch[i][j] ) );
	    NMRpbpt.push_back( 0.8 * Period.getOfflineTPol( Epoch[i][j] ));
	    NMRpbptErr.push_back( 0.8 * Period.getOfflineTPolErr( Epoch[i][j] ));
	}
     }
    }
    if( runs.size() > 0 ){
	//string name = to_string(runs[i]) + "-" + to_string(runs[ runs.size()-1 ]);
	TGraphErrors* PbPt_Vals = new TGraphErrors( runs.size(), runs.data(), pbpt.data(), zeros.data(), pbptErr.data() );
	PbPt_Vals->SetMarkerColor(kBlue); PbPt_Vals->SetMarkerStyle(kFullCircle);
	//TGraphErrors* NMR_PbPt_Vals = new TGraphErrors( runs.size(), runs.data(), NMRpbpt.data(), zeros.data(), NMRpbptErr.data() );
	//NMR_PbPt_Vals->SetMarkerColor(kViolet); NMR_PbPt_Vals->SetMarkerStyle(kFullSquare);
	TMultiGraph* mg = new TMultiGraph();
	mg->Add( PbPt_Vals, "p");
	//mg->Add( NMR_PbPt_Vals, "p");
	
	return mg;
    }
    else{
	cout <<"ERROR: Couldn't calculate any PbPt for this run range. Returning empty TMultiGraph...\n";
	TMultiGraph* empty = new TMultiGraph();
	return empty;
    }
}
*/

// This function takes in a DataSet object and makes plots for Q2 bins
TMultiGraph* Make_NH3_Bin_Plots( DataSet& AllData, string sector ){
 
    vector<int> palette = { 1, 632, 800, 400, 416, 600, 840, 880, 900, 616, 820, 432, 920, 860 };
    int pint = 0;

    TMultiGraph* mgNH3 = new TMultiGraph();

    for(size_t i=0; i<Q2_Bin_Bounds.size()-1; i++){
	double qmid = (Q2_Bin_Bounds[i]+Q2_Bin_Bounds[i+1])/2.0;
	vector<double> xbins, dfnh3vals, zeros, dfnh3errs;
	for(size_t j=0; j<X_Bin_Bounds.size()-1; j++){
	    double xmid = (X_Bin_Bounds[j]+X_Bin_Bounds[j+1])/2.0;
	    double thisDFNH3 = AllData.getDF_NH3(qmid, xmid);
	    double thisErrDFNH3 = AllData.getErrDF_NH3(qmid, xmid);
	    if( thisDFNH3 != 0.0 ){
		xbins.push_back(xmid); zeros.push_back(0.0);
		dfnh3vals.push_back( thisDFNH3 );
		dfnh3errs.push_back( thisErrDFNH3 );
	    }
	}
	if( xbins.size() > 0 ){
	    TGraphErrors* gr = new TGraphErrors(xbins.size(), xbins.data(), dfnh3vals.data(), zeros.data(), dfnh3errs.data() );
	    gr->SetMarkerColor( palette[pint] );
	    if(pint < 4 ) gr->SetMarkerStyle(kFullCircle); 
	    else if(pint >= 4 && pint < 8) gr->SetMarkerStyle(kFullTriangleUp);
	    else if(pint >= 8) gr->SetMarkerStyle(kFullSquare);
	    mgNH3->Add( gr, "p" );
	    pint++;
	}
    }

    string title = "DF_{NH3} for "+sector+"; X; DF_{NH3}";
    mgNH3->SetTitle(title.c_str());
    mgNH3->GetXaxis()->SetLimits(0,0.8); mgNH3->GetYaxis()->SetRangeUser(0.0,0.4);

    return mgNH3;
}

TCanvas* NH3_Legend_Plot( DataSet& AllData, string sector ){
 
    vector<int> palette = { 1, 632, 800, 400, 416, 600, 840, 880, 900, 616, 820, 432, 920, 860 };
    int pint = 0;

    TMultiGraph* mgNH3 = new TMultiGraph();
    TLegend* mgLeg = new TLegend(0.1,0.7,0.4,0.9);
    mgLeg->SetNColumns(3);
    mgLeg->SetHeader("Q^{2} Bins (GeV^{2})","C");

    for(size_t i=0; i<Q2_Bin_Bounds.size()-1; i++){
	double qmid = (Q2_Bin_Bounds[i]+Q2_Bin_Bounds[i+1])/2.0;
	vector<double> xbins, dfnh3vals, zeros, dfnh3errs;
	for(size_t j=0; j<X_Bin_Bounds.size()-1; j++){
	    double xmid = (X_Bin_Bounds[j]+X_Bin_Bounds[j+1])/2.0;
	    double thisDFNH3 = AllData.getDF_NH3(qmid, xmid);
	    double thisErrDFNH3 = AllData.getErrDF_NH3(qmid, xmid);
	    if( thisDFNH3 != 0.0 ){
		xbins.push_back(xmid); zeros.push_back(0.0);
		dfnh3vals.push_back( thisDFNH3 );
		dfnh3errs.push_back( thisErrDFNH3 );
	    }
	}
	if( xbins.size() > 0 ){
	    TGraphErrors* gr = new TGraphErrors(xbins.size(), xbins.data(), dfnh3vals.data(), zeros.data(), dfnh3errs.data() );
	    gr->SetMarkerColor( palette[pint] );
	    if(pint < 4 ) gr->SetMarkerStyle(kFullCircle); 
	    else if(pint >= 4 && pint < 8) gr->SetMarkerStyle(kFullTriangleUp);
	    else if(pint >= 8) gr->SetMarkerStyle(kFullSquare);
	    mgNH3->Add( gr, "p" ); 
	    //string legTitle = "Q^{2}=" + to_string( trunc(qmid*1000)/1000 );
	    stringstream legTitle; legTitle << "Q^{2}=" << fixed << setprecision(3) << qmid;
	    mgLeg->AddEntry( gr, legTitle.str().c_str(), "p");
	    pint++;
	}
    }

    string title = "DF_{NH3} for "+sector+"; X; DF_{NH3}";
    mgNH3->SetTitle(title.c_str());
    mgNH3->GetXaxis()->SetLimits(0,0.8); mgNH3->GetYaxis()->SetRangeUser(0.0,0.4);

    string pltTitle = "NH3_Plot_"+sector;
    TCanvas* NH3_Plot = new TCanvas(pltTitle.c_str(),pltTitle.c_str(),800,600);
    NH3_Plot->cd();
    mgNH3->Draw("AP");
    mgLeg->Draw("same");
    return NH3_Plot;
}


// Makes a plot of the numerator term from the DF through PF calculation
/*
TCanvas* DF_Numerator_Legend_Plot( DataSet& AllData, string sector, string targetType ){

    if( !(targetType == "NH3" || targetType == "ND3") ){
	    cout <<"ERROR: enter a valid target type to plot: 'NH3' or 'ND3'\n";
	    TCanvas* emptyC = new TCanvas("emptyC","emptyC",800,600);
	    return emptyC;
    }
 
    vector<int> palette = { 1, 632, 800, 400, 416, 600, 840, 880, 900, 616, 820, 432, 920, 860 };
    int pint = 0;

    TMultiGraph* mgDF = new TMultiGraph();
    TLegend* mgLeg = new TLegend(0.1,0.7,0.4,0.9);
    mgLeg->SetNColumns(3);
    mgLeg->SetHeader("Q^{2} Bins (GeV^{2})","C");

    // Check to see if we're using the average PF across all bins (thruPF = "YES") or if the
    // PF is being calculated for each bin separately (thruPF = else)
    //if( thruPF == "YES" ){
	AllData.CalculateDFThruPF( targetType );
	cout << "Calculating DF thru PF for " << targetType << endl;
    //}
    //else
	//AllData.CalculateDF();

    for(size_t i=0; i<Q2_Bin_Bounds.size()-1; i++){
	double qmid = (Q2_Bin_Bounds[i]+Q2_Bin_Bounds[i+1])/2.0;
	vector<double> xbins, dfvals, zeros, dferrs;
	for(size_t j=0; j<X_Bin_Bounds.size()-1; j++){
	    double xmid = (X_Bin_Bounds[j]+X_Bin_Bounds[j+1])/2.0;
	    double thisDFnumerator = 0; double thisErrDFnumerator = 0;
	    
	    thisDFnumerator = AllData.getNumeratorDFThruPF( qmid, xmid, targetType );
	    thisErrDFnumerator = AllData.getNumeratorDFThruPFError( qmid, xmid, targetType );
	    
            if( thisDFnumerator != 0.0 ){
		    xbins.push_back(xmid); zeros.push_back(0.0);
		    dfvals.push_back( thisDFnumerator );
		    dferrs.push_back( thisErrDFnumerator );
	    }
	}
	if( xbins.size() > 0 ){
	    TGraphErrors* gr = new TGraphErrors(xbins.size(), xbins.data(), dfvals.data(), zeros.data(), dferrs.data() );
	    gr->SetMarkerColor( palette[pint] );
	    if(pint < 4 ){ gr->SetMarkerStyle(kFullCircle); }
	    else if(pint >= 4 && pint < 8){ gr->SetMarkerStyle(kFullTriangleUp); }
	    else if(pint >= 8){ gr->SetMarkerStyle(kFullSquare); }
	    mgDF->Add( gr, "p" );
	    stringstream legTitle; legTitle << "Q^{2}=" << fixed << setprecision(3) << qmid;
	    mgLeg->AddEntry( gr, legTitle.str().c_str(), "p");
	    pint++;
	}
    }

    string title = targetType+" DF_{"+targetType+"} 'Numerator' Term for "+sector+"; X; Numerator Term";
    mgDF->SetTitle(title.c_str());
    mgDF->GetXaxis()->SetLimits(0,0.8); //mgDF->GetYaxis()->SetRangeUser(0.0,0.5);
    
    string pltTitle = targetType+"_DF_Num_Plots_"+sector;
    TCanvas* DF_Plot = new TCanvas(pltTitle.c_str(),pltTitle.c_str(),800,600);
    DF_Plot->cd();
    mgDF->Draw("AP");
    mgLeg->Draw("same");

    return DF_Plot;
}

// Makes a plot of the numerator term from the DF through PF calculation
TCanvas* DF_NA_Counts_Legend_Plot( DataSet& AllData, string sector, string targetType ){

    if( !(targetType == "NH3" || targetType == "ND3") ){
	    cout <<"ERROR: enter a valid target type to plot: 'NH3' or 'ND3'\n";
	    TCanvas* emptyC = new TCanvas("emptyC","emptyC",800,600);
	    return emptyC;
    }
 
    vector<int> palette = { 1, 632, 800, 400, 416, 600, 840, 880, 900, 616, 820, 432, 920, 860 };
    int pint = 0;

    TMultiGraph* mgDF = new TMultiGraph();
    TLegend* mgLeg = new TLegend(0.1,0.7,0.4,0.9);
    mgLeg->SetNColumns(3);
    mgLeg->SetHeader("Q^{2} Bins (GeV^{2})","C");

    // Check to see if we're using the average PF across all bins (thruPF = "YES") or if the
    // PF is being calculated for each bin separately (thruPF = else)
    //if( thruPF == "YES" ){
	AllData.CalculateDFThruPF( targetType );
	cout << "Calculating DF thru PF for " << targetType << endl;
    //}
    //else
	//AllData.CalculateDF();

    for(size_t i=0; i<Q2_Bin_Bounds.size()-1; i++){
	double qmid = (Q2_Bin_Bounds[i]+Q2_Bin_Bounds[i+1])/2.0;
	vector<double> xbins, dfvals, zeros, dferrs;
	for(size_t j=0; j<X_Bin_Bounds.size()-1; j++){
	    double xmid = (X_Bin_Bounds[j]+X_Bin_Bounds[j+1])/2.0;
	    double thisDFnacounts = 0; double thisErrDFnacounts = 0;
	    
	    thisDFnacounts = AllData.getNACountsDFThruPF( qmid, xmid, targetType );
	    thisErrDFnacounts = AllData.getNACountsDFThruPFError( qmid, xmid, targetType );
	    
            if( thisDFnacounts != 0.0 ){
		    xbins.push_back(xmid); zeros.push_back(0.0);
		    dfvals.push_back( thisDFnacounts );
		    dferrs.push_back( thisErrDFnacounts );
	    }
	}
	if( xbins.size() > 0 ){
	    TGraphErrors* gr = new TGraphErrors(xbins.size(), xbins.data(), dfvals.data(), zeros.data(), dferrs.data() );
	    gr->SetMarkerColor( palette[pint] );
	    if(pint < 4 ){ gr->SetMarkerStyle(kFullCircle); }
	    else if(pint >= 4 && pint < 8){ gr->SetMarkerStyle(kFullTriangleUp); }
	    else if(pint >= 8){ gr->SetMarkerStyle(kFullSquare); }
	    mgDF->Add( gr, "p" );
	    stringstream legTitle; legTitle << "Q^{2}=" << fixed << setprecision(3) << qmid;
	    mgLeg->AddEntry( gr, legTitle.str().c_str(), "p");
	    pint++;
	}
    }

    string title = targetType+" DF_{"+targetType+"} 'n_{A} Counts' Term for "+sector+"; X; n_{A} Counts Term";
    mgDF->SetTitle(title.c_str());
    mgDF->GetXaxis()->SetLimits(0,0.8); //mgDF->GetYaxis()->SetRangeUser(0.0,0.5);
    
    string pltTitle = targetType+"_DF_NA_Counts_Plots_"+sector;
    TCanvas* DF_Plot = new TCanvas(pltTitle.c_str(),pltTitle.c_str(),800,600);
    DF_Plot->cd();
    mgDF->Draw("AP");
    mgLeg->Draw("same");

    return DF_Plot;
}
*/
//TCanvas* Counts_Q2_Legend_Plot( DataSet& AllData, string sector, string numTarg, string denomTarg ){
TMultiGraph* Counts_Q2_Legend_Plot( DataSet& AllData, string sector, string numTarg, string denomTarg, string Period ){

    if( !( numTarg == "NH3" || numTarg == "ND3" || numTarg == "C" || numTarg == "CH2" || numTarg == "CD2" || numTarg == "F" || numTarg == "ET" || 
	   denomTarg == "NH3" || denomTarg == "ND3" || denomTarg == "C" || denomTarg == "CH2" || denomTarg == "CD2" || denomTarg == "F" || denomTarg == "ET" ) ){
	    cout <<"ERROR: enter a valid target type to plot: 'NH3' or 'ND3'\n";
	    //TCanvas* emptyC = new TCanvas("emptyC","emptyC",800,600);
	    TMultiGraph* emptyC;
	    return emptyC;
    }
 
    vector<int> palette = { 1, 632, 800, 400, 416, 600, 840, 880, 900, 616, 820, 432, 920, 860 };
    int pint = 0;

    TMultiGraph* mgRatio = new TMultiGraph();
    TLegend* mgLeg = new TLegend(0.1,0.7,0.4,0.9);
    mgLeg->SetNColumns(3);
    mgLeg->SetHeader("Q^{2} Bins (GeV^{2})","C");

    for(size_t i=0; i<Q2_Bin_Bounds.size()-1; i++){
	double qmid = (Q2_Bin_Bounds[i]+Q2_Bin_Bounds[i+1])/2.0;
	vector<double> xbins, ratiovals, zeros, ratioerrs;
	for(size_t j=0; j<X_Bin_Bounds.size()-1; j++){
	    double xmid = (X_Bin_Bounds[j]+X_Bin_Bounds[j+1])/2.0;
	    double xavg1 = AllData.getAvgX( numTarg, qmid, xmid );
	    double xavg2 = AllData.getAvgX( denomTarg, qmid, xmid );
	    double xavg = (xavg1 + xavg2)/2.0;
	    double thisRatio = 0; double thisErrRatio = 0;
	    thisRatio = AllData.getCountRatio( qmid, xmid, numTarg, denomTarg, Period );
	    thisErrRatio = AllData.getCountRatioErr( qmid, xmid, numTarg, denomTarg, Period );
            if( thisRatio != 0.0 && thisErrRatio){
		    xbins.push_back(xavg); zeros.push_back(0.0);
		    ratiovals.push_back( thisRatio );
		    ratioerrs.push_back( thisErrRatio );
	    }
	}
	if( xbins.size() > 0 ){
	    TGraphErrors* gr = new TGraphErrors(xbins.size(), xbins.data(), ratiovals.data(), zeros.data(), ratioerrs.data() );
	    gr->SetMarkerColor( palette[pint] );
	    if(pint < 4 ){ gr->SetMarkerStyle(kFullCircle); }
	    else if(pint >= 4 && pint < 8){ gr->SetMarkerStyle(kFullTriangleUp); }
	    else if(pint >= 8){ gr->SetMarkerStyle(kFullSquare); }
	    mgRatio->Add( gr, "p" );
	    stringstream legTitle; legTitle << "Q^{2}=" << fixed << setprecision(3) << qmid;
	    mgLeg->AddEntry( gr, legTitle.str().c_str(), "p");
	    pint++;
	}
    }

    string title = "Ratio of n_{"+numTarg+"}/n_{"+denomTarg+"}; X; n_{"+numTarg+"}/n_{"+denomTarg+"}";
    mgRatio->SetTitle(title.c_str());
    //mgRatio->GetXaxis()->SetLimits(0,0.8); mgRatio->GetYaxis()->SetRangeUser(0.0,0.5);
    //return Ratio_Plot;
    return mgRatio;
}

// This makes a plot of the count ratios for all target types.
TCanvas* Counts_Q2_Legend_All_Plots( DataSet& AllData, string Runperiod, string denomTarg, string Period ){

    if( !( denomTarg == "NH3" || denomTarg == "ND3" || denomTarg == "C" || denomTarg == "CH2" || denomTarg == "CD2" || denomTarg == "F" || denomTarg == "ET" ) ){
	    cout <<"ERROR: enter a valid target type to plot: 'NH3' or 'ND3'\n";
	    TCanvas* emptyC = new TCanvas("emptyC","emptyC",800,600);
	    //TMultiGraph* emptyC;
	    return emptyC;
    }

    string title = Runperiod+"_Data";
    TCanvas* All_Plots = new TCanvas(title.c_str(),title.c_str(),1400,700);
    All_Plots->Divide(3,2);
    vector<string> AllTargs = {"NH3","ND3","C","CH2","CD2","ET","F"};
    int it = 1; // iterator used for plotting
    for( size_t i=0; i<AllTargs.size(); i++ ){
	string numTarg = AllTargs[i];
	if( numTarg != denomTarg ){
	    TMultiGraph* thisRatio = Counts_Q2_Legend_Plot( AllData, Runperiod, numTarg, denomTarg, Period );
	    if( numTarg == "C" || numTarg == "CD2" || numTarg == "CH2" ){
		if( denomTarg == "ND3" ) thisRatio->GetYaxis()->SetRangeUser(0.8,1.2);
		else if( denomTarg == "NH3" ) thisRatio->GetYaxis()->SetRangeUser(1.0,1.4);
	    }
	    All_Plots->cd(it); thisRatio->Draw("AP");
	    it++;
	}
    }
    return All_Plots;
}

void DF_Legend_Plot( DataSet& AllData, string sector, string targetType, vector<double>& DFs, TVirtualPad* pad, bool useScaling=true, bool usePseudoData=true ){
// TCanvas* DF_Legend_Plot( DataSet& AllData, string sector, string targetType, string Period, bool useScaling=true, bool usePseudoData=true ){

    if( !(targetType == "NH3" || targetType == "ND3") ){
	    cout <<"ERROR: enter a valid target type to plot: 'NH3' or 'ND3'\n";
	    return;
    }
 
    vector<int> palette = { 1, 632, 800, 400, 416, 600, 840, 880, 900, 616, 820, 432, 920, 860 };
    int pint = 0;

    TMultiGraph* mgDF = new TMultiGraph();
    TLegend* mgLeg = new TLegend(0.1,0.7,0.5,0.9);
    mgLeg->SetNColumns(3);
    mgLeg->SetHeader("Q^{2} Bins (GeV^{2})","C");

    //AllData.CalculateDF( Period, useScaling, usePseudoData );

    for(size_t i=0; i<Q2_Bin_Bounds.size()-1; i++){
	double qmid = (Q2_Bin_Bounds[i]+Q2_Bin_Bounds[i+1])/2.0;
	// These are for calculating the statistics-weighted
	vector<double> xbins, dfvals, dferrs;
	for(size_t j=0; j<X_Bin_Bounds.size()-1; j++){
	    double xmid = (X_Bin_Bounds[j]+X_Bin_Bounds[j+1])/2.0;
	    double thisDF = 0; 
	    double thisErrDF = 0;
	    if( targetType == "NH3" ){
	        thisDF = AllData.getDF_NH3(qmid, xmid);
	        thisErrDF = AllData.getErrDF_NH3(qmid, xmid);
		if( qmid > 2.7 && qmid < 3 && xmid > 0.2 && xmid < 0.225 ){
		    DFs[0] = thisDF;
		    DFs[1] = thisErrDF;
		}
	    }
	    else if( targetType == "ND3" ){
	        thisDF = AllData.getDF_ND3(qmid, xmid);
	        thisErrDF = AllData.getErrDF_ND3(qmid, xmid);
		if( qmid > 2.7 && qmid < 3 && xmid > 0.2 && xmid < 0.225 ){
		    DFs[0] = thisDF;
		    DFs[1] = thisErrDF;
		}
	    }
	    /*
	    else if( targetType == "NH3" && thruPF == "SYS" ){
	        thisDF = AllData.getDF_NH3(qmid, xmid);
	        thisErrDF = AllData.getTotalErrorDF( targetType, qmid, xmid);
	    }    
	    else if( targetType == "ND3" && thruPF == "SYS" ){
	        thisDF = AllData.getDF_NH3(qmid, xmid);
	        thisErrDF = AllData.getTotalErrorDF( targetType, qmid, xmid);
	    }
	    */
	    double avgX = AllData.getAvgX("All", qmid, xmid ); 
            if( thisDF != 0.0 && avgX != 0.0 ){
		    xbins.push_back(avgX);
		    dfvals.push_back( thisDF );
		    dferrs.push_back( thisErrDF );
	    }
	}
	if( xbins.size() > 0 ){
	    TGraphErrors* gr = new TGraphErrors(xbins.size(), xbins.data(), dfvals.data(), nullptr, dferrs.data() );
	    gr->SetMarkerColor( palette[pint] );
	    if(pint < 4 ){ gr->SetMarkerStyle(kFullCircle); }
	    else if(pint >= 4 && pint < 8){ gr->SetMarkerStyle(kFullTriangleUp); }
	    else if(pint >= 8){ gr->SetMarkerStyle(kFullSquare); }
	    mgDF->Add( gr, "p" );
	    stringstream legTitle; legTitle << "Q^{2}=" << fixed << setprecision(3) << qmid;
	    mgLeg->AddEntry( gr, legTitle.str().c_str(), "p");
	    pint++;
	}
    }

    string title = "DF_{"+targetType+"} for "+sector+"; X; DF_{"+targetType+"}";
    mgDF->SetTitle(title.c_str());
    mgDF->GetXaxis()->SetLimits(0,0.8); 
    if( targetType == "NH3" ) mgDF->GetYaxis()->SetRangeUser(0.1,0.3);
    else mgDF->GetYaxis()->SetRangeUser(0.0,0.5);
    
    //string pltTitle = targetType+"_DF_Plots_"+sector;
    //TCanvas* DF_Plot = new TCanvas(pltTitle.c_str(),pltTitle.c_str(),800,600);
    //DF_Plot->cd();
    pad->cd();
    mgDF->Draw("AP");
    mgLeg->Draw("same");

    //return DF_Plot;
}

// This function assumes that CalculateDF(), which also calculates the PF, has also been called
void PF_Legend_Plot( DataSet& AllData, string sector, string targetType, vector<double>& PFs, TVirtualPad* pad, string BathOrCell="Bath", bool useScaling =true, bool usePseudoData=true,
double ymin=0.4, double ymax=0.7 ){
//TCanvas* PF_Legend_Plot( DataSet& AllData, string sector, string targetType, vector<double>& PFs, string BathOrCell="Bath", bool useScaling =true, bool usePseudoData=true,
//double ymin=0.4, double ymax=0.7 ){


    if( !(targetType == "NH3" || targetType == "ND3") ){
	    cout <<"ERROR: enter a valid target type to plot: 'NH3' or 'ND3'\n";
	    // TCanvas* emptyC = new TCanvas("emptyC","emptyC",800,600);
	    return;
    }

    // This is just used for plotting different colors on the TMultiGraph Objects 
    vector<int> palette = { 1, 632, 800, 400, 416, 600, 840, 880, 900, 616, 820, 432, 920, 860 };
    int pint = 0;

    // These TMultiGraph objects hold the packing fraction information for each x, Q2 bin separately
    // for the "bath" values (which is only what we'll use) and the "cell" values (mostly ignored)
    TMultiGraph* mgPFbath = new TMultiGraph();
    TMultiGraph* mgPFcell = new TMultiGraph();
    TLegend* mgLegbath = new TLegend(0.1,0.7,0.5,0.9);
    TLegend* mgLegcell = new TLegend(0.1,0.7,0.5,0.9);
    mgLegbath->SetNColumns(3);
    mgLegbath->SetHeader("Q^{2} Bins (GeV^{2})","C");
    mgLegcell->SetNColumns(3);
    mgLegcell->SetHeader("Q^{2} Bins (GeV^{2})","C");

    //AllData.CalculatePF( Period, useScaling, usePseudoData ); // Calculates PF in each bin
    //AllData.CalculateDF( Period, useScaling, usePseudoData ); // Calculates PF in each bin

    // These are used to calculate the weighted average of all the PF values in the data set
    double totalPFbath = 0; double denominatorbath = 0;
    double totalPFcell = 0; double denominatorcell = 0;

    // Go over the kinematic bins. Creates a TGraphErrors object for each Q2 bin. The PFs are plotted against
    // x for a given Q2 bin and added to the multigraph.
    for(size_t i=0; i<Q2_Bin_Bounds.size()-1; i++){
	double qmid = (Q2_Bin_Bounds[i]+Q2_Bin_Bounds[i+1])/2.0;
	vector<double> xbinsbath, xbinscell, pfcellvals, pfcellerrs, pfbathvals, pfbatherrs;
	for(size_t j=0; j<X_Bin_Bounds.size()-1; j++){
	    double xmid = (X_Bin_Bounds[j]+X_Bin_Bounds[j+1])/2.0;
	    double thisPFbath = 0; double thisErrPFbath = 0;
	    double thisPFcell = 0; double thisErrPFcell = 0;
	    if( targetType == "NH3" ){
	        thisPFbath = AllData.getPF_bath_NH3(qmid, xmid);
	        thisErrPFbath = AllData.getErrPF_bath_NH3(qmid, xmid);
	        thisPFcell = AllData.getPF_cell_NH3(qmid, xmid);
	        thisErrPFcell = AllData.getErrPF_cell_NH3(qmid, xmid);
	    }
	    else if( targetType == "ND3" ){
	        thisPFbath = AllData.getPF_bath_ND3(qmid, xmid);
	        thisErrPFbath = AllData.getErrPF_bath_ND3(qmid, xmid);
	        thisPFcell = AllData.getPF_cell_ND3(qmid, xmid);
	        thisErrPFcell = AllData.getErrPF_cell_ND3(qmid, xmid);
	    }
            if( thisPFbath > 0.0 && thisErrPFbath > 0.0 ){
		    xbinsbath.push_back(xmid);
		    pfbathvals.push_back( thisPFbath );
		    pfbatherrs.push_back( thisErrPFbath );
		    //if( thisPF > 0.0 && thisErrPF != 0.0 && abs(thisErrPF) < (0.1*thisPF) ){ 
		    totalPFbath += thisPFbath / pow(thisErrPFbath,2);
		    denominatorbath += 1.0 / pow(thisErrPFbath,2);
		    //AllPFBathVals.push_back( thisPFbath );
	    	    //}
	    }
	    if( thisPFcell > 0.0 && thisErrPFcell > 0.0 ){
		    xbinscell.push_back(xmid);
		    pfcellvals.push_back( thisPFcell );
		    pfcellerrs.push_back( thisErrPFcell );
		    //if( thisPF > 0.0 && thisErrPF != 0.0 && abs(thisErrPF) < (0.1*thisPF) ){ 
		    totalPFcell += thisPFcell / pow(thisErrPFcell,2) ;
		    denominatorcell += 1.0 / pow(thisErrPFcell,2);
	    	    //}
	    }	    
	}
	if( xbinsbath.size() > 0 && xbinscell.size() > 0 ){
	    TGraphErrors* gr_cell = new TGraphErrors(xbinscell.size(), xbinscell.data(), pfcellvals.data(), nullptr, pfcellerrs.data() );
	    TGraphErrors* gr_bath = new TGraphErrors(xbinsbath.size(), xbinsbath.data(), pfbathvals.data(), nullptr, pfbatherrs.data() );
	    gr_cell->SetMarkerColor( palette[pint] );
	    gr_bath->SetMarkerColor( palette[pint] );
	    if(pint < 4 ){ gr_cell->SetMarkerStyle(kFullCircle); gr_bath->SetMarkerStyle(kFullCircle); }
	    else if(pint >= 4 && pint < 8){ gr_cell->SetMarkerStyle(kFullTriangleUp); gr_bath->SetMarkerStyle(kFullTriangleUp); }
	    else if(pint >= 8){ gr_cell->SetMarkerStyle(kFullSquare); gr_bath->SetMarkerStyle(kFullSquare); }
	    mgPFcell->Add( gr_cell, "p" );
	    mgPFbath->Add( gr_bath, "p" );
	    //string legTitle = "Q^{2}=" + to_string( trunc(qmid*1000)/1000 );
	    stringstream legTitle; legTitle << "Q^{2}=" << fixed << setprecision(3) << qmid;
	    mgLegcell->AddEntry( gr_cell, legTitle.str().c_str(), "p");
	    mgLegbath->AddEntry( gr_bath, legTitle.str().c_str(), "p");

	    pint++;
	}
    }

    string title = targetType+" PF_{bath} for "+sector+"; X; "+targetType+" PF_{bath}";
    mgPFbath->SetTitle(title.c_str());
    mgPFbath->GetXaxis()->SetLimits(0,0.8); 
    mgPFbath->GetYaxis()->SetRangeUser( ymin, ymax );
    //mgPFbath->GetYaxis()->SetRangeUser( pfbath_min, pfbath_max );
    
    title = targetType+" PF_{cell} for "+sector+"; X; "+targetType+" PF_{cell}";
    mgPFcell->SetTitle(title.c_str());
    mgPFcell->GetXaxis()->SetLimits(0,0.8);
    mgPFcell->GetYaxis()->SetRangeUser( ymin+0.1, ymax+0.1 );

    // Now calculate the weighted average for the bath values
    TFitResultPtr r;
    r = mgPFbath->Fit("pol0","S","Q");
    cout << "PF_bath from Minuit is "<< r->Value(0) <<" +- "<< r->Error(0) << endl;
    cout << "Chi2/NDF from Minuit is "<< r->Chi2()/r->Ndf() << endl;

    // Create a second legend for the reduced chi2
    TLegend* chi2Bath = new TLegend(0.5,0.7,0.79,0.89);
    chi2Bath->SetFillStyle(0);
    chi2Bath->SetLineColorAlpha(0,0.5);
    chi2Bath->SetLineWidth(0);

    stringstream chi2BathTitle; chi2BathTitle << "#chi^{2}_{red}=" << setprecision(3) << r->Chi2()/r->Ndf();
    chi2Bath->AddEntry( chi2BathTitle.str().c_str(), chi2BathTitle.str().c_str(), "");
    chi2Bath->SetTextSize(0.04);

    stringstream avgPFTitle; avgPFTitle << "P_{F}="<< setprecision(4) << r->Value(0) <<" #pm "<< setprecision(2) << r->Error(0)*sqrt( r->Chi2()/r->Ndf() );
    chi2Bath->AddEntry( avgPFTitle.str().c_str(), avgPFTitle.str().c_str(), "");
    chi2Bath->SetTextSize(0.04);

    // Set the PFbath and systematic error for this epoch
    PFs[0] = r->Value(0);
    PFs[1] = r->Error(0)*sqrt( r->Chi2()/r->Ndf() );
    PFs[2] = r->Chi2()/r->Ndf();

    // Now calculate the weighted average for the cell values
    r = mgPFcell->Fit("pol0","S","Q");
    cout << "PF_cell from Minuit is "<< r->Value(0) <<" +- "<< r->Error(0) << endl;
    cout << "Chi2/NDF from Minuit is "<< r->Chi2()/r->Ndf() << endl;

    // Create a second legend for the reduced chi2
    TLegend* chi2Cell = new TLegend(0.5,0.7,0.79,0.89);
    chi2Cell->SetFillStyle(0);
    chi2Cell->SetLineColorAlpha(0,0.5);
    chi2Cell->SetLineWidth(0);

    stringstream chi2CellTitle; chi2CellTitle << "#chi^{2}_{red}=" << setprecision(3) << r->Chi2()/r->Ndf();
    chi2Cell->AddEntry( chi2CellTitle.str().c_str(), chi2CellTitle.str().c_str(), "");
    chi2Cell->SetTextSize(0.04);

    stringstream avgPFcellTitle; avgPFcellTitle << "P_{F}="<< setprecision(4) << r->Value(0) <<" #pm "<< setprecision(2) << r->Error(0)*sqrt( r->Chi2()/r->Ndf() );
    chi2Cell->AddEntry( avgPFcellTitle.str().c_str(), avgPFcellTitle.str().c_str(), "");
    chi2Cell->SetTextSize(0.04);

    // Set the PFcell and systematic error for this epoch
    PFs[3] = r->Value(0);
    PFs[4] = r->Error(0)*sqrt( r->Chi2()/r->Ndf() );
    PFs[5] = r->Chi2()/r->Ndf();

    string pltTitle = targetType+"_PF_Plots_"+sector;
    //TCanvas* PF_Plot = new TCanvas(pltTitle.c_str(),pltTitle.c_str(),800,600);
/*
    PF_Plot->Divide(2,1);
    PF_Plot->cd(1);
    //PF_Plot->cd();
    mgPFbath->Draw("AP");
    mgLegbath->Draw("same");
    chi2Bath->Draw("same");

    PF_Plot->cd(2);
    mgPFcell->Draw("AP");
    mgLegcell->Draw("same");
    chi2Cell->Draw("same");
*/
    //PF_Plot->cd();
    pad->cd();
    if( BathOrCell == "Bath" ){ 
        mgPFbath->Draw("AP");
        mgLegbath->Draw("same");
        chi2Bath->Draw("same");
    }
    else if( BathOrCell == "Cell" ){
        mgPFcell->Draw("AP");
        mgLegcell->Draw("same");
        chi2Cell->Draw("same");
    }

    return;
}

// Plots the ratio of nx/nA counts for every bin in X and Q2 for one epoch
TCanvas* CountRatio_Legend_Plot( DataSet& AllData, string sector, string targNum, string targDen, string Period ){

    vector<int> palette = { 1, 632, 800, 400, 416, 600, 840, 880, 900, 616, 820, 432, 920, 860 };
    int pint = 0;

    TMultiGraph* mgRatio = new TMultiGraph();
    TLegend* mgLeg = new TLegend(0.1,0.7,0.5,0.9);
    mgLeg->SetNColumns(3);
    mgLeg->SetHeader("Q^{2} Bins (GeV^{2})","C");

    for(int i=0; i<Q2_Bin_Bounds.size()-1; i++){
	double qmid = (Q2_Bin_Bounds[i] + Q2_Bin_Bounds[i+1])/2.0;
	vector<double> ratioCounts, ratioCountsErr, zeros, xVals;
	for(int j=0; j<X_Bin_Bounds.size()-1; j++){
	    double xmid = (X_Bin_Bounds[j] + X_Bin_Bounds[j+1])/2.0;
	    double ratio = AllData.getCountRatio( qmid, xmid, targNum, targDen, Period );
	    double ratioErr = AllData.getCountRatioErr( qmid, xmid, targNum, targDen, Period );
	    if( ratio > 0 ){
		ratioCounts.push_back( ratio ); 
		ratioCountsErr.push_back( ratioErr );
		zeros.push_back( 0 );
		xVals.push_back( xmid );
	    }
	}
	if( ratioCounts.size() > 0 ){
	    TGraphErrors* gr = new TGraphErrors(xVals.size(), xVals.data(), ratioCounts.data(), zeros.data(), ratioCountsErr.data() );
	    gr->SetMarkerStyle(kFullCircle);
	    gr->SetMarkerColor( palette[i] );
	    mgRatio->Add( gr, "p" );
	    stringstream q2label; q2label << Form("Q^{2} = %.2f", qmid);
	    mgLeg->AddEntry( gr, q2label.str().c_str(), "p" );
	}
    }

    // Set the Y-axis bounds based on target type
    double ymin = 0; double ymax = 1;
    if( (targNum == "C" || targNum == "CH2") && (targDen == "NH3" || targDen == "ND3") ){ ymin = 1; ymax = 1.4; }
    else if( targNum == "ET" && (targDen == "NH3" || targDen == "ND3") ){ ymin = 0.1; ymax = 0.4; }
    else if( targNum == "F"  && (targDen == "NH3" || targDen == "ND3") ){ ymin = 0;   ymax = 0.02; }
    
    // Edge cases for the background ratio plots
    if( targNum == "CH2" && targDen == "C" ){ ymin = 0.85; ymax = 1.2; }
    else if( targNum == "ET" && targDen == "C"  ){ ymin = 0.1;  ymax = 0.4; }
    else if( targNum == "F"  && targDen == "C"  ){ ymin =0.002; ymax =0.015;}
    else if( targNum == "F"  && targDen == "CH2"){ ymin =0.002; ymax =0.015;}
    else if( targNum == "ET" && targDen == "CH2"){ ymin = 0.15; ymax = 0.4; }
    else if( targNum == "F"  && targDen == "ET" ){ ymin = 0.02; ymax = 0.05;}

    string mgTitle = "Count Ratio of "+targNum+"/"+targDen+" for "+sector+"; X;"+targNum+"/"+targDen;
    mgRatio->SetTitle( mgTitle.c_str() );
    mgRatio->GetYaxis()->SetRangeUser(ymin, ymax);

    string pltTitle = "CR_"+targNum+"_"+targDen+"_Plots_"+sector;
    TCanvas* CR_Plot = new TCanvas(pltTitle.c_str(),pltTitle.c_str(),800,600);
    CR_Plot->cd();
    mgRatio->Draw("AP");
    mgLeg->Draw("same");
    return CR_Plot;

}

// Plots the raw double-spin asymmetries for NH3 targets
TCanvas* AllRaw_Legend_Plot( DataSet& AllData, string sector, string targetType, string FCcorr="No" ){

    vector<int> palette = { 1, 632, 800, 400, 416, 600, 840, 880, 900, 616, 820, 432, 920, 860 };
    int pint = 0;

    if( !( targetType == "NH3" || targetType == "ND3" ) ){
	cout << "Valid target types are 'NH3' or 'ND3'\n";
	TCanvas* empty = new TCanvas();
	return empty;
    }

    TMultiGraph* mgAllRaw = new TMultiGraph();
    TLegend* mgLeg = new TLegend(0.1,0.7,0.4,0.9);
    mgLeg->SetNColumns(3);
    mgLeg->SetHeader("Q^{2} Bins (GeV^{2})","C");

    AllData.CalculateAllRaw(); // Calculate the raw double-spin asymmetries for NH3

    for(size_t i=0; i<Q2_Bin_Bounds.size()-1; i++){
	double qmid = (Q2_Bin_Bounds[i]+Q2_Bin_Bounds[i+1])/2.0;
	vector<double> xbins, arawvals, zeros, arawerrs;
	for(size_t j=0; j<X_Bin_Bounds.size()-1; j++){
	    double xmid = (X_Bin_Bounds[j]+X_Bin_Bounds[j+1])/2.0;
	    double thisAllRaw = 0.0; double thisErrAllRaw = 0.0;
	    if( FCcorr == "Yes" ){
	        thisAllRaw = AllData.getAllRawNoFC(qmid, xmid, targetType ); 
	        thisErrAllRaw = AllData.getAllRawNoFCErr(qmid, xmid, targetType);//0.05*thisAllRaw; // Set this to 5% error as a placeholder
	    }
	    else{
	        thisAllRaw = AllData.getAllRaw(qmid, xmid, targetType ); 
	        thisErrAllRaw = AllData.getAllRawErr(qmid, xmid, targetType);//0.05*thisAllRaw; // Set this to 5% error as a placeholder
	    }
            if( thisAllRaw != 0.0 ){
		    xbins.push_back(xmid); zeros.push_back(0.0);
		    arawvals.push_back( thisAllRaw );
		    arawerrs.push_back( thisErrAllRaw );
	    }
	}
	if( xbins.size() > 0 ){
	    TGraphErrors* gr = new TGraphErrors(xbins.size(), xbins.data(), arawvals.data(), zeros.data(), arawerrs.data() );
	    gr->SetMarkerColor( palette[pint] );
	    if(pint < 4 ){ gr->SetMarkerStyle(kFullCircle); }
	    else if(pint >= 4 && pint < 8){ gr->SetMarkerStyle(kFullTriangleUp); }
	    else if(pint >= 8){ gr->SetMarkerStyle(kFullSquare); }
	    mgAllRaw->Add( gr, "p" );
	    stringstream legTitle; legTitle << "Q^{2}=" << fixed << setprecision(3) << qmid;
	    mgLeg->AddEntry( gr, legTitle.str().c_str(), "p");
	    pint++;
	}
    }

    string title = targetType+" A_{||,raw} for "+sector+"; X; A_{||}";
    mgAllRaw->SetTitle(title.c_str());
    mgAllRaw->GetXaxis()->SetLimits(0,0.8); mgAllRaw->GetYaxis()->SetRangeUser(-0.07,0.07);
    
    string pltTitle = "AllRaw"+targetType+" Plots "+sector;
    TCanvas* AllRaw_Plot = new TCanvas(pltTitle.c_str(),pltTitle.c_str(),800,600);
    AllRaw_Plot->cd();
    //gPad->SetLogx();
    mgAllRaw->Draw("AP");
    mgLeg->Draw("same");

    return AllRaw_Plot;
}

// Plots the raw double-spin asymmetries for NH3 targets
TCanvas* AllRawNH3_Legend_Plot( DataSet& AllData, string sector ){

    vector<int> palette = { 1, 632, 800, 400, 416, 600, 840, 880, 900, 616, 820, 432, 920, 860 };
    int pint = 0;

    TMultiGraph* mgAllRaw = new TMultiGraph();
    TLegend* mgLeg = new TLegend(0.1,0.7,0.4,0.9);
    mgLeg->SetNColumns(3);
    mgLeg->SetHeader("Q^{2} Bins (GeV^{2})","C");

    AllData.CalculateAllRaw(); // Calculate the raw double-spin asymmetries for NH3

    for(size_t i=0; i<Q2_Bin_Bounds.size()-1; i++){
	double qmid = (Q2_Bin_Bounds[i]+Q2_Bin_Bounds[i+1])/2.0;
	vector<double> xbins, arawvals, zeros, arawerrs;
	for(size_t j=0; j<X_Bin_Bounds.size()-1; j++){
	    double xmid = (X_Bin_Bounds[j]+X_Bin_Bounds[j+1])/2.0;
	    double thisAllRaw = AllData.getAllRaw(qmid, xmid, "NH3" ); 
	    double thisErrAllRaw = AllData.getAllRawErr(qmid, xmid, "NH3");//0.05*thisAllRaw; // Set this to 5% error as a placeholder
            if( thisAllRaw != 0.0 ){
		    xbins.push_back(xmid); zeros.push_back(0.0);
		    arawvals.push_back( thisAllRaw );
		    arawerrs.push_back( thisErrAllRaw );
	    }
	}
	if( xbins.size() > 0 ){
	    TGraphErrors* gr = new TGraphErrors(xbins.size(), xbins.data(), arawvals.data(), zeros.data(), arawerrs.data() );
	    gr->SetMarkerColor( palette[pint] );
	    if(pint < 4 ){ gr->SetMarkerStyle(kFullCircle); }
	    else if(pint >= 4 && pint < 8){ gr->SetMarkerStyle(kFullTriangleUp); }
	    else if(pint >= 8){ gr->SetMarkerStyle(kFullSquare); }
	    mgAllRaw->Add( gr, "p" );
	    stringstream legTitle; legTitle << "Q^{2}=" << fixed << setprecision(3) << qmid;
	    mgLeg->AddEntry( gr, legTitle.str().c_str(), "p");
	    pint++;
	}
    }

    string title = "NH3 A_{||,raw} for "+sector+"; X; A_{||}";
    mgAllRaw->SetTitle(title.c_str());
    mgAllRaw->GetXaxis()->SetLimits(0,0.8); mgAllRaw->GetYaxis()->SetRangeUser(-0.07,0.07);
    
    string pltTitle = "AllRawNH3 Plots "+sector;
    TCanvas* AllRawNH3_Plot = new TCanvas(pltTitle.c_str(),pltTitle.c_str(),800,600);
    AllRawNH3_Plot->cd();
    //gPad->SetLogx();
    mgAllRaw->Draw("AP");
    mgLeg->Draw("same");

    return AllRawNH3_Plot;
}
// Plots the raw double-spin asymmetries for ND3 targets
TCanvas* AllRawND3_Legend_Plot( DataSet& AllData, string sector ){

    vector<int> palette = { 1, 632, 800, 400, 416, 600, 840, 880, 900, 616, 820, 432, 920, 860 };
    int pint = 0;

    TMultiGraph* mgAllRaw = new TMultiGraph();
    TLegend* mgLeg = new TLegend(0.1,0.7,0.4,0.9);
    mgLeg->SetNColumns(3);
    mgLeg->SetHeader("Q^{2} Bins (GeV^{2})","C");

    AllData.CalculateAllRaw(); // Calculate the raw double-spin asymmetries for ND3

    for(size_t i=0; i<Q2_Bin_Bounds.size()-1; i++){
	double qmid = (Q2_Bin_Bounds[i]+Q2_Bin_Bounds[i+1])/2.0;
	vector<double> xbins, arawvals, zeros, arawerrs;
	for(size_t j=0; j<X_Bin_Bounds.size()-1; j++){
	    double xmid = (X_Bin_Bounds[j]+X_Bin_Bounds[j+1])/2.0;
	    double thisAllRaw = AllData.getAllRaw(qmid, xmid, "ND3" ); 
	    double thisErrAllRaw = AllData.getAllRawErr(qmid, xmid, "ND3");//0.05*thisAllRaw; // Set this to 5% error as a placeholder
            if( thisAllRaw != 0.0 ){
		    xbins.push_back(xmid); zeros.push_back(0.0);
		    arawvals.push_back( thisAllRaw );
		    arawerrs.push_back( thisErrAllRaw );
	    }
	}
	if( xbins.size() > 0 ){
	    TGraphErrors* gr = new TGraphErrors(xbins.size(), xbins.data(), arawvals.data(), zeros.data(), arawerrs.data() );
	    gr->SetMarkerColor( palette[pint] );
	    if(pint < 4 ){ gr->SetMarkerStyle(kFullCircle); }
	    else if(pint >= 4 && pint < 8){ gr->SetMarkerStyle(kFullTriangleUp); }
	    else if(pint >= 8){ gr->SetMarkerStyle(kFullSquare); }
	    mgAllRaw->Add( gr, "p" );
	    stringstream legTitle; legTitle << "Q^{2}=" << fixed << setprecision(3) << qmid;
	    mgLeg->AddEntry( gr, legTitle.str().c_str(), "p");
	    pint++;
	}
    }

    string title = "ND3 A_{||,raw} for "+sector+"; X; A_{||}";
    mgAllRaw->SetTitle(title.c_str());
    mgAllRaw->GetXaxis()->SetLimits(0,0.8); mgAllRaw->GetYaxis()->SetRangeUser(-0.07,0.07);
    
    string pltTitle = "AllRawND3 Plots "+sector;
    TCanvas* AllRawND3_Plot = new TCanvas(pltTitle.c_str(),pltTitle.c_str(),800,600);
    AllRawND3_Plot->cd();
    //gPad->SetLogx();
    mgAllRaw->Draw("AP");
    mgLeg->Draw("same");

    return AllRawND3_Plot;
}

// Plots the raw double-spin asymmetries for NH3 targets
vector<TCanvas*> AllRawNH3_Q2_Bin_Plot( DataSet& AllData, string sector ){
    
    string title1 = sector + " 1";
    string title2 = sector + " 2";
    TCanvas* q2BinPlots1 = new TCanvas(title1.c_str(),title1.c_str(),1900,1000); // Plots the first 15 x bins
    TCanvas* q2BinPlots2 = new TCanvas(title2.c_str(),title2.c_str(),1900,1000); // Plots the last 15 x bins
    q2BinPlots1->Divide(5,3,0,0); q2BinPlots2->Divide(5,3,0,0);
    // Iterators used for drawing to the correct canvas
    int plotIt = 1;

    AllData.CalculateAllRaw(); // Calculate the raw double-spin asymmetries for NH3

    for(size_t i=0; i<X_Bin_Bounds.size()-1; i++){
	double xmid = (X_Bin_Bounds[i]+X_Bin_Bounds[i+1])/2.0;
	vector<double> q2bins, arawvals, zeros, arawerrs;
	for(size_t j=0; j<Q2_Bin_Bounds.size()-1; j++){
	    double q2mid = (Q2_Bin_Bounds[j]+Q2_Bin_Bounds[j+1])/2.0;
	    double thisAllRaw = AllData.getAllRaw(q2mid, xmid, "NH3" ); 
	    double thisErrAllRaw = AllData.getAllRawErr(q2mid, xmid, "NH3");//0.05*thisAllRaw; // Set this to 5% error as a placeholder
            if( thisAllRaw != 0.0 ){
		    q2bins.push_back(q2mid); zeros.push_back(0.0);
		    arawvals.push_back( thisAllRaw );
		    arawerrs.push_back( thisErrAllRaw );
	    }
	}
	if( q2bins.size() > 0 ){
	    TGraphErrors* gr = new TGraphErrors(q2bins.size(), q2bins.data(), arawvals.data(), zeros.data(), arawerrs.data() ); gr->SetTitle("");
	    gr->SetMarkerStyle(kFullCircle);
	    gr->GetXaxis()->SetLimits(0,11); // 0 to 11 in GeV^2
	    gr->GetYaxis()->SetRangeUser(0,0.2);// Raw asymmetry from 0 to 5%
	    string title = "A_{||} at X = "+to_string(xmid);
	    TLegend* leg = new TLegend(0.1,0.5,0.7,0.9);
	    leg->AddEntry(gr, title.c_str(), ""); leg->SetFillStyle(0); leg->SetBorderSize(0);
	    if( plotIt > 15 ){
		q2BinPlots2->cd(plotIt-15); //gr->SetTitle(title.c_str());
		//gPad->SetLogx();
		gr->Draw("AP"); leg->Draw("same"); plotIt++;
	    }		
	    else{
		q2BinPlots1->cd(plotIt); //gr->SetTitle(title.c_str());
		//gPad->SetLogx();
		gr->Draw("AP"); leg->Draw("same"); plotIt++;
	    }
	}
    }

    vector<TCanvas*> Plots = { q2BinPlots1, q2BinPlots2 };
    return Plots;
}

// Makes a plot of the PF for the target cell only, excluding the LHe bath
/*
TCanvas* PFcell_Legend_Plot( DataSet& AllData, string sector ){
 
    vector<int> palette = { 1, 632, 800, 400, 416, 600, 840, 880, 900, 616, 820, 432, 920, 860 };
    int pint = 0;

    TMultiGraph* mgPF = new TMultiGraph();
    TLegend* mgLeg = new TLegend(0.1,0.7,0.4,0.9);
    mgLeg->SetNColumns(3);
    mgLeg->SetHeader("Q^{2} Bins (GeV^{2})","C");
    AllData.CalculatePF();

    // These are used to calculate the weighted average of all the PF values in the data set
    double totalPF = 0; double denominator = 0;

    for(size_t i=0; i<Q2_Bin_Bounds.size()-1; i++){
	double qmid = (Q2_Bin_Bounds[i]+Q2_Bin_Bounds[i+1])/2.0;
	vector<double> xbins, dfnh3vals, zeros, dfnh3errs;
	for(size_t j=0; j<X_Bin_Bounds.size()-1; j++){
	    double xmid = (X_Bin_Bounds[j]+X_Bin_Bounds[j+1])/2.0;
	    // Multiplying by LHe / 5 cm gives the PF for the cell alone
	    double thisPF = (LHe/5.0)*AllData.getPF_bath_NH3(qmid, xmid);
	    double thisErrPF = (LHe/5.0)*AllData.getErrPF_bath_NH3(qmid, xmid);
	    if( thisPF != 0.0 ){
		xbins.push_back(xmid); zeros.push_back(0.0);
		dfnh3vals.push_back( thisPF );
		dfnh3errs.push_back( thisErrPF );
		if( thisPF > 0.0 && thisErrPF != 0.0 && abs(thisErrPF) < (0.1*thisPF) ){ 
		    totalPF += thisPF / pow(thisErrPF,2) ;
		    denominator += 1.0 / pow(thisErrPF,2);
	    	}
	    }
	}
	if( xbins.size() > 0 ){
	    TGraphErrors* gr = new TGraphErrors(xbins.size(), xbins.data(), dfnh3vals.data(), zeros.data(), dfnh3errs.data() );
	    gr->SetMarkerColor( palette[pint] );
	    if(pint < 4 ) gr->SetMarkerStyle(kFullCircle); 
	    else if(pint >= 4 && pint < 8) gr->SetMarkerStyle(kFullTriangleUp);
	    else if(pint >= 8) gr->SetMarkerStyle(kFullSquare);
	    mgPF->Add( gr, "p" ); 
	    //string legTitle = "Q^{2}=" + to_string( trunc(qmid*1000)/1000 );
	    stringstream legTitle; legTitle << "Q^{2}=" << fixed << setprecision(3) << qmid;
	    mgLeg->AddEntry( gr, legTitle.str().c_str(), "p");
	    pint++;
	}
    }
    if( denominator != 0.0 ){
        double avgPF = totalPF / denominator;
        cout << "Average PF_cell for data set is: "<< avgPF <<" +/- "<< 1.0/sqrt(denominator) << endl;
    }
    string title = "NH3 PF_{cell} for "+sector+"; X; NH3 PF_{cell}";
    mgPF->SetTitle(title.c_str());
    mgPF->GetXaxis()->SetLimits(0,0.8); mgPF->GetYaxis()->SetRangeUser(0.0,1.0);

    string pltTitle = "PF_{cell} Plot "+sector;
    TCanvas* PF_Plot = new TCanvas(pltTitle.c_str(),pltTitle.c_str(),800,600);
    PF_Plot->cd();
    mgPF->Draw("AP");
    mgLeg->Draw("same");
    return PF_Plot;
}
*/


// This function creates a single plot for a given epoch, creating DF values that are averaged over
// all Q^2 bins for a fixed bin in Bjorken x. Assumes DF and PF have already been calculated.
//TCanvas* DF_Q2_Averaged_Plot( vector<DataSet>& Epochs, string epoch, string targetType, vector<string> Legends={}, string thruPF="NO", double ymin=0.1, double ymax=0.3 ){
void DF_Q2_Averaged_Plot( vector<DataSet>& Epochs, string epoch, string targetType, TVirtualPad* pad, double ymin=0.1, double ymax=0.3 ){

    vector<int> palette = { 1, 632, 800, 416, 600, 840, 880, 900, 616, 820, 432, 920, 860,
   			   40,  41,  42,  46,  28,  30,  49 };
    int pint = 0;

    if( !(targetType == "NH3" || targetType == "ND3") ){
	    cout <<"ERROR: enter a valid target type to plot: 'NH3' or 'ND3'\n";
	    return;
    }

    // Now read in the data
 TLegend* leg = new TLegend(0.1,0.65,0.7,0.9);
 //leg->SetHeader("Epochs","C"); 
 leg->SetNColumns(3);
 TMultiGraph* mg = new TMultiGraph();

 string mgTitle = "All RG-C "+targetType+" Data; X; D_{F}";
 if( epoch == "Summer" ) mgTitle = epoch+" 2022 "+ targetType+" Data; X; D_{F}";
 else if( epoch == "Fall" ) mgTitle = epoch+" 2022 "+ targetType+" Data; X; D_{F}";
 else if( epoch == "Spring" ) mgTitle = epoch+" 2023 "+ targetType+" Data; X; D_{F}";
 else mgTitle = epoch+" D_{F} Variations; X; D_{F}";

 mg->SetTitle(mgTitle.c_str());
 vector<TGraphErrors*> DF_Epoch_Plots;
 // First, loop over all the epochs and calculate the weighted average DF across all Q2 bins for a given X bin
 for(size_t k=0; k<Epochs.size(); k++){

    DataSet AllData = Epochs[k]; // Very inefficient!!! Needs improvement!!!!11!
    vector<double> xbins, zeros, avg_dfs, avg_df_errs;

    for( size_t i=0; i<X_Bin_Bounds.size()-1; i++ ){

	double xmid = (X_Bin_Bounds[i]+X_Bin_Bounds[i+1])/2.0;
	vector<double> q2bins, dfvals, dferrs, modelq2bins, modeldfvals, modeldferrs, xavgs;
	for( size_t j=0; j<Q2_Bin_Bounds.size()-1; j++ ){

	    double qmid = (Q2_Bin_Bounds[j]+Q2_Bin_Bounds[j+1])/2.0;
	    double xavg = AllData.getAvgX("All", qmid, xmid );
	    double thisDF = 0; double thisErrDF = 0;

	    if( targetType == "NH3" ){
	        thisDF = AllData.getDF_NH3(qmid, xmid);
	        thisErrDF = AllData.getErrDF_NH3(qmid, xmid);
	    }
	    else if( targetType == "ND3" ){
	        thisDF = AllData.getDF_ND3(qmid, xmid);
	        thisErrDF = AllData.getErrDF_ND3(qmid, xmid);
	    }
            if( thisDF > 0.0 && thisErrDF > 0.0 ){
		q2bins.push_back(qmid);
		dfvals.push_back( thisDF );
		dferrs.push_back( thisErrDF );
		xavgs.push_back( xavg );
	    }
	}
	if( q2bins.size() > 0 ){
	    // Calculate the weighted average of all the DFs in the given xbin
	    double num = 0.0; double den = 0.0;
	    double xavgSum = 0.0;
	    //FIXME: This doesn't plot the values at the true stat.-weighted mean of the X values, it
	    //       just uses the regular mean. This should be ok since it's just meant to show the
	    //       overall behavior...
	    for(size_t z=0; z<q2bins.size(); z++){
		num += dfvals[z] / pow(dferrs[z],2);
		den += 1.0 / pow(dferrs[z],2);
		xavgSum += xavgs[z];
	    }
	    if( den != 0.0 ){
		//xbins.push_back(xmid); 
		xbins.push_back( xavgSum / (xavgs.size()*1.0) );
		zeros.push_back(0.0);
		avg_dfs.push_back( num / den ); avg_df_errs.push_back( 1.0 / sqrt(den) );
	    }
	}
    }
    //cout << "Epoch "<< epoch <<" first bin value: "<< avg_dfs[0] << endl;
    TGraphErrors* thisGR = new TGraphErrors(xbins.size(), xbins.data(), avg_dfs.data(), zeros.data(), avg_df_errs.data() );
    DF_Epoch_Plots.push_back(thisGR);
 }
 // If I'm looking at the special case of the summer, fall, and spring data sets, change the palette vector
 if( DF_Epoch_Plots.size() == 3 ) palette = { 4, 3, 2 }; // blue, green, red

    for( size_t i=0; i<DF_Epoch_Plots.size(); i++){
        DF_Epoch_Plots[i]->SetMarkerColor(palette[i]); 
        DF_Epoch_Plots[i]->SetMarkerSize(1);
        if( i<2 ) DF_Epoch_Plots[i]->SetMarkerStyle(kFullCircle);
        else if( i<4 ) DF_Epoch_Plots[i]->SetMarkerStyle(kFullTriangleUp);
        else if( i<6 ) DF_Epoch_Plots[i]->SetMarkerStyle(kFullTriangleDown);
        else if( i<10) DF_Epoch_Plots[i]->SetMarkerStyle(kFullSquare);
	else if( i<13) DF_Epoch_Plots[i]->SetMarkerStyle(33); // full diamond
	else if( i<16) DF_Epoch_Plots[i]->SetMarkerStyle(34); // full cross
	else DF_Epoch_Plots[i]->SetMarkerStyle(29); // full star
    }

    for(size_t i=0; i<DF_Epoch_Plots.size(); i++){
	string thisTitle = "Data Set "+to_string(i+1);
	if( targetType == "NH3" && epoch == "Summer" ) thisTitle = "Epoch "+to_string(i+2);
	else if( targetType == "NH3" && epoch == "Fall" ) thisTitle = "Epoch "+to_string(i+11);
	else if( targetType == "NH3" && epoch == "Spring" ) thisTitle = "Epoch "+to_string(i+20);
	
	else if( targetType == "ND3" && epoch == "Summer" ) thisTitle = "Epoch "+to_string(i+1);
	else if( targetType == "ND3" && epoch == "Fall" ) thisTitle = "Epoch "+to_string(i+5);
	else if( targetType == "ND3" && epoch == "Spring" ) thisTitle = "Epoch "+to_string(i+11);

	else{
	    //thisTitle = Legends[i];
	}

        //DF_Epoch_Plots[i]->SetMarkerColor( palette[i+1] );
        //DF_Epoch_Plots[i]->SetMarkerStyle(i);
        mg->Add( DF_Epoch_Plots[i] , "p");
        leg->AddEntry( DF_Epoch_Plots[i] , thisTitle.c_str(), "p");
    }
    if( targetType == "NH3" ) mg->GetYaxis()->SetRangeUser(ymin,ymax);
    else if( targetType == "ND3" ) mg->GetYaxis()->SetRangeUser(ymin,ymax);
    mg->GetXaxis()->SetLimits(0.12,0.8);
    epoch += "_DFs";

    //TCanvas* mgplt = new TCanvas(epoch.c_str(),epoch.c_str(),800,600);
    //mgplt->cd();
    pad->cd();
    mg->Draw("ap"); leg->Draw("same");

    return;

}

// This function creates a single plot for a set of given epochs, creating PF values that are averaged over
// all Q^2 bins for a fixed bin in Bjorken x. Assumes DF and PF have already been calculated.
//void DF_Q2_Averaged_Plot( vector<DataSet>& Epochs, string epoch, string targetType, TVirtualPad* pad, double ymin=0.1, double ymax=0.3 ){
//TCanvas* PF_Q2_Averaged_Plot( vector<DataSet>& Epochs, string epoch, string targetType, vector<string> Legends={}, double ymin=0.1, double ymax=0.3 ){
void PF_Q2_Averaged_Plot( vector<DataSet>& Epochs, string epoch, string targetType, string BathOrCell, TVirtualPad* pad, double ymin=0.4, double ymax=0.7 ){

    vector<int> palette = { 1, 632, 800, 416, 600, 840, 880, 900, 616, 820, 432, 920, 860,
   			   40,  41,  42,  46,  28,  30,  49 };
    int pint = 0;

    if( !(targetType == "NH3" || targetType == "ND3") ){
	cout <<"ERROR: enter a valid target type to plot: 'NH3' or 'ND3'\n";
	return;
    }
    if( !(BathOrCell == "Bath" || BathOrCell == "Cell") ){
	cout <<"ERROR: You must choose either 'Bath' or 'Cell' as the draw options.\n";
	return;
    }

    // Now read in the data
 double yLegMin = 0.65; double yLegMax = 0.9;
 if( BathOrCell == "Cell" ){ yLegMin = 0.1; yLegMax = 0.35; }

 TLegend* leg = new TLegend(0.1,yLegMin,0.7,yLegMax);
 //leg->SetHeader("Epochs","C"); 
 leg->SetNColumns(3);
 TMultiGraph* mg = new TMultiGraph();

 string mgTitle = "All RG-C "+targetType+" Data; X; PF_{"+ BathOrCell +"}";
 if( epoch == "Summer" ) mgTitle = epoch+" 2022 "+ targetType+" Data; X; PF_{"+ BathOrCell +"}";
 else if( epoch == "Fall" ) mgTitle = epoch+" 2022 "+ targetType+" Data; X; PF_{"+ BathOrCell +"}";
 else if( epoch == "Spring" ) mgTitle = epoch+" 2023 "+ targetType+" Data; X; PF_{"+ BathOrCell +"}";
 else mgTitle = epoch+" PF_{"+ BathOrCell +"} Variations; X; PF_{"+ BathOrCell +"}";

 mg->SetTitle(mgTitle.c_str());
 vector<TGraphErrors*> PF_Epoch_Plots;
 for(size_t k=0; k<Epochs.size(); k++){

    DataSet AllData = Epochs[k]; // Very inefficient!!! Needs improvement!!!!11!
    vector<double> xbins, zeros, avg_pfs, avg_pf_errs;

    for(size_t i=0; i<X_Bin_Bounds.size()-1; i++){
	double xmid = (X_Bin_Bounds[i]+X_Bin_Bounds[i+1])/2.0;
	vector<double> q2bins, pfvals, pferrs, xavgs;
	for(size_t j=0; j<Q2_Bin_Bounds.size()-1; j++){
	    double qmid = (Q2_Bin_Bounds[j]+Q2_Bin_Bounds[j+1])/2.0;
	    double xavg = AllData.getAvgX("All", qmid, xmid );
	    double thisPF = 0; double thisErrPF = 0;
	    // double modelDF = ModelDFs.getDF_FixedPF_NH3( qmid, xmid );
	    // double modelDFerr = 0;
	    if( targetType == "NH3" ){
		if( BathOrCell == "Bath" ){
	            thisPF = AllData.getPF_bath_NH3(qmid, xmid);
	            thisErrPF = AllData.getErrPF_bath_NH3(qmid, xmid);
		}
		else if( BathOrCell == "Cell" ){
	            thisPF = AllData.getPF_cell_NH3(qmid, xmid);
	            thisErrPF = AllData.getErrPF_cell_NH3(qmid, xmid);
		}
	    }
	    else if( targetType == "ND3" ){
		if( BathOrCell == "Bath" ){
	            thisPF = AllData.getPF_bath_ND3(qmid, xmid);
	            thisErrPF = AllData.getErrPF_bath_ND3(qmid, xmid);
		}
		else if( BathOrCell == "Cell" ){
	            thisPF = AllData.getPF_cell_ND3(qmid, xmid);
	            thisErrPF = AllData.getErrPF_cell_ND3(qmid, xmid);
		}
	    }
            if( thisPF > 0.0 && thisErrPF > 0.0 ){
		q2bins.push_back(qmid);
		pfvals.push_back( thisPF );
		pferrs.push_back( thisErrPF );
		xavgs.push_back( xavg );
	    }
	}
	if( q2bins.size() > 0 ){
	    // Calculate the weighted average of all the PFs in the given xbin
	    double num = 0.0; double den = 0.0;
	    double xavgSum = 0.0;
	    for(size_t z=0; z<q2bins.size(); z++){
		num += pfvals[z] / pow(pferrs[z],2);
		den += 1.0 / pow(pferrs[z],2);
		xavgSum += xavgs[z];
	    }
	    if( den != 0.0 ){
		//xbins.push_back(xmid); 
		xbins.push_back( xavgSum / (xavgs.size()*1.0) );
		zeros.push_back(0.0);
		avg_pfs.push_back( num / den ); avg_pf_errs.push_back( 1.0 / sqrt(den) );
	    }
	}
    }
    //cout << "Epoch "<< epoch <<" first bin value: "<< avg_pfs[0] << endl;
    TGraphErrors* thisGR = new TGraphErrors(xbins.size(), xbins.data(), avg_pfs.data(), zeros.data(), avg_pf_errs.data() );
    PF_Epoch_Plots.push_back(thisGR);
 }
 // If I'm looking at the special case of the summer, fall, and spring data sets, change the palette vector
 if( PF_Epoch_Plots.size() == 3 ) palette = { 4, 3, 2 }; // blue, green, red

    for( size_t i=0; i<PF_Epoch_Plots.size(); i++){
        PF_Epoch_Plots[i]->SetMarkerColor(palette[i]); 
        PF_Epoch_Plots[i]->SetMarkerSize(1);
        if( i<2 ) PF_Epoch_Plots[i]->SetMarkerStyle(kFullCircle);
        else if( i<4 ) PF_Epoch_Plots[i]->SetMarkerStyle(kFullTriangleUp);
        else if( i<6 ) PF_Epoch_Plots[i]->SetMarkerStyle(kFullTriangleDown);
        else if( i<10) PF_Epoch_Plots[i]->SetMarkerStyle(kFullSquare);
	else if( i<13) PF_Epoch_Plots[i]->SetMarkerStyle(33); // full diamond
	else if( i<16) PF_Epoch_Plots[i]->SetMarkerStyle(34); // full cross
	else PF_Epoch_Plots[i]->SetMarkerStyle(29); // full star

    }

    for(size_t i=0; i<PF_Epoch_Plots.size(); i++){
	string thisTitle = "Data Set "+to_string(i+1);
	if( targetType == "NH3" && epoch == "Summer" ) thisTitle = "Epoch "+to_string(i+2);
	else if( targetType == "NH3" && epoch == "Fall" ) thisTitle = "Epoch "+to_string(i+11);
	else if( targetType == "NH3" && epoch == "Spring" ) thisTitle = "Epoch "+to_string(i+20);
	
	else if( targetType == "ND3" && epoch == "Summer" ) thisTitle = "Epoch "+to_string(i+1);
	else if( targetType == "ND3" && epoch == "Fall" ) thisTitle = "Epoch "+to_string(i+5);
	else if( targetType == "ND3" && epoch == "Spring" ) thisTitle = "Epoch "+to_string(i+11);

	else{
	    //thisTitle = Legends[i];
	}

        //PF_Epoch_Plots[i]->SetMarkerColor( palette[i+1] );
        //PF_Epoch_Plots[i]->SetMarkerStyle(i);
        mg->Add( PF_Epoch_Plots[i] , "p");
        leg->AddEntry( PF_Epoch_Plots[i] , thisTitle.c_str(), "p");
    }
    if( targetType == "NH3" ) mg->GetYaxis()->SetRangeUser(ymin,ymax);
    else if( targetType == "ND3" ) mg->GetYaxis()->SetRangeUser(ymin,ymax);
    mg->GetXaxis()->SetLimits(0.12,0.8);
    epoch += "_PFs";
    //TCanvas* mgplt = new TCanvas(epoch.c_str(),epoch.c_str(),800,600);
    //mgplt->cd();
    pad->cd();
    mg->Draw("ap"); leg->Draw("same");

    return;
    //return mgplt;

}


// This function is similar to the one above, "DF_Q2_Averaged_Plot", except it doesn't plot the
// average DFs. For each x, Q2 bin, I calculate the weighted average of the DF across all the variations
// in the "Epochs" vector. Then, using this average, I calculate the largest magnitude percent
// difference for each kinematic bin. I then find the average value of the magnitude percent differences.
// The "Latest" DataSet object is the value of DF that has the latest verion of the dilution factor
// systematic errors written to it to be saved later.
//double Systematic_Error_Calc( vector<DataSet>& Epochs, DataSet& Latest, string variation="" ){
//TCanvas* Systematic_Error_Calc( vector<DataSet>& Epochs, DataSet& Latest, string variation="", string epoch="" ){
/*
vector<TCanvas*> Systematic_Error_Calc( vector<DataSet>& Epochs, DataSet& Latest, string variation="", string epoch="" ){

    //string outName = "OutputDataTXT/Systematic_Errors_"+variation+".txt";
    //ofstream fout(outName.c_str());

    // These vectors are used to plot the percent differences vs x for each kinematic bin
    vector<double> maxPerDiffs, xVals, dfStatErr;

    // These vectors are used for plotting the means and standard deviations
    // "stdDevs" holds the standard deviation from 
    vector<double> xBins, stdDevs, noExpStatErr;

    // This tracks the average magnitude of the systematic error for each bin
    double avgSysSum = 0.0; double avgSysDen = 0.0; // Used for weighed average percent difference

    for(int i=0; i<Q2_Bin_Bounds.size()-1; i++){
	double qmid = (Q2_Bin_Bounds[i]+Q2_Bin_Bounds[i+1])/2.0;

	for(int j=0; j<X_Bin_Bounds.size()-1; j++){
	    double xmid = (X_Bin_Bounds[j]+X_Bin_Bounds[j+1])/2.0;

	    double weightedAvgDF = 0; double weightedAvgDFErr = 0;
	    double avgDF = 0; // Just the normal average without any weighting
	    double avgNum = 0.0; double avgDen = 0.0;
	    // I start with the "k=1" to skip the "NoExp" value, since that shouldn't be a benchmark
	    // for the systematics.
	    for(int k=1; k<Epochs.size(); k++){

		double df = Epochs[k].getDF_NH3( qmid, xmid );
		double dferr = Epochs[k].getErrDF_NH3( qmid, xmid );

		if( df > 0.0 && dferr > 0.0 ){
		    avgDF  += df;
		    avgNum += df / pow(dferr,2);
		    avgDen += 1.0 / pow(dferr,2);
		}
	
	    }
	    if( avgNum > 0.0 && avgDen > 0.0 ){
		avgDF = avgDF / ( 1.0 * (Epochs.size()-1) );
		weightedAvgDF = avgNum / avgDen;
		weightedAvgDFErr = 1.0 / sqrt( avgDen );
		// Because all epochs cover the same data sets, just with the DF calculated
		// differently, the weighted avg. X, Q2 should all be the same, so I just
		// grab them from the first epoch in the list
	        double xavg = Epochs[0].getAvgX( "All", qmid, xmid );
	        double q2avg= Epochs[0].getAvgQ2("All", qmid, xmid );
		// Now loop over the Epochs vector again to find the largest systematic difference
		double maxPercentDiff = 0.0;
		// Track the sum of the difference squared between the "avgDF" and each individual DF value
		double sumDiffSq = 0.0;
		//double largestDFErr = 0.0; // Largest statistical error on the set of DFs
		for(int k=1; k<Epochs.size(); k++){
		    double df = Epochs[k].getDF_NH3( qmid, xmid );
		    double dferr = Epochs[k].getErrDF_NH3( qmid, xmid );
		    double thisDiff = abs( (weightedAvgDF - df)/weightedAvgDF  );
	    	    if( maxPercentDiff < thisDiff && df > 0.0 ) maxPercentDiff = thisDiff;
		    if( df > 0.0 ) sumDiffSq += pow( ( df - avgDF ), 2 );
		    //if( largestDFErr < dferr ) largestDFErr = dferr;
		}
		// If there's a valid value of the percent difference and dferr...
		if( maxPercentDiff > 0.0001 ){
		    double percentErrOnDF = weightedAvgDFErr / weightedAvgDF; // Percent error on the average value of DF for this bin
		    avgSysSum += maxPercentDiff / pow( percentErrOnDF, 2 );
		    avgSysDen += 1.0 / pow( percentErrOnDF, 2);
		    maxPerDiffs.push_back( maxPercentDiff );
		    xVals.push_back( xavg ); // statistics-weighted value of x for this bin
		    dfStatErr.push_back( percentErrOnDF );

		    // Now calculate the standard deviation...
		    double stddev = sqrt( sumDiffSq /(1.0 * (Epochs.size()-2) ) ); // Minus 2 from the "n-1" term combined with needing only
										   // the size()-1 data from the vector.
		    stdDevs.push_back( stddev / avgDF );

		    // Now write the systematic error to the given bin in x, Q2 in the "Latest" configuration...
		    double latestDF = Latest.getDF_NH3( qmid, xmid );
		    Latest.SetSysErrDF( maxPercentDiff * latestDF, "NH3", qmid, xmid );
		}
		
	    }
	}

    }

    // Now calculate the average systematic error across all bins:
    
    double averageSysErr = 0.0;
    if( avgSysDen > 0.0 ) averageSysErr = avgSysSum / avgSysDen;

    cout << "Average systematic error for "<< epoch <<" in variation "<< variation <<": "<< averageSysErr << endl;

    // Draws a plot of the max percent differences
    TMultiGraph* mg    = new TMultiGraph();
    TMultiGraph* mgDev = new TMultiGraph();
    TLegend* leg    = new TLegend(0.1,0.7,0.4,0.9);
    TLegend* legDev = new TLegend(0.1,0.7,0.4,0.9);

    TGraph* PltPerDiff = new TGraph( xVals.size(), xVals.data(), maxPerDiffs.data() );
    //PltPerDiff->GetYaxis()->SetRangeUser(0,2*averageSysErr);
    //string title = variation+" Largest \% Difference; X; (D_{F,avg} - D_{F,largest})/D_{F,avg}";
    //PltPerDiff->SetTitle( title.c_str() );
    PltPerDiff->SetMarkerStyle( kFullCircle );
    PltPerDiff->SetMarkerColor( kRed );
    mg->Add( PltPerDiff, "p" );
    leg->AddEntry( PltPerDiff, "Systematic Error", "p" );

    TGraph* PltStdDev = new TGraph( xVals.size(), xVals.data(), stdDevs.data() );
    //PltStdDev->GetYaxis()->SetRangeUser(0,2*averageSysErr);
    //PltStdDev->SetTitle( title.c_str() );
    PltStdDev->SetMarkerStyle( kFullCircle );
    PltStdDev->SetMarkerColor( kRed );
    mgDev->Add( PltStdDev, "p" );
    legDev->AddEntry( PltStdDev, "Systematic Error", "p" );

    TGraph* PltdfStatErr = new TGraph( xVals.size(), xVals.data(), dfStatErr.data() );
    //PltdfStatErr->GetYaxis()->SetRangeUser(0,2*averageSysErr);
    //string title = variation+" Largest \% Difference; X; (D_{F,avg} - D_{F,largest})/D_{F,avg}";
    //PltdfStatErr->SetTitle( title.c_str() );
    PltdfStatErr->SetMarkerStyle( kFullSquare );
    PltdfStatErr->SetMarkerColor( kBlue );
    mg->Add( PltdfStatErr, "p" );
    mgDev->Add( PltdfStatErr, "p" );
    leg->AddEntry( PltdfStatErr, "Statistical Error", "p" );
    legDev->AddEntry( PltdfStatErr, "Statistical Error", "p" );

    string mgTitle = variation+" Statistical vs. Systematic Errors; X; Error";
    mg->SetTitle( mgTitle.c_str() );
    //mg->GetYaxis()->SetRangeUser(0,2*averageSysErr);
    mgDev->SetTitle( mgTitle.c_str() );
    mgDev->GetYaxis()->SetRangeUser(0,0.15);

    // Create a second legend for the reduced chi2
    TLegend* avgSysErrLeg = new TLegend(0.4,0.7,0.79,0.89);
    //TLegend* avgSysErrLeg = new TLegend(0.4,0.1,0.79,0.29);
    avgSysErrLeg->SetFillStyle(0);
    //avgSysErrLeg->SetFillColorAlpha( kWhite, 0.5);
    avgSysErrLeg->SetLineColorAlpha(0,0.5);
    avgSysErrLeg->SetLineWidth(0);

    //stringstream avgSysErrLegTitle; avgSysErrLegTitle << "Weighted Avg. Systematic Error";
    //avgSysErrLeg->AddEntry( avgSysErrLegTitle.str().c_str(), avgSysErrLegTitle.str().c_str(), "");
    //avgSysErrLeg->SetTextSize(0.04);
    stringstream avgSysErrLegNum; avgSysErrLegNum << "#delta Sys_{avg} =" << setprecision(3) << 100*averageSysErr <<"\%";
    avgSysErrLeg->AddEntry( avgSysErrLegNum.str().c_str(), avgSysErrLegNum.str().c_str(), "");
    avgSysErrLeg->SetTextSize(0.04);
//
    string outTitle = variation+"_plt";
    TCanvas* c = new TCanvas( outTitle.c_str(), outTitle.c_str(), 800, 600 );
    c->cd();
    mg->Draw("AP");
    leg->Draw("same");
    avgSysErrLeg->Draw("same");
//
    string outTitle = variation+"_stdDev_plt";
    TCanvas* c2 = new TCanvas( outTitle.c_str(), outTitle.c_str(), 800, 600 );
    c2->cd();
    mgDev->Draw("AP");
    legDev->Draw("same");

    //return averageSysErr;
    return c2;

}
*/
/*
// This version is convoluted lol
vector<TCanvas*> Systematic_Error_Calc( vector<DataSet>& Epochs, DataSet& Latest, string variation="", string epoch="" ){

    //string outName = "OutputDataTXT/Systematic_Errors_"+variation+".txt";
    //ofstream fout(outName.c_str());

    // These vectors are used to plot the percent differences vs x for each kinematic bin
    vector<double> maxPerDiffs, xVals;

    // These vectors are used for plotting the means and standard deviations
    // "stdDevs" holds the standard deviation from the average dilution factor
    vector<double> stdDevs;

    // Loop over the data
    for(int i=0; i<X_Bin_Bounds.size()-1; i++){
	double xmid = (X_Bin_Bounds[i]+X_Bin_Bounds[i+1])/2.0;
        
	double avgDF = 0; // Just the average value of DF for this Q2 bin, without including the "NoExp" value
	double wAvgDF = 0; // The weighted average of all DF for a Q2 bin, without the "NoExp" data
	double sumWeightSq = 0; // Sum of the weights squared; used in calculation of std. dev. of weighted average
	double meanXavg = 0; // The average of the weighted averages of x for each Q2 bin
        double counter = 0; // Used to track how many good values of DF there are
        bool goodDfBin = true; // Tracks if all the dilution factor values were reasonable (i.e. above 0.1)
        // I start with the "k=1" to skip the "NoExp" value, since that shouldn't be a benchmark
        // for the systematics.

	for(int j=0; j<Q2_Bin_Bounds.size()-1; j++){
	    double qmid = (Q2_Bin_Bounds[j]+Q2_Bin_Bounds[j+1])/2.0;
	    
	    for(int k=1; k<Epochs.size(); k++){

		double df = Epochs[k].getDF_NH3( qmid, xmid );
		double dferr = Epochs[k].getErrDF_NH3( qmid, xmid );
	        double xavg = Epochs[k].getAvgX( "All", qmid, xmid );

		//if( df < 0.1 || dferr <= 0.0 ) goodDfBin = false;
		
		if( df >= 0.1 && dferr > 0.0 ){
		    double weight = 1.0 / (dferr * dferr)
		    avgDF += df;
		    wAvgDF += df / weight;
		    counter++;
		    meanXavg += xavg;
		}
	    }
	}
        if( counter > 0 && avgDF > 0 && goodDfBin ){
	    avgDF = avgDF / ( 1.0 * counter );
	    meanXavg = meanXavg / ( 1.0 * counter );

	    //weightedAvgDFErr = 1.0 / sqrt( avgDen );
	    // Because all epochs cover the same data sets, just with the DF calculated
	    // differently, the weighted avg. X, Q2 should all be the same, so I just
	    // grab them from the first epoch ("NoExp") in the list
            //double xavg = Epochs[0].getAvgX( "All", qmid, xmid );
            //double q2avg= Epochs[0].getAvgQ2("All", qmid, xmid );

	    // Track the sum of the difference squared between the "avgDF" and each individual DF value
	    double sumDiffSq = 0.0;

	    // Tracks the largest percent difference between the given variation and the average value
	    double maxPercentDiff = 0.0;

	    //double largestDFErr = 0.0; // Largest statistical error on the set of DFs
	    //cout <<"For kinematic bin x = "<< xmid <<", Q2 = "<< qmid <<":\n";
	    //cout <<" Average DF = "<< avgDF << endl;
	    double backupCounter = 0.0; // Just to make sure nothing went wrong
	    for(int j=0; j<Q2_Bin_Bounds.size()-1; j++){
	        double qmid = (Q2_Bin_Bounds[j]+Q2_Bin_Bounds[j+1])/2.0;

		for(int k=1; k<Epochs.size(); k++){
		    double df = Epochs[k].getDF_NH3( qmid, xmid );
		    double dferr = Epochs[k].getErrDF_NH3( qmid, xmid );
		    double thisDiff = abs( avgDF - df ) / avgDF;
		    if( df > 0.1 && dferr > 0.0 ){
		        sumDiffSq += pow( ( df - avgDF ), 2 );
		        if( thisDiff > maxPercentDiff ) maxPercentDiff = thisDiff;
		        backupCounter++;
		    }
		}
	        //double percentErrOnDF =  Epochs[0].getErrDF_NH3( qmid, xmid ) / Epochs[0].getDF_NH3( qmid, xmid ); // Percent error on the average value 
													           // of the "NoExp" dilution factor
	    }
	    //double percentErrOnDF =  Epochs[0].getErrDF_NH3( qmid, xmid ) / Epochs[0].getDF_NH3( qmid, xmid ); // Percent error on the average value 
													       // of the "NoExp" dilution factor
	    if( counter != backupCounter ) cout << "ERROR: Misalignment between expected number of entries. DF value invalid.\n";

	    maxPerDiffs.push_back( maxPercentDiff );

	    xVals.push_back( meanXavg ); // average value of the stat.-weighted mean values of x for each Q2 bin
	    //noExpStatErr.push_back( percentErrOnDF );

	    // Now calculate the standard deviation...
	    double stddev = sqrt( sumDiffSq /( (1.0 * counter) - 1.0 ) );

	    //if( stddev / avgDF > 0.1 ){
		//cout << "Large error at x = "<< xmid <<", Q2 = "<< qmid <<": dDF_sys = "<< 100 * stddev / avgDF <<" \%\n";
	    //}

	    stdDevs.push_back( stddev / avgDF );
	    // Now write the systematic error to the given bin in x, Q2 in the "Latest" configuration...
	    //double latestDF = Latest.getDF_NH3( qmid, xmid );
	    //Latest.SetSysErrDF( maxPercentDiff * latestDF, "NH3", qmid, xmid );
	}
    }
    
    // Get all the statistical errors for the "NoExp" data for plotting
    vector<double> noExpStatErr, noExpXvals;
    for(int i=0; i<Q2_Bin_Bounds.size()-1; i++){
	double qmid = (Q2_Bin_Bounds[i] + Q2_Bin_Bounds[i+1]) / 2.0;

	for(int j=0; j<X_Bin_Bounds.size()-1; j++){
	    double xmid = (X_Bin_Bounds[j] + X_Bin_Bounds[j+1]) / 2.0;
	    double df = Epochs[0].getDF_NH3( qmid, xmid );
	    double dferr = Epochs[0].getErrDF_NH3( qmid, xmid );
	    if( df > 0.1 && dferr > 0.0 ){
		noExpStatErr.push_back( dferr / df );
		noExpXvals.push_back( Epochs[0].getAvgX( "All", qmid, xmid ) );
	    }
	}
    }

    // Now calculate the average systematic error across all bins:
    
    // Draws a plot of the max percent differences
    TMultiGraph* mg    = new TMultiGraph();
    TMultiGraph* mgDev = new TMultiGraph();
    TLegend* leg    = new TLegend(0.1,0.7,0.4,0.9);
    TLegend* legDev = new TLegend(0.1,0.7,0.4,0.9);

    // Plot the percent differences (as the systematic errors)
    TGraph* PltPerDiff = new TGraph( xVals.size(), xVals.data(), maxPerDiffs.data() );
    PltPerDiff->SetMarkerStyle( kFullCircle );
    PltPerDiff->SetMarkerColor( kRed );
    mg->Add( PltPerDiff, "p" );
    leg->AddEntry( PltPerDiff, "Systematic Error", "p" );

    // Plot the standard deviations (as the systematic errors)
    TGraph* PltStdDev = new TGraph( xVals.size(), xVals.data(), stdDevs.data() );
    PltStdDev->SetMarkerStyle( kFullCircle );
    PltStdDev->SetMarkerColor( kRed );
    mgDev->Add( PltStdDev, "p" );
    legDev->AddEntry( PltStdDev, "Systematic Error", "p" );

    // Plot the statistical errors from the "NoExp" plots
    TGraph* PltdfStatErr = new TGraph( noExpXvals.size(), noExpXvals.data(), noExpStatErr.data() );
    PltdfStatErr->SetMarkerStyle( kFullSquare );
    PltdfStatErr->SetMarkerColor( kBlue );
    mg->Add( PltdfStatErr, "p" );
    mgDev->Add( PltdfStatErr, "p" );
    leg->AddEntry( PltdfStatErr, "Statistical Error", "p" );
    legDev->AddEntry( PltdfStatErr, "Statistical Error", "p" );

    // Set the titles for the multigraph plots
    string mgTitle = variation+" Errors (\%-Difference Method); X; Error";
    mg->SetTitle( mgTitle.c_str() );
    mg->GetYaxis()->SetRangeUser(0,0.2);
    //mg->GetYaxis()->SetRangeUser(0,2*averageSysErr);
    string mgTitleDev = variation+" Errors (Standard Deviation Method); X; Error";
    mgDev->SetTitle( mgTitleDev.c_str() );
    mgDev->GetYaxis()->SetRangeUser(0,0.2);

    // Create the canvas with the percent difference systematic errors
    string outTitle = variation+"_plt";
    TCanvas* c = new TCanvas( outTitle.c_str(), outTitle.c_str(), 800, 600 );
    c->cd();
    mg->Draw("AP");
    leg->Draw("same");
    //avgSysErrLeg->Draw("same");

    // Create the canvas with the 
    outTitle = variation+"_stdDev_plt";
    TCanvas* c2 = new TCanvas( outTitle.c_str(), outTitle.c_str(), 800, 600 );
    c2->cd();
    mgDev->Draw("AP");
    legDev->Draw("same");

    vector<TCanvas*> Plots = {c, c2};
    //return averageSysErr;
    return Plots;

}*/

// This version gets rid of all the weighted average crap.
//vector<TCanvas*> Systematic_Error_Calc( vector<DataSet>& Epochs, string variation="", string epoch="", bool useFirstEntry=false ){
vector<TCanvas*> Systematic_Error_Calc( vector<DataSet>& Epochs, string variation="", int EpochNum=0, bool useFirstEntry=false ){

    //string outName = "OutputDataTXT/Systematic_Errors_"+variation+".txt";
    //ofstream fout(outName.c_str());

    // These vectors are used to plot the percent differences vs x for each kinematic bin
    vector<double> maxPerDiffs, xVals, noExpStatErr, q2Vals;
    vector<double> PFxVals, PFnoExpStatErr, PFq2Vals; // values for plotting the packing fraction values

    // These vectors are used for plotting the means and standard deviations
    // "stdDevs" holds the standard deviation from the average dilution factor
    vector<double> xBins, stdDevs;
    vector<double> PFxBins, PFstdDevs; // values for the packing fraction values

    // Loop over the data
    for(int i=0; i<Q2_Bin_Bounds.size()-1; i++){
	double qmid = (Q2_Bin_Bounds[i]+Q2_Bin_Bounds[i+1])/2.0;

	for(int j=0; j<X_Bin_Bounds.size()-1; j++){
	    double xmid = (X_Bin_Bounds[j]+X_Bin_Bounds[j+1])/2.0;

	    double avgDF = 0; // Just the average value of DF for this x, Q2 bin, without including the "NoExp" value
	    double avgPF = 0; // Same as above, but using the PF values
	    double counter = 0; // Used to track how many good values of DF there are
	    double counterPF=0; // Same as above, but using PF values
	    bool goodDfBin = true; // Tracks if all the dilution factor values were reasonable (i.e. above 0.1)
	    bool goodPfBin = true; // Same, but tracking the PF values
	    // I start with the "k=1" to skip the "NoExp" value, since that shouldn't be a benchmark
	    // for the systematics.
	    // HOWEVER: There is an exception for the carbon data, where "NoExp" is included as a parameter
	    int StartingIt = 1; // Iterator to loop over the Epochs vector
	    if( useFirstEntry ){ /*cout <<"Using NoExp data for this calculation!\n";*/ StartingIt = 0; }
	    
	    for(int k=StartingIt; k<Epochs.size(); k++){

		double df = Epochs[k].getDF_NH3( qmid, xmid );
		double dferr = Epochs[k].getErrDF_NH3( qmid, xmid );
		if( df < 0.1 ) goodDfBin = false;
		if( df > 0.0 && dferr > 0.0 ){
		    avgDF  += df;
		    counter++;
		}
	
		double pf = Epochs[k].getPF_bath_NH3( qmid, xmid );
		double pferr = Epochs[k].getErrPF_bath_NH3( qmid, xmid );
		if( pf < 0.1 ) goodPfBin = false;
		if( pf > 0.0 && pferr > 0.0 ){
		    avgPF  += pf;
		    counterPF++;
		}

	    }
	    // Calculate the differences for DF if all values are valid
	    if( counter > 0 && avgDF > 0 && goodDfBin ){
		avgDF = avgDF / ( 1.0 * counter );
		// Because all epochs cover the same data sets, just with the DF calculated
		// differently, the weighted avg. X, Q2 should all be the same, so I just
		// grab them from the first epoch ("NoExp") in the list
	        double xavg = Epochs[0].getAvgX( "All", qmid, xmid );
	        double q2avg= Epochs[0].getAvgQ2("All", qmid, xmid );

		// Now loop over the Epochs vector again to find the largest systematic difference
		double maxPercentDiff = 0.0;

		// Track the sum of the difference squared between the "avgDF" and each individual DF value
		double sumDiffSq = 0.0;

		for(int k=StartingIt; k<Epochs.size(); k++){
		    double df = Epochs[k].getDF_NH3( qmid, xmid );
		    double thisDiff = abs( avgDF - df ) / avgDF;

	    	    if( maxPercentDiff < thisDiff && df > 0.0 ) maxPercentDiff = thisDiff;
		    if( df > 0.0 ) sumDiffSq += pow( ( df - avgDF ), 2 );
		    else maxPercentDiff = 100000; // Used as an error condition; if there's an invalid DF, the calculation for this bin is skipped (very sloppy)
		}
		// If there's a valid value of the percent difference and dferr...
		if( maxPercentDiff > 0.0001 && maxPercentDiff != 100000 ){
		    double percentErrOnDF =  Epochs[0].getErrDF_NH3( qmid, xmid ) / Epochs[0].getDF_NH3( qmid, xmid ); // Percent error on the average value 
														       // of the "NoExp" dilution factor
		    maxPerDiffs.push_back( maxPercentDiff );
		    xVals.push_back( xavg ); // statistics-weighted value of x for this bin from the "NoExp" data set
		    q2Vals.push_back( q2avg );// stat-weighted value of q2 for this bin from the "NoExp" data set
		    noExpStatErr.push_back( percentErrOnDF );

		    // Now calculate the standard deviation...
		    double stddev = sqrt( sumDiffSq /( (1.0 * counter) - 1.0 ) ); // Minus 2 from the "n-1" term combined with needing only
										   // the size()-1 data from the vector.
		    if( stddev / avgDF > 0.1 ){
			cout << "Large error at x = "<< xmid <<", Q2 = "<< qmid <<": dDF_sys = "<< 100 * stddev / avgDF <<" \%\n";
		    }
		    stdDevs.push_back( stddev / avgDF );

		    // Now write the systematic error to the given bin in x, Q2 in the "NoExp" configuration...
		    //double noexpDF = Epochs[0].getDF_NH3( qmid, xmid );
		    //Epochs[0].SetSysErrDF( (stddev / avgDF) * noexpDF, "NH3", qmid, xmid );
		}	
	    } // End of "if" statement for DFs

	    // Calculate the differences for PF if all values are valid
	    if( counterPF > 0 && avgPF > 0 && goodPfBin ){
		avgPF = avgPF / ( 1.0 * counterPF );
		// Because all epochs cover the same data sets, just with the DF calculated
		// differently, the weighted avg. X, Q2 should all be the same, so I just
		// grab them from the first epoch ("NoExp") in the list
	        double xavg = Epochs[0].getAvgX( "All", qmid, xmid );
	        double q2avg= Epochs[0].getAvgQ2("All", qmid, xmid );

		// Now loop over the Epochs vector again to find the largest systematic difference
		double maxPercentDiff = 0.0;

		// Track the sum of the difference squared between the "avgDF" and each individual DF value
		double sumDiffSq = 0.0;

		for(int k=StartingIt; k<Epochs.size(); k++){
		    double pf = Epochs[k].getPF_bath_NH3( qmid, xmid );
		    double thisDiff = abs( avgPF - pf ) / avgPF;

	    	    if( maxPercentDiff < thisDiff && pf > 0.0 ) maxPercentDiff = thisDiff;
		    if( pf > 0.0 ) sumDiffSq += pow( ( pf - avgPF ), 2 );
		    else maxPercentDiff = 100000; // Used as an error condition; if there's an invalid DF, the calculation for this bin is skipped (very sloppy)
		}
		// If there's a valid value of the percent difference and dferr...
		if( maxPercentDiff > 0.0001 && maxPercentDiff != 100000 ){
		    double percentErrOnPF =  Epochs[0].getErrPF_bath_NH3( qmid, xmid ) / Epochs[0].getPF_bath_NH3( qmid, xmid ); // Percent error on the average value 
														                 // of the "NoExp" packing fraction
		    //maxPerDiffs.push_back( maxPercentDiff );
		    PFxVals.push_back( xavg ); // statistics-weighted value of x for this bin from the "NoExp" data set
		    PFq2Vals.push_back( q2avg ); // same thing as above but for Q2 values
		    PFnoExpStatErr.push_back( percentErrOnPF );

		    // Now calculate the standard deviation...
		    double stddev = sqrt( sumDiffSq /( (1.0 * counterPF) - 1.0 ) ); // Minus 2 from the "n-1" term combined with needing only
										   // the size()-1 data from the vector.
		    if( stddev / avgPF > 0.1 ){
			cout << "Large error at x = "<< xmid <<", Q2 = "<< qmid <<": dPF_sys = "<< 100 * stddev / avgPF <<" \%\n";
		    }
		    PFstdDevs.push_back( stddev / avgPF );

		    // Now write the systematic error to the given bin in x, Q2 in the "NoExp" configuration...
		    //double noexpPF = Epochs[0].getPF_bath_NH3( qmid, xmid );
		    //Epochs[0].SetSysErrDF( (stddev / avgDF) * noexpDF, "NH3", qmid, xmid );
		}	
	    } // End of "if" statement for DFs

	}

    }
    
    
    // Draws a plot of the systematics for the DFs
    TMultiGraph* mg     = new TMultiGraph();
    TMultiGraph* mgDev  = new TMultiGraph();
    TMultiGraph* mgFrac = new TMultiGraph();
    TLegend* leg     = new TLegend(0.1,0.7,0.4,0.9);
    TLegend* legDev  = new TLegend(0.1,0.7,0.4,0.9);
    TLegend* legFrac = new TLegend(0.1,0.8,0.5,0.9);

    // Draws plots for the PF systematics
    TMultiGraph* mgPfDev = new TMultiGraph();
    TLegend* legPfDev = new TLegend(0.1,0.7,0.4,0.9);

    // Plot the percent differences (as the systematic errors)
    TGraph* PltPerDiff = new TGraph( xVals.size(), xVals.data(), maxPerDiffs.data() );
    PltPerDiff->SetMarkerStyle( kFullCircle );
    PltPerDiff->SetMarkerColor( kRed );
    mg->Add( PltPerDiff, "p" );
    leg->AddEntry( PltPerDiff, "Systematic Error", "p" );

    // Plot the standard deviations for the DFs (as the systematic errors)
    TGraph* PltStdDev = new TGraph( xVals.size(), xVals.data(), stdDevs.data() );

    // Fit the standard deviation plot using a third degree polynomial
    string cuFitName = variation+"_fit";
    TF1* cubicFit = new TF1( cuFitName.c_str(), "pol2", 0, 1 );
    cubicFit->SetLineColor( kBlack );
    cubicFit->SetLineWidth(1);
    TFitResultPtr r;
    r = PltStdDev->Fit( cubicFit, "S", "Q");
    // Plot the fit parameters
    TLegend* legFit = new TLegend(0.4,0.7,0.9,0.9);
    string legFitName = cuFitName+"_leg";
    legFit->SetNColumns(2); 
    //legFit->SetHeader("Systematic Error Fit #delta(x) = A + Bx + Cx^{2} + Dx^{3}","C");
    legFit->SetHeader("Systematic Error Fit #delta(x) = A + Bx + Cx^{2}","C");
    legFit->AddEntry( cubicFit, "Sys. Err. Fit #delta(x)", "l");
    //legFit->SetFillColorAlpha(0, 0); // Make background transparent
    stringstream p1; p1 << "A = "<< setprecision(2) << r->Value(0) <<" #pm "<< r->Error(0);
    legFit->AddEntry( p1.str().c_str(), p1.str().c_str(), "");
    stringstream p2; p2 << "B = "<< setprecision(2) << r->Value(1) <<" #pm "<< r->Error(1);
    legFit->AddEntry( p2.str().c_str(), p2.str().c_str(), "");
    stringstream p3; p3 << "C = "<< setprecision(2) << r->Value(2) <<" #pm "<< r->Error(2);
    legFit->AddEntry( p3.str().c_str(), p3.str().c_str(), "");
    //stringstream p4; p4 << "D = "<< setprecision(2) << r->Value(2) <<" #pm "<< r->Error(2);
    //legFit->AddEntry( p4.str().c_str(), p4.str().c_str(), "");

    PltStdDev->SetMarkerStyle( kFullCircle );
    PltStdDev->SetMarkerColor( kRed );
    mgDev->Add( PltStdDev, "p" );
    legDev->AddEntry( PltStdDev, "Systematic Error", "p" );

    // Plot the statistical errors from the "NoExp" plots
    TGraph* PltdfStatErr = new TGraph( xVals.size(), xVals.data(), noExpStatErr.data() );
    PltdfStatErr->SetMarkerStyle( kFullSquare );
    PltdfStatErr->SetMarkerColor( kBlue );
    mg->Add( PltdfStatErr, "p" );
    mgDev->Add( PltdfStatErr, "p" );
    leg->AddEntry( PltdfStatErr, "Statistical Error", "p" );
    legDev->AddEntry( PltdfStatErr, "Statistical Error", "p" );

    // Set the titles for the multigraph plots
    string mgTitle = variation+" D_{F} Errors (\%-Difference Method); X; Error D_{F}";
    mg->SetTitle( mgTitle.c_str() );
    mg->GetYaxis()->SetRangeUser(0,0.2);
    mg->GetXaxis()->SetLimits(0.12,0.8);
    //mg->GetYaxis()->SetRangeUser(0,2*averageSysErr);
    string mgTitleDev = variation+" D_{F} Errors; X; Error D_{F}"; // Used for the standard deviation method
    mgDev->SetTitle( mgTitleDev.c_str() );
    mgDev->GetYaxis()->SetRangeUser(0,0.2);
    mgDev->GetXaxis()->SetLimits(0.12,0.8);

    // Now make plots for the PFs
    // Plot the standard deviations for the PFs (as the systematic errors)
    TGraph* PltPFStdDev = new TGraph( PFxVals.size(), PFxVals.data(), PFstdDevs.data() );

    // Fit the standard deviation plot using a third degree polynomial
    cuFitName = variation+"_PFfit";
    TF1* cubicFitPF = new TF1( cuFitName.c_str(), "[0]", 0, 1 );
    cubicFitPF->SetLineColor( kBlack );
    cubicFitPF->SetLineWidth(1);
    //TFitResultPtr r;
    r = PltPFStdDev->Fit( cubicFitPF, "S", "Q");
    // Plot the fit parameters
    TLegend* legFitPF = new TLegend(0.4,0.7,0.9,0.9);
    string legFitPFName = cuFitName+"_pfleg";
    //legFitPF->SetNColumns(2); 
    //legFitPF->SetHeader("Systematic Error Fit #delta(x) = A + Bx","C");
    legFitPF->SetHeader("Systematic Error Fit = A","C");

    //legFitPF->SetFillColorAlpha(0, 0); // Make background transparent
    legFitPF->AddEntry( cubicFitPF, "Sys. Err. Fit #delta(x)","l");
    stringstream P1; P1 << "A = "<< setprecision(2) << r->Value(0) <<" #pm "<< r->Error(0);
    legFitPF->AddEntry( P1.str().c_str(), P1.str().c_str(), "");
    //stringstream P2; P2 << "B = "<< setprecision(2) << r->Value(1) <<" #pm "<< r->Error(1);
    //legFitPF->AddEntry( P2.str().c_str(), P2.str().c_str(), "");

    PltPFStdDev->SetMarkerStyle( kFullCircle );
    PltPFStdDev->SetMarkerColor( kViolet );
    mgPfDev->Add( PltPFStdDev, "p" );
    legPfDev->AddEntry( PltPFStdDev, "Systematic Error", "p" );

    // Plot the PF statistical errors from the "NoExp" plots
    TGraph* pltPFStatErr = new TGraph( PFxVals.size(), PFxVals.data(), PFnoExpStatErr.data() );
    pltPFStatErr->SetMarkerStyle( kFullSquare );
    pltPFStatErr->SetMarkerColor( kGreen );
    mgPfDev->Add( pltPFStatErr, "p" );
    legPfDev->AddEntry( pltPFStatErr, "Statistical Error", "p" );

    string mgTitlePFDev = variation+" P_{F} Errors; X; Error P_{F}"; // Used for the standard deviation method
    mgPfDev->SetTitle( mgTitlePFDev.c_str() );
    mgPfDev->GetYaxis()->SetRangeUser(0,0.2);
    mgPfDev->GetXaxis()->SetLimits(0.12,0.8);

    // Create the canvas with the percent difference systematic errors
    string outTitle = variation+"_plt";
    TCanvas* c = new TCanvas( outTitle.c_str(), outTitle.c_str(), 800, 600 );
    c->cd();
    mg->Draw("AP");
    leg->Draw("same");
    //avgSysErrLeg->Draw("same");

    // Create the canvas with the standard deviations
    outTitle = variation+"_stdDev_plt";
    TCanvas* c2 = new TCanvas( outTitle.c_str(), outTitle.c_str(), 800, 600 );
    c2->cd();
    mgDev->Draw("AP");
    cubicFit->Draw("same");
    legDev->Draw("same");
    legFit->Draw("same");

    // Create the canvas for the PF systematics
    outTitle = variation+"_PFstdDev_plt";
    TCanvas* c4 = new TCanvas( outTitle.c_str(), outTitle.c_str(), 800, 600 );
    c4->cd();
    mgPfDev->GetYaxis()->SetRangeUser(0,0.1);
    mgPfDev->Draw("AP");
    cubicFitPF->Draw("same");
    legPfDev->Draw("same");
    legFitPF->Draw("same");

    // Now create a plot where the systematic and statistical errors are multiplied by the values of DF
    // for that given bin in x, Q2. All values of DF used here are the errors from the "NoExp" data.
    // First, generate the scaled values using the standard error systematics and statistical errors
    vector<double> FracSysErr, FracStatErr, FracX;
    for(int i=0; i<X_Bin_Bounds.size()-1; i++){
	double xmid = (X_Bin_Bounds[i] + X_Bin_Bounds[i+1]) / 2.0;

	for(int j=0; j<Q2_Bin_Bounds.size()-1; j++){
	    double qmid = (Q2_Bin_Bounds[j] + Q2_Bin_Bounds[j+1]) / 2.0;
	    // Get the systematic and statistical errors for this bin
	    double df     = Epochs[0].getDF_NH3( qmid, xmid );
	    double dfStat = Epochs[0].getErrDF_NH3( qmid, xmid );
	    double dfSys  = Epochs[0].getSysErrDF_NH3( qmid, xmid );
	    if( df > 0 && dfStat > 0 && dfSys > 0 ){
		FracStatErr.push_back( dfStat );
		FracSysErr.push_back( dfSys );
		FracX.push_back( Epochs[0].getAvgX("All", qmid, xmid) );
	    }
	}
    }

    // Plot the fractional statistical errors
    TGraph* PltFracStatErrs = new TGraph( FracX.size(), FracX.data(), FracStatErr.data() );
    PltFracStatErrs->SetMarkerStyle(kFullSquare);
    PltFracStatErrs->SetMarkerColor(kBlue);
    string fstatName  = variation+"_pltfracstat"; PltFracStatErrs->SetName( fstatName.c_str() );
    mgFrac->Add( PltFracStatErrs, "p" );
    legFrac->AddEntry( PltFracStatErrs, "Total Error D_{F}*#delta_{stat}", "p" );

    TGraph* PltFracSysErrs = new TGraph( FracX.size(), FracX.data(), FracSysErr.data() );
    PltFracSysErrs->SetMarkerStyle(kFullCircle);
    PltFracSysErrs->SetMarkerColor(kRed);
    fstatName  = variation+"_pltfracstat"; PltFracSysErrs->SetName( fstatName.c_str() );
    mgFrac->Add( PltFracSysErrs, "p" );
    legFrac->AddEntry( PltFracSysErrs, "Total Error D_{F}*#delta_{sys}", "p" );

    mgTitle = variation+" Dilution Factor Error Sizes; X; #delta D_{F}";
    mgFrac->SetTitle( mgTitle.c_str() );

    // Create the canvas for this plot
    outTitle = variation+"_frac_plt";
    TCanvas* c3 = new TCanvas( outTitle.c_str(), outTitle.c_str(), 800, 600 );
    c3->cd();
    mgFrac->GetYaxis()->SetRangeUser(0,0.05);
    mgFrac->GetYaxis()->SetMaxDigits(1);
    mgFrac->GetXaxis()->SetRangeUser(0.1,0.8);
    mgFrac->Draw("ap");
    legFrac->Draw("same");

    // Hol' up: before you go, write the standard deviations and statistical errors to an output file 
    // to be used later for finalizing the systematics
    replace( variation.begin(), variation.end(), ' ', '_' ); // replace the space with an underscore
    string outFileDF = "Final_Plots/Text_Files/"+variation+"_DF_StdDev.txt";
    ofstream foutDF( outFileDF.c_str() );
    foutDF <<"# x_Avg   q2_Avg   DF_StdDev   DF_StatErr\n";
    for( size_t i=0; i<xVals.size(); i++ ){
	foutDF << xVals[i] <<"  "<< q2Vals[i] <<"  "<< stdDevs[i] <<"  "<< noExpStatErr[i] << endl;
    }
    foutDF.close();

    string outFilePF = "Final_Plots/Text_Files/"+variation+"_PF_StdDev.txt";
    ofstream foutPF( outFilePF.c_str() );
    foutPF <<"# x_Avg   q2_Avg   PF_StdDev   PF_StatErr\n";
    for( size_t i=0; i<xVals.size(); i++ ){
	foutPF << PFxVals[i] <<"  "<< PFq2Vals[i] <<"  "<< PFstdDevs[i] <<"  "<< PFnoExpStatErr[i] << endl;
    }
    foutPF.close();

    vector<TCanvas*> Plots = {c, c2, c3, c4};
    //return averageSysErr;
    return Plots;

}


/*
TCanvas* AllRaw_Q2_Averaged_Plot( vector<DataSet>& Epochs, string epoch, string targetType, string thruPF ){

    vector<int> palette = { 1, 632, 800, 400, 416, 600, 840, 880, 900, 616, 820, 432, 920, 860 };
    int pint = 0;

    if( !(targetType == "NH3" || targetType == "ND3") ){
	    cout <<"ERROR: enter a valid target type to plot: 'NH3' or 'ND3'\n";
	    TCanvas* emptyC = new TCanvas();
	    return emptyC;
    }

    // Check to see if we're using the average PF across all bins (thruPF = "YES") or if the
    // PF is being calculated for each bin separately (thruPF = else)
    //vector<double> xbins, zeros, avg_dfs, avg_df_errs;
    vector<double> modXbin, modZeros, modAvgdfs; // Model values for comparison
    string toFile = THIS_DIR + "Darren_Model_DF.txt";
    // Loop over the model data first
    DataSet ModelDFs;
    ifstream fin(toFile.c_str()); string line;
    while(getline(fin,line)){
	stringstream sin(line);
	double q2, x, qbin, xbin, df, dferr;
	sin >> q2 >> x >> qbin >> xbin >> df >> dferr;
        ModelDFs.SetDFs( q2+0.1, x+0.01, df, 0 );
    }
    fin.close();
    // Now read in the data
 TLegend* leg = new TLegend(0.1,0.7,0.3,0.9);
 leg->SetHeader("Epochs","C"); leg->SetNColumns(2);
 TMultiGraph* mg = new TMultiGraph();
 mg->SetTitle(epoch.c_str());
 vector<TGraphErrors*> DF_Epoch_Plots;
 for(size_t k=0; k<Epochs.size(); k++){
    Epochs[k].CalculateAllRaw();
    if( thruPF == "YES" ){
	Epochs[k].CalculateDFThruPF( targetType );
    }
    else{
	Epochs[k].CalculateDF();
    }

    DataSet AllData = Epochs[k]; // Very inefficient!!! Needs improvement!!!!11!
    vector<double> xbins, zeros, avg_dfs, avg_df_errs;

    for(size_t i=0; i<X_Bin_Bounds.size()-1; i++){
	double xmid = (X_Bin_Bounds[i]+X_Bin_Bounds[i+1])/2.0;
	vector<double> q2bins, dfvals, dferrs, modelq2bins, modeldfvals, modeldferrs;
	for(size_t j=0; j<Q2_Bin_Bounds.size()-1; j++){
	    double qmid = (Q2_Bin_Bounds[j]+Q2_Bin_Bounds[j+1])/2.0;
	    double thisDF = 0; double thisErrDF = 0;
	    double modelDF = ModelDFs.getDF_FixedPF_NH3( qmid, xmid );
	    double modelDFerr = 0;
	    //if( targetType == "NH3" && thruPF != "YES" ){
	        thisDF = AllData.getAllRaw(qmid, xmid, targetType);//AllData.getDF_NH3(qmid, xmid);
	        thisErrDF = AllData.getAllRawErr(qmid, xmid, targetType);//getErrDF_NH3(qmid, xmid);
	    //}
            if( thisDF != 0.0 ){
		    q2bins.push_back(qmid);
		    dfvals.push_back( thisDF );
		    dferrs.push_back( thisErrDF );
	    }
	    if( modelDF != 0.0 && k==0 ){
		    modelq2bins.push_back(qmid);
		    modeldfvals.push_back( modelDF );
		    modeldferrs.push_back( modelDFerr );	
	    }
	}
	if( q2bins.size() > 0 ){
	    // Calculate the weighted average of all the DFs in the given xbin
	    double num = 0.0; double den = 0.0;
	    for(size_t k=0; k<q2bins.size(); k++){
		num += dfvals[k] / pow(dferrs[k],2);
		den += 1.0 / pow(dferrs[k],2);
	    }
	    if( den != 0.0 ){
		xbins.push_back(xmid); zeros.push_back(0.0);
		avg_dfs.push_back( num / den ); avg_df_errs.push_back( 1.0 / sqrt(den) );
	    }
	}
	if( modelq2bins.size() > 0 && k==0 ){
	    // Calculate the weighted average of all the DFs in the given xbin
	    double num = 0.0; double den = 0.0;
	    for(size_t k=0; k<modelq2bins.size(); k++){
		num += modeldfvals[k] ;/// pow(modeldferrs[k],2);
		den += 1.0 ;/// pow(modeldferrs[k],2);
	    }
	    if( den != 0.0 ){
		modXbin.push_back(xmid); modAvgdfs.push_back( num/den ); modZeros.push_back(0.0);
	    }
	}

    }
    //cout << "Epoch "<< epoch <<" first bin value: "<< avg_dfs[0] << endl;
    TGraphErrors* thisGR = new TGraphErrors(xbins.size(), xbins.data(), avg_dfs.data(), zeros.data(), avg_df_errs.data() );
    DF_Epoch_Plots.push_back(thisGR);
 }
    DF_Epoch_Plots[0]->SetMarkerColor(kRed);   DF_Epoch_Plots[0]->SetMarkerStyle(kFullCircle);
    DF_Epoch_Plots[1]->SetMarkerColor(kOrange);DF_Epoch_Plots[1]->SetMarkerStyle(kFullCircle);
    DF_Epoch_Plots[2]->SetMarkerColor(kGreen); DF_Epoch_Plots[2]->SetMarkerStyle(kFullTriangleUp);
    DF_Epoch_Plots[3]->SetMarkerColor(kBlue);  DF_Epoch_Plots[3]->SetMarkerStyle(kFullTriangleUp);
    DF_Epoch_Plots[4]->SetMarkerColor(kViolet);DF_Epoch_Plots[4]->SetMarkerStyle(kFullTriangleDown);
    DF_Epoch_Plots[5]->SetMarkerColor(kTeal);  DF_Epoch_Plots[5]->SetMarkerStyle(kFullTriangleDown);
    DF_Epoch_Plots[6]->SetMarkerColor(kBlack);  DF_Epoch_Plots[6]->SetMarkerStyle(kFullSquare);


    for(size_t i=0; i<DF_Epoch_Plots.size()-1; i++){
        string thisTitle = "Epoch "+to_string(i+2); //DF_Epoch_Plots[i]->SetTitle(thisTitle.c_str());
        //DF_Epoch_Plots[i]->SetMarkerColor( palette[i+1] );
        //DF_Epoch_Plots[i]->SetMarkerStyle(i);
        mg->Add( DF_Epoch_Plots[i] , "p");
        leg->AddEntry( DF_Epoch_Plots[i] , thisTitle.c_str(), "p");
    }
    leg->AddEntry( DF_Epoch_Plots[6], "All Su22 Data","p");
    mg->Add( DF_Epoch_Plots[6], "p");
    mg->GetYaxis()->SetRangeUser(0,0.4);
    mg->GetXaxis()->SetLimits(0,1);
    TCanvas* mgplt = new TCanvas("mgplt2","mgplt2",800,600);
    mgplt->cd();
    mg->Draw("ap"); leg->Draw("same");

    return mgplt;

}

// Using the outputs from the proceeding "DF_Q2_Averaged_Plot" function, this function takes
// a vector of pointers to the TGraphErrors objects produced by the previous function and
// plots them together.
TCanvas* Averaged_DF_By_Sector( vector<TGraphErrors*>& Graphs, string epoch ){

    vector<int> palette = { 1, 632, 800, 400, 416, 600, 840, 880, 900, 616, 820, 432, 920, 860 };
    int pint = 0;

    TMultiGraph* mg = new TMultiGraph();
    TLegend* mgLeg = new TLegend(0.1,0.7,0.4,0.9);
    mgLeg->SetNColumns(3);
    mgLeg->SetHeader("Run Range Epochs","C");

    for(size_t i=0; i<Graphs.size(); i++){
        Graphs[i]->SetMarkerColor( palette[pint] );
	//Graphs[i]->GetXaxis()->SetRangeUser(0,0.8);
	if(pint < 4 ){ Graphs[i]->SetMarkerStyle(kFullCircle); }
	else if(pint >= 4 && pint < 8){ Graphs[i]->SetMarkerStyle(kFullTriangleUp); }
	else if(pint >= 8){ Graphs[i]->SetMarkerStyle(kFullSquare); }
	Graphs[i]->GetYaxis()->SetRangeUser(0,0.4);
	mg->Add( Graphs[i], "p" );
	string legTitle = Graphs[i]->GetTitle();
        mgLeg->AddEntry( Graphs[i], legTitle.c_str(), "p");
	pint++;
    }

    string title = "Average D_{F} for "+ epoch +"; X; D_{F}";
    //mg->SetTitle("Average D_{F} for 2022 Epochs; X; D_{F}");
    mg->SetTitle(title.c_str());
    
    TCanvas* AvgDF = new TCanvas("AvgDF","AvgDF",800,600);
    AvgDF->cd();
    mg->Draw("ap");
    mgLeg->Draw("same");

    return AvgDF;
}

// This function creates plots of the ratio of other target counts normalized to the NH3 target counts.
// Each target type is first normalized to the respective total gated FC charge.
TCanvas* PlotNH3NormalizedCounts( DataSet& AllData ){

    vector<int> palette = { 1, 632, 800, 400, 416, 600, 840, 880, 900, 616, 820, 432, 920, 860 };
    int pint = 0;

    TMultiGraph* mgNormCounts = new TMultiGraph();
    TLegend* mgLeg = new TLegend(0.1,0.7,0.4,0.9);
    //mgLeg->SetNColumns(3);
    mgLeg->SetHeader("Counts Normalized to NH3","C");

    // These vectors hold the total number normalized counts for each bin in Bjorken X.	
    vector<double> NH3_Xbins, ND3_Xbins, CH2_Xbins, C_Xbins, F_Xbins, ET_Xbins, CD2_Xbins, Xbins, Zeros;
    for(size_t i=0; i<X_Bin_Bounds.size()-1; i++){
	double xmid = (X_Bin_Bounds[i]+Q2_Bin_Bounds[i+1])/2.0;
        // Get the normalized counts for each bin
	double n_NH3 = 0; double n_ND3 = 0; double n_CH2 = 0; double n_C = 0; double n_F = 0; double n_ET = 0; double n_CD2 = 0;
	
	for(size_t j=0; j<Q2_Bin_Bounds.size()-1; j++){
	    double qmid = (Q2_Bin_Bounds[j]+Q2_Bin_Bounds[j+1])/2.0;
	    Bin thisBin = AllData.getThisBin( qmid, xmid ); // Get the current bin and bin contents
	    n_NH3 += thisBin.getNormNt("NH3");
	    n_ND3 += thisBin.getNormNt("ND3");
	    n_CH2 += thisBin.getNormNt("CH2");
	    n_C   += thisBin.getNormNt("C");
	    n_F   += thisBin.getNormNt("F");
	    n_ET  += thisBin.getNormNt("ET");
	    n_CD2 += thisBin.getNormNt("CD2");
	}
	if( n_NH3 > 0 ){
	    NH3_Xbins.push_back( n_NH3/n_NH3 ); // Should just be one, but this is a sanity check.
	    ND3_Xbins.push_back( n_ND3/n_NH3 ); CH2_Xbins.push_back( n_CH2/n_NH3 ); C_Xbins.push_back( n_C/n_NH3 ); 
	    F_Xbins.push_back( n_F/n_NH3 ); ET_Xbins.push_back( n_ET/n_NH3 ); CD2_Xbins.push_back( n_CD2/n_NH3 );
	    Xbins.push_back( xmid ); Zeros.push_back(0.0);
	}
    }
    if( Xbins.size() > 0 ){
	TGraphErrors* gr_NH3 = new TGraphErrors(Xbins.size(), Xbins.data(), NH3_Xbins.data(), Zeros.data(), Zeros.data() );
	gr_NH3->SetMarkerColor( palette[0] ); gr_NH3->SetMarkerStyle(kFullCircle);
	mgNormCounts->Add( gr_NH3, "p" ); mgLeg->AddEntry( gr_NH3, "NH3", "p" ); 
	//gr_NH3->GetYaxis()->SetLogy();
 	TGraphErrors* gr_ND3 = new TGraphErrors(Xbins.size(), Xbins.data(), ND3_Xbins.data(), Zeros.data(), Zeros.data() );
	gr_ND3->SetMarkerColor( palette[1] ); gr_ND3->SetMarkerStyle(kFullCircle);
	mgNormCounts->Add( gr_ND3, "p" ); mgLeg->AddEntry( gr_ND3, "ND3", "p" ); 
	//gr_ND3->GetYaxis()->SetLogy();
 	TGraphErrors* gr_CH2 = new TGraphErrors(Xbins.size(), Xbins.data(), CH2_Xbins.data(), Zeros.data(), Zeros.data() );
	gr_CH2->SetMarkerColor( palette[2] ); gr_CH2->SetMarkerStyle(kFullCircle);
	mgNormCounts->Add( gr_CH2, "p" ); mgLeg->AddEntry( gr_CH2, "CH2", "p" ); 
	//gr_CH2->GetYaxis()->SetLogy();
 	TGraphErrors* gr_C = new TGraphErrors(Xbins.size(), Xbins.data(), C_Xbins.data(), Zeros.data(), Zeros.data() );
	gr_C->SetMarkerColor( palette[3] ); gr_C->SetMarkerStyle(kFullCircle);
	mgNormCounts->Add( gr_C, "p" ); mgLeg->AddEntry( gr_C, "C", "p" ); 
	//gr_C->GetYaxis()->SetLogy();
 	TGraphErrors* gr_F = new TGraphErrors(Xbins.size(), Xbins.data(), F_Xbins.data(), Zeros.data(), Zeros.data() );
	gr_F->SetMarkerColor( palette[4] ); gr_F->SetMarkerStyle(kFullCircle);
	mgNormCounts->Add( gr_F, "p" ); mgLeg->AddEntry( gr_F, "F", "p" ); 
	//gr_F->GetYaxis()->SetLogy();
 	TGraphErrors* gr_ET = new TGraphErrors(Xbins.size(), Xbins.data(), ET_Xbins.data(), Zeros.data(), Zeros.data() );
	gr_ET->SetMarkerColor( palette[5] ); gr_ET->SetMarkerStyle(kFullCircle);
	mgNormCounts->Add( gr_ET, "p" ); mgLeg->AddEntry( gr_ET, "ET", "p" ); 
	//gr_ET->GetYaxis()->SetLogy();
 	TGraphErrors* gr_CD2 = new TGraphErrors(Xbins.size(), Xbins.data(), CD2_Xbins.data(), Zeros.data(), Zeros.data() );
	gr_CD2->SetMarkerColor( palette[6] ); gr_CD2->SetMarkerStyle(kFullCircle);
	mgNormCounts->Add( gr_CD2, "p" ); mgLeg->AddEntry( gr_CD2, "CD2", "p" ); 
	//gr_CD2->GetYaxis()->SetLogy();
	TCanvas* NH3_Norm_Plots = new TCanvas("NH3_Norm_Plots","NH3_Norm_Plots",800,600);
	NH3_Norm_Plots->cd();
	gPad->SetLogy();
	mgNormCounts->Draw("ap");
	mgLeg->Draw("same");
	return NH3_Norm_Plots;
    }
    else{
	cout << "Couldn't read in any counts; check inputs and try again.\n";
	TCanvas* nullPlt = new TCanvas();
	return nullPlt;
    }
}
*/

#endif
