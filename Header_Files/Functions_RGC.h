/************************************************************************************************
 *
 * Author: Derek Holmberg
 *
 * Date Created: 8/29/25
 *
 * Last Modified: 9/18/26
 *
 * The purpose of this program is to make fiducial plots of the DC and ECAL for RGC data in CLAS12.
 * This contains functions and classes to make looping over input HIPO files easier. I moved this
 * information to a separate header file, since I need to make separate programs for Monte Carlo (MC)
 * data and regular HIPO data (the MC files lack certain data banks, which causes issues with my file
 * loop).
 *
************************************************************************************************/


#include "QADB.h"
#include "/w/hallb-scshelf2102/clas12/holmberg/Spin_Structure/Header_Files/Particle_Masses.h"

#ifndef FUNCTIONS_RGC_H
#define FUNCTIONS_RGC_H

string THISDIR = "/w/hallb-scshelf2102/clas12/holmberg/Spin_Structure/Header_Files/";

using namespace std;
using namespace QA;

// This function calculates the energy based on a particle's PID and momentum
double Energy(int pid, double p2){
	int id = abs(pid);
	if(id==2212){return sqrt(p2 + p_mass*p_mass);}
	else if(id==2112){return sqrt(p2 + n_mass*n_mass);}
	else if(id==321 ){return sqrt(p2 + k_mass*k_mass);}
	else if(id==211 ){return sqrt(p2 + ppm_mass*ppm_mass);}
	else if(id==11 || id==13 ){return sqrt(p2);} // Ignore electron and muon mass
	else if(id==22  ){return sqrt(p2);}
	else if(id==45  ){return sqrt(p2 + d_mass*d_mass);}
	return 0.0;
}

// Returns angle in degrees
double getPhi(double px, double py){
    if( px > 0 && py > 0 ) return 180.0 * atan(py / px) / M_PI;       // first quadrant
    else if( px<0 && py>0 ) return 180.0+ 180.0*atan(py / px) / M_PI;// second quadrant
    else if( px<0 && py<0 ) return 180.0+ 180.0*atan(py / px) / M_PI; // third quadrant
    else if( px>0 && py<0 ) return 360.0+ 180.0*atan(py / px) / M_PI;// fourth quadrant
    else{ // either px or py is zero 
	//std::cout << "Crap dammit there was a zero! I put the angle at -20000 degrees for error purposes lol.\n";
	return -20000;
    }
}
// Returns the phi angle local to the DC sector
double getLocalPhi(double x, double y, int sector){
    double phi = getPhi(x,y);
    if( sector == 1 ){
	if( phi > 180 ) return phi - 360;
	else return phi;
    }
    else return phi - 60*(sector-1);
}


// Holds information for calculating various kinematic quantities.
// Electron mass is ignored for all calculations in this program.
TLorentzVector BeamElectron1 = TLorentzVector( 0, 0, beam_energy1, beam_energy1 ); // Neglects electron mass
TLorentzVector BeamElectron2 = TLorentzVector( 0, 0, beam_energy2, beam_energy2 ); // Neglects electron mass
TLorentzVector BeamElectron3 = TLorentzVector( 0, 0, beam_energy3, beam_energy3 ); // Neglects electron mass

TLorentzVector TargetProton = TLorentzVector( 0, 0, 0, nucleon_mass );


// Simple class to hold some information from programs more easily
class Particle{

  public:
	int Run = 0; // Run number of this particular particle
	int PID = 0; // Monte-carlo ID of particle
	int Status = 0; // Location of hit within CLAS12 geometry
	double Chi2PID = 0; // How closely PID matches time-of-flight information based on particle momentum
	int Helicity = 0; // HWP and target polarization-corrected helicity state of particle
	int Charge = 0; // Charge of the particle
	double Beta = 0; // velocity over speed of light, taken directly from the REC::Particle bank
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
	// Kinematic quantities that can be calculated
	double E = 0; // Total energy of particle (electrons, muons neglect mass in this calculation)
	double W = 0; // Missing mass of event
	double Theta = 0; // Polar scattering angle IN DEGREES!!
	double Phi = 0; // Azimuthal scattering angle of event IN DEGREES!!
	double XposDC[3] = {0,0,0}; // X, Y, Z positions in each of the three regions of DC
	double YposDC[3] = {0,0,0}; // = {R1, R2, R3}
	double ZposDC[3] = {0,0,0};
	double LocalTheta[3] = {0,0,0}; // Polar angle of particle track relative to sector geometry (for each region R1, R2, R3)
	double LocalPhi[3] = {0,0,0}; // Azumuthal angle relative to sector geometry (for each region R1, R2, R3)
	double Q2 = 0; // Four-momentum transfer of the event
	double X = 0; // Bjorken X of the event
	double SF = 0; // Sampling fraction of this particle
	double P = 0; // Momentum of particle
	double Y = 0; // Energy transfer divided by beam energy
	double RasterX = 0; double RasterY = 0; // The raster position for this given particle
	double M2_PCAL = 0; // Second moment of the shower in the PCAL
	double M2_ECIN = 0; // Second moment of the shower in the ECIN
	double M2_EOUT = 0; // Second moment of the shower in the EOUT

	// Mutator
	void SetData( stringstream& sin ){
	    //stringstream sin(inputString);
	    sin >> PID >> Status >> Chi2PID >> Helicity >> Charge >> Px >> Py >> Pz >> Vz >> TrackChi2NDF >> EdgeR1 >> EdgeR2
		>> EdgeR3 >> SectorDC >> SectorECAL >> SectorHTCC
		>> PcalX >> PcalY >> PcalE >> PcalLu >> PcalLv >> PcalLw
		>> EcinX >> EcinY >> EcinE >> EcinLu >> EcinLv >> EcinLw
		>> EoutX >> EoutY >> EoutE >> EoutLu >> EoutLv >> EoutLw >> nphe
		>> XposDC[0] >> YposDC[0] >> ZposDC[0]
		>> XposDC[1] >> YposDC[1] >> ZposDC[1]
		>> XposDC[2] >> YposDC[2] >> ZposDC[2]
		>> Vx >> Vy >> RasterX >> RasterY >> Run >> Beta
		>> M2_PCAL >> M2_ECIN >> M2_EOUT;
	    // Now calculate the kinematic terms
	    E = Energy( PID, (Px*Px + Py*Py + Pz*Pz) ); // Electron, muon mass ignored
	    TLorentzVector ScatteredE = TLorentzVector( Px, Py, Pz, E );
	    TLorentzVector VirtualPhoton;
	    if( Run >= 16128 && Run <= 17067 ){      VirtualPhoton = BeamElectron1 - ScatteredE; Y = (beam_energy1 - E) / beam_energy1; }
	    else if( Run >= 17067 && Run <= 17704 ){ VirtualPhoton = BeamElectron2 - ScatteredE; Y = (beam_energy2 - E) / beam_energy2; }
	    else if( Run >= 17720 && Run <= 17811 ){ VirtualPhoton = BeamElectron3 - ScatteredE; Y = (beam_energy3 - E) / beam_energy3; }
	    else{ VirtualPhoton = BeamElectron1 - ScatteredE; Y = (beam_energy1 - E) / beam_energy1; } // Default to this for invalid runs or MC data for now...

	    // If you're looking at any other particle besides the trigger electron,
	    // then W, Q2, and X are going to be meaningless...
	    W = sqrt( (TargetProton+VirtualPhoton)*(TargetProton+VirtualPhoton) );
	    Q2 = abs( VirtualPhoton*VirtualPhoton );
	    X = Q2 / (2.0*TargetProton*VirtualPhoton);

	    Theta = acos( Pz/sqrt(Px*Px + Py*Py + Pz*Pz)) * 180.0/M_PI;
	    Phi = getPhi( Px, Py );
	    //LocalPhi = Phi - 60*(SectorDC -1);
	    SF = (PcalE + EcinE + EoutE)/sqrt(Px*Px + Py*Py + Pz*Pz);
	    Edep = PcalE + EcinE + EoutE;
	    P = sqrt(Px*Px + Py*Py + Pz*Pz);
	    //Y = (beam_energy - E) / beam_energy;
	   
	    if( XposDC[0] != 0 && YposDC[0] != 0 && ZposDC[0] != 0 && SectorDC != 0 ){
	        LocalTheta[0] = atan( sqrt( XposDC[0]*XposDC[0] + YposDC[0]*YposDC[0] )/ ZposDC[0] )*180.0/M_PI;
	        LocalPhi[0] = getLocalPhi( XposDC[0], YposDC[0], SectorDC );// - 360 + 60*(SectorDC-1);//atan( XposDC[0] / YposDC[0] )*180.0/M_PI;
	    }
	    if( XposDC[1] != 0 && YposDC[1] != 0 && ZposDC[1] != 0 && SectorDC != 0 ){
	        LocalTheta[1] = atan( sqrt( XposDC[1]*XposDC[1] + YposDC[1]*YposDC[1] )/ ZposDC[1] )*180.0/M_PI;
	        LocalPhi[1] = getLocalPhi( XposDC[1], YposDC[1], SectorDC );// - 360 + 60*(SectorDC-1);// atan( XposDC[1] / YposDC[1] )*180.0/M_PI;
	    }
	    if( XposDC[2] != 0 && YposDC[2] != 0 && ZposDC[2] != 0 && SectorDC != 0 ){
	        LocalTheta[2] = atan( sqrt( XposDC[2]*XposDC[2] + YposDC[2]*YposDC[2] )/ ZposDC[2] )*180.0/M_PI;
	        LocalPhi[2] = getLocalPhi( XposDC[2], YposDC[2], SectorDC );// - 360 + 60*(SectorDC-1);// atan( XposDC[2] / YposDC[2] )*180.0/M_PI;
	    }
	}

	void Print(){
	    cout << PID <<" "<< Status <<" "<< Chi2PID <<" "<< Helicity <<" "<< Charge <<" "<< Px <<" "<< Py <<" "<< Pz <<" "<< Vz <<" "<< TrackChi2NDF <<" "<< EdgeR1 <<" "<< EdgeR2
	    <<" "<< EdgeR3 <<" "<< SectorDC <<" "<< SectorECAL
	    <<" "<< PcalX <<" "<< PcalY <<" "<< PcalE <<" "<< PcalLu <<" "<< PcalLv <<" "<< PcalLw
	    <<" "<< EcinX <<" "<< EcinY <<" "<< EcinE <<" "<< EcinLu <<" "<< EcinLv <<" "<< EcinLw
	    <<" "<< EoutX <<" "<< EoutY <<" "<< EoutE <<" "<< EoutLu <<" "<< EoutLv <<" "<< EoutLw <<" "<< nphe << endl;
	    cout << "E = "<< E << endl;
	    cout << "W = "<< W << endl;
	    cout << "Q2= "<< Q2 << endl;
	    cout << "X = "<< X << endl;
	    cout << "Theta = "<< Theta << endl;
	    cout << "Phi = "<< Phi << endl;
	}
	// This function takes a filestream object and writes the minimum amount of information to the file
	void Write( ofstream& fout ){
	    fout << PID <<" "<< Chi2PID <<" "<< Helicity <<" "<< Px <<" "<< Py <<" "<< Pz 
	    <<" "<< TrackChi2NDF <<" "<< EdgeR1 <<" "<< EdgeR2
	    <<" "<< EdgeR3 <<" "<< SectorDC <<" "<< nphe;
	    if( PcalE == 0 && EcinE == 0 && EoutE == 0) fout << endl;
	    else if( PcalE != 0 && EcinE == 0 && EoutE == 0) fout <<" "<< PcalE <<" "<< PcalLu <<" "<< PcalLv <<" "<< PcalLw << endl;
	    else if( PcalE != 0 && EcinE != 0 && EoutE == 0){ 
		fout <<" "<< PcalE <<" "<< PcalLu <<" "<< PcalLv <<" "<< PcalLw
	             <<" "<< EcinE <<" "<< EcinLu <<" "<< EcinLv <<" "<< EcinLw << endl;
	    }
	    else{
	        fout <<" "<< PcalE <<" "<< PcalLu <<" "<< PcalLv <<" "<< PcalLw
	        <<" "<< EcinE <<" "<< EcinLu <<" "<< EcinLv <<" "<< EcinLw
	        <<" "<< EoutE <<" "<< EoutLu <<" "<< EoutLv <<" "<< EoutLw  << endl;
	    }
	}
	void Reset(){ // Resets the values of the particle to zero if there's an invalid input from the HIPO file
	    PID = 0; // Monte-carlo ID of particle
	    Status = 0; // Location of hit within CLAS12 geometry
	    Chi2PID = 0; // How closely PID matches time-of-flight information based on particle momentum
	    Helicity = 0; // HWP and target polarization-corrected helicity state of particle
	    Charge = 0; // Charge of the particle
	    Beta = 0;
	    Px = 0;     Py = 0;     Pz = 0; // Particle momenta
	    Vx = 0;	Vy = 0;     Vz = 0; // Vertex location in CLAS12 geometry
	    TrackChi2NDF = 0; // Quality of track reconstruction in the DC
	    EdgeR1 = 0;     EdgeR2 = 0;     EdgeR3 = 0; // Edge variable in each of the three regions of DC
	    SectorDC = 0;     SectorECAL = 0;     SectorHTCC = 0; // Location of track/hit in DC, ECAL, and HTCC, respectively (should be equal)
	    PcalX = 0;     PcalY = 0; // Hit locations in the Pcal
	    PcalE = 0; // Energy deposited in Pcal
	    PcalLu = 0;     PcalLv = 0;     PcalLw = 0; // Lu, Lv, Lw hit positions in the Pcal
	    EcinX = 0;     EcinY = 0; // Hit locations in the Ecin
	    EcinE = 0; // Energy deposited in Ecin
	    EcinLu = 0;     EcinLv = 0;     EcinLw = 0; // Lu, Lv, Lw hit positions in the Ecin
	    EoutX = 0;     EoutY = 0; // Hit locations in the Eout
	    EoutE = 0; // Energy deposited in Eout
	    EoutLu = 0;     EoutLv = 0;     EoutLw = 0; // Lu, Lv, Lw hit positions in the Eout
	    Edep = 0; // Total energy deposited    o the ECAL system
	    nphe = 0; // Number of photoelectrons produced in the HTCC
	    // Kinematic quantities that can be calculated
	    E = 0; // Total energy of particle (electrons, muons neglect mass in this calculation)
	    W = 0; // Missing mass of event
	    Theta = 0; // Polar scattering angle 
	    Phi = 0; // Azimuthal scattering angle of event
	    XposDC[0] = 0; XposDC[1] = 0; XposDC[2] = 0; // X, Y, Z positions in each of the three regions of DC
	    YposDC[0] = 0; YposDC[1] = 0; YposDC[2] = 0; // X, Y, Z positions in each of the three regions of DC
	    ZposDC[0] = 0; ZposDC[1] = 0; ZposDC[2] = 0; // X, Y, Z positions in each of the three regions of DC
	    LocalTheta[0] = 0; LocalTheta[1] = 0; LocalTheta[2] = 0; // X, Y, Z positions in each of the three regions of DC
	    LocalPhi[0] = 0; LocalPhi[1] = 0; LocalPhi[2] = 0; // X, Y, Z positions in each of the three regions of DC
	    Q2 = 0; // Four-momentum transfer of the event
	    X = 0; // Bjorken X of the event
	    SF = 0; // Sampling fraction of this particle
	    P = 0; // Momentum of particle
	    Y = 0; // Energy transfer divided by beam energy
	    RasterX = 0; RasterY = 0; // Beam raster positions
	}
	// These member functions check if the given particle passes given cuts
	bool IsTriggerElectron() const{
	    return Status < -2000 && Status > -4000 && PID == 11 /*&& abs(Helicity) == 1*/;
	}
	bool IsTriggerPositron() const{
	    return Status < -2000 && Status > -4000 && PID == -11 /*&& abs(Helicity) == 1*/;
	}
	bool IsValidSectorDC() const{
	    return SectorDC >= 1 && SectorDC <= 6;
	}
	bool IsValidSectorECAL() const{
	    return SectorECAL >= 1 && SectorECAL <= 6;
	}
	bool IsValidSectorHTCC() const{
	    return SectorHTCC >= 1 && SectorHTCC <= 6;
	}
	bool AllValidSectorHits() const{
	    return IsValidSectorHTCC() && IsValidSectorECAL() && IsValidSectorDC();
	}
	bool IsHitECAL(string ecal) const{
	    if( ecal == "Pcal" && PcalX != 0 && PcalY != 0 ) return true;
	    else if( ecal == "Ecin" && EcinX != 0 && EcinY != 0 ) return true;
    	    else if( ecal == "Eout" && EoutX != 0 && EoutY != 0 ) return true;
	    return false;
	}
	bool IsValidEnergyECAL() const{
	    return (PcalE + EcinE + EoutE) > 0;
	}
	bool IsValidTrackChi2() const{
	    return TrackChi2NDF > 0;
	}
	bool IsValidDCEdge(string region) const{
	    if( region == "R1" && EdgeR1 > 0 ) return true;
	    else if( region == "R2" && EdgeR2 > 0 ) return true;
	    else if( region == "R3" && EdgeR3 > 0 ) return true;
	    else return false;
	}
	bool BasicEdgeCut(string region) const{
	    if( region == "R1" && EdgeR1 > 3 ) return true;
	    else if( region == "R2" && EdgeR2 > 3 ) return true;
	    else if( region == "R3" && EdgeR3 > 10 ) return true;
	    else return false;
	}
	bool PassDCFiducialCut() const{
	    return P > 2 && PID == 11 && Status < -2000 && Status >= -4000 && abs(Helicity) == 1;
	}
	bool IsInForwardDetector() const{ // Makes sure there was no CD hit
	    return abs(Status) < 4000 && abs(Status) > 2000;
	}
	// Checks if the particle is a good ECAL hit
	bool IsGoodECAL() const{
	    return PcalE>0 && EcinE>0 && EoutE>0 && IsValidSectorECAL();
	}
	// Checks if the particle is a good DC track (goes through all three sectors)
	bool IsGoodDC() const{
	    return EdgeR1 > 0 && EdgeR2 > 0 && EdgeR3 > 0 && IsValidSectorDC();
	}


};


// This is a pretty hefty function that takes the input banks from a HIPO file and writes the contents of the first particle from the REC::Particle bank to 
// the Particle class object
//vector<Particle> SetParticleData( int PID, QADB* qa, const int corr_factor, bank& CONF, bank& PART, bank& HEL, bank& TRACK, bank& TRAJ, bank& CAL, bank& CHKV ){
vector<Particle> SetParticleData( int PID, QADB* qa, const int corr_factor, bank& CONF, bank& PART, bank& HEL, bank& TRACK, bank& TRAJ, bank& CAL, bank& CHKV, bank& RASTER ){

        int evnum = CONF.getInt("event",0);
        int runnum= CONF.getInt("run",0); // Already read in, but this is a double-check

    vector<Particle> Particles;

    // Require that the first particle in the bank be a trigger electron
    int trigPID = PART.getInt("pid",0);
    int trigStatus = PART.getInt("status",0);
    bool isTrigElec = trigPID == 11 && trigStatus > -4000 && trigStatus < -2000;

    //if(qa->Pass(runnum,evnum)) {
        // First, check the PID
        for(int i=0; i<PART.getRows(); i++){ //#1
	    Particle P;
	    int charge = PART.getInt("charge",i); // Gets the charge of the particle (1, -1, or 0)
            int pid = PART.getInt("pid",i);
            int status = PART.getInt("status",i);
	    bool pidCut = false; // Used for trigger electrons
	    if( PID == 11 ) pidCut = status < -2000 && status > -4000 && pid == 11;
	    //else if( PID == -11 ) pidCut = true; // this just returns the entire vector of particles...
	    else if( abs(PID) == 1 ) pidCut = PID == charge; // That is, if the particle PID is +/-1 for the charge of the particle to look at
	    else if( PID == 2 ) pidCut = isTrigElec; // This is for elastic ep scattering events; it assumes that this event has only one trigger electron and proton
	    else pidCut = pid == PID; // This just causes the program to effectively ignore this flag if it's not an electron the user's looking for
	    // This check is for a very specific case of getting all available electrons, pi+, and pi-

	    bool qaPass = qa->Pass(runnum,evnum); // Of course, ignore the qa if looking at MC data
            //if( pid == PID && qaPass ){ //#2
            if( pidCut && qaPass ){ //#2
                // For each particle, get its kinematic and detector geometry info
                double partP[3] = { PART.getFloat("px",i), PART.getFloat("py",i), PART.getFloat("pz",i) };
                double vx = PART.getFloat("vx",i); // the x vertex
                double vy = PART.getFloat("vy",i); // the y vertex
                double vz = PART.getFloat("vz",i); // the z vertex
		double beta = PART.getDouble("beta",i); // velocity divided by the speed of light
                int corrHelicity = HEL.getInt("helicity",i) * corr_factor; // Corrects the helicity state for HWP and target polarization for asymmetry
                double chi2pid = PART.getFloat("chi2pid",i); // How closely PID matches time-of-flight information based on particle momentum
		//int charge = PART.getInt("charge",i); // Gets the charge of the particle (1, -1, or 0)
		
                // Now loop through and get the relevant calorimeter information
                double pcalXY[2] = {0,0}; double ecinXY[2] = {0,0}; double eoutXY[2] = {0,0}; // x,y positions for each section of the ECAL
                double pcalL[3] = {0,0,0}; // Lu, Lv, Lw in PCAL
                double ecinL[3] = {0,0,0};
                double eoutL[3] = {0,0,0};
                double calEnergy[3] = {0,0,0}; // energy deposited in PCAL, ECIN, EOUT
		double secondMomentsPCAL[3] = {0,0,0}; // m2u, m2v, m2w = second moments for each PCAL view
		double secondMomentsECIN[3] = {0,0,0}; // m2u, m2v, m2w = second moments for each ECIN view
		double secondMomentsEOUT[3] = {0,0,0}; // m2u, m2v, m2w = second moments for each EOUT view

                int sectorECAL = 0; // Tracks the ECAL sector
		
                for(int j=0; j<CAL.getRows(); j++){
                    //cout <<" Check and make sure this is the correct particle for the loop"<<endl;
                    if( CAL.getInt("pindex",j) == i ){
                        int layer = CAL.getInt("layer",j);
                        if( layer == 1 ){ // PCAL
                            pcalXY[0] = CAL.getFloat("x",j); pcalXY[1] = CAL.getFloat("y",j);
                            calEnergy[0] = CAL.getFloat("energy",j); sectorECAL = CAL.getInt("sector",j);
                            pcalL[0] = CAL.getFloat("lu",j); pcalL[1] = CAL.getFloat("lv",j); pcalL[2] = CAL.getFloat("lw",j);
			    secondMomentsPCAL[0] = CAL.getFloat("m2u",j);
			    secondMomentsPCAL[1] = CAL.getFloat("m2v",j);
			    secondMomentsPCAL[2] = CAL.getFloat("m2w",j);
                        }
                        else if( layer == 4 ){ // ECIN
                            ecinXY[0] = CAL.getFloat("x",j); ecinXY[1] = CAL.getFloat("y",j);
                            calEnergy[1] = CAL.getFloat("energy",j); sectorECAL = CAL.getInt("sector",j);
                            ecinL[0] = CAL.getFloat("lu",j); ecinL[1] = CAL.getFloat("lv",j); ecinL[2]= CAL.getFloat("lw",j);
			    secondMomentsECIN[0] = CAL.getFloat("m2u",j);
			    secondMomentsECIN[1] = CAL.getFloat("m2v",j);
			    secondMomentsECIN[2] = CAL.getFloat("m2w",j);
                        }
                        else if( layer == 7 ){ // EOUT
                            eoutXY[0] = CAL.getFloat("x",j); eoutXY[1] = CAL.getFloat("y",j);
                            calEnergy[2] = CAL.getFloat("energy",j); sectorECAL = CAL.getInt("sector",j);
                            eoutL[0] = CAL.getFloat("lu",j); eoutL[1] = CAL.getFloat("lv",j); eoutL[2] = CAL.getFloat("lw",j);
			    secondMomentsEOUT[0] = CAL.getFloat("m2u",j);
			    secondMomentsEOUT[1] = CAL.getFloat("m2v",j);
			    secondMomentsEOUT[2] = CAL.getFloat("m2w",j);
                        }
                    }
                } // end of CAL bank loop
		
                // Now get the DC information
                double chi2 = 0; double ndf = 0;
                int sectorDC = 0;
                //cout <<" Loop through TRACK and get the reduced chi2 for the DC track"<<endl;
                for(int j=0; j<TRACK.getRows(); j++){
                    if( TRACK.getInt("pindex",j) == i ){
                        chi2 = TRACK.getFloat("chi2",j); ndf = 1.0*TRACK.getInt("NDF",j);
                        sectorDC = TRACK.getInt("sector",j);
                    }
                }
                //cout <<" Now go throught the TRAJ bank and get the edge variables for each region of DC\n";
                double edgeRegion[3] = {-1,-1,-1}; // Tracks edge for R1, R2, R3
		double X_POS_DC[3] = {0,0,0}; double Y_POS_DC[3] = {0,0,0}; double Z_POS_DC[3] = {0,0,0};
                for(int j=0; j<TRAJ.getRows(); j++){
                    if( TRAJ.getInt("pindex",j) == i ){
                        int detector = TRAJ.getInt("detector",j); // 6 corresponds to DC
                        int layer = TRAJ.getInt("layer",j); // 6 is R1, 18 is R2, 36 is R3
                        if( detector == 6 && layer == 6 ){
                            edgeRegion[0] = TRAJ.getFloat("edge",j);
			    X_POS_DC[0] = TRAJ.getFloat("x",j);
			    Y_POS_DC[0] = TRAJ.getFloat("y",j);
			    Z_POS_DC[0] = TRAJ.getFloat("z",j);
                        }
                        else if( detector == 6 && layer == 18 ){
                            edgeRegion[1] = TRAJ.getFloat("edge",j);
			    X_POS_DC[1] = TRAJ.getFloat("x",j);
			    Y_POS_DC[1] = TRAJ.getFloat("y",j);
			    Z_POS_DC[1] = TRAJ.getFloat("z",j);
                        }
                        else if( detector == 6 && layer == 36 ){
                            edgeRegion[2] = TRAJ.getFloat("edge",j);
			    X_POS_DC[2] = TRAJ.getFloat("x",j);
			    Y_POS_DC[2] = TRAJ.getFloat("y",j);
			    Z_POS_DC[2] = TRAJ.getFloat("z",j);
                        }
                    }
                }
                //cout <<" Now get the number of photoelectrons in the HTCC\n";
                double nphe = 0; int HTCCsector = 0;
                for(int j=0; j<CHKV.getRows(); j++){
                    if( CHKV.getInt("pindex",j) == i && CHKV.getInt("detector",j) == 15 ){
                        nphe = CHKV.getFloat("nphe",j);
			HTCCsector = CHKV.getInt("sector",j);
                    }
                }
		double rx = RASTER.getDouble("x",0); double ry = RASTER.getDouble("y",0);
		stringstream sout;
                if( ndf != 0 && abs(charge) == 1 ){
		    // Get the raster positions
		    //evCount++;
                    sout << pid<<" "<< status<<" "<< chi2pid<<" "<< corrHelicity<<" "<< charge <<" "<< partP[0]<<" "<< partP[1]<<" "<< partP[2]<<" "<< vz<<" "<< chi2/ndf<<" "<< edgeRegion[0]<<" "<< edgeRegion[1]<<" ";
                    sout<<edgeRegion[2]<<" "<<sectorDC<<" "<<sectorECAL<<" "<< HTCCsector <<" "<< pcalXY[0]<<" "<< pcalXY[1]<<" "<< calEnergy[0]<<" "<< pcalL[0]<<" "<< pcalL[1]<<" "<< pcalL[2]<<" ";
                    sout << ecinXY[0]<<" "<< ecinXY[1]<<" "<< calEnergy[1]<<" "<< ecinL[0]<<" "<< ecinL[1]<<" "<< ecinL[2]<<" ";
                    sout << eoutXY[0]<<" "<< eoutXY[1]<<" "<< calEnergy[2]<<" "<< eoutL[0]<<" "<< eoutL[1]<<" "<< eoutL[2]<<" "<< nphe <<" ";// << endl;
		    sout << X_POS_DC[0]<<" "<< Y_POS_DC[0]<<" "<< Z_POS_DC[0] <<" ";
		    sout << X_POS_DC[1]<<" "<< Y_POS_DC[1]<<" "<< Z_POS_DC[1] <<" ";
		    sout << X_POS_DC[2]<<" "<< Y_POS_DC[2]<<" "<< Z_POS_DC[2] <<" ";
		    sout << vx <<" "<< vy <<" "<< rx <<" "<< ry <<" "<< runnum <<" "<< beta <<" ";
		    sout << ( secondMomentsPCAL[0] + secondMomentsPCAL[1] + secondMomentsPCAL[2] ) / 3.0 <<" ";
		    sout << ( secondMomentsECIN[0] + secondMomentsECIN[1] + secondMomentsECIN[2] ) / 3.0 <<" ";
		    sout << ( secondMomentsEOUT[0] + secondMomentsEOUT[1] + secondMomentsEOUT[2] ) / 3.0 ;
		    P.SetData( sout );
		    Particles.push_back( P );
		} // end of if checking ndf and charge
/*
		else if( PID == 22 || PID == 2112 ){

                    sout << pid<<" "<< status<<" "<< chi2pid<<" "<< corrHelicity<<" "<< charge <<" "<< partP[0]<<" "<< partP[1]<<" "<< partP[2]<<" "<< vz<<" "<< 0 <<" "<< edgeRegion[0]<<" "<< edgeRegion[1]<<" ";
                    sout<<edgeRegion[2]<<" "<<sectorDC<<" "<<sectorECAL<<" "<< HTCCsector <<" "<< pcalXY[0]<<" "<< pcalXY[1]<<" "<< calEnergy[0]<<" "<< pcalL[0]<<" "<< pcalL[1]<<" "<< pcalL[2]<<" ";
                    sout << ecinXY[0]<<" "<< ecinXY[1]<<" "<< calEnergy[1]<<" "<< ecinL[0]<<" "<< ecinL[1]<<" "<< ecinL[2]<<" ";
                    sout << eoutXY[0]<<" "<< eoutXY[1]<<" "<< calEnergy[2]<<" "<< eoutL[0]<<" "<< eoutL[1]<<" "<< eoutL[2]<<" "<< nphe <<" ";// << endl;
		    sout << X_POS_DC[0]<<" "<< Y_POS_DC[0]<<" "<< Z_POS_DC[0] <<" ";
		    sout << X_POS_DC[1]<<" "<< Y_POS_DC[1]<<" "<< Z_POS_DC[1] <<" ";
		    sout << X_POS_DC[2]<<" "<< Y_POS_DC[2]<<" "<< Z_POS_DC[2] <<" ";
		    sout << vx <<" "<< vy <<" "<< rx <<" "<< ry <<" "<< runnum; // Passing in the zeros just keeps the raster position turned off
		    P.SetData( sout );
		    Particles.push_back( P );
		}
*/
		//else P->Reset(); // That way, the rest of the input loop knows it's an invalid read in
	    } //end of #2
	    //else P->Reset();
	} //end of #1
	
	// Before returning the particle, make sure it's a valid hit
	//if( !P->IsGoodECAL() || !P->IsGoodDC() ) P->Reset();

    //}

	//if( Particles.size() == 0 ) cout << "No particles in event.\n";	

	return Particles;
	

}


// This is a pretty hefty function that takes the input banks from a HIPO file and writes the contents of the first particle from the REC::Particle bank to 
// the Particle class object for MC data
//void SetParticleData( int PID, Particle* P, QADB* qa, const int corr_factor, bank& CONF, bank& PART, bank& HEL, bank& TRACK, bank& TRAJ, bank& CAL, bank& CHKV ){
vector<Particle> SetMCParticleData( int PID, bank& CONF, bank& PART, bank& TRACK, bank& TRAJ, bank& CAL, bank& CHKV ){

        int evnum = CONF.getInt("event",0);
        int runnum= CONF.getInt("run",0); // Already read in, but this is a double-check

    vector<Particle> Particles;

    // Require that the first particle in the bank be a trigger electron
    int trigPID = PART.getInt("pid",0);
    int trigStatus = PART.getInt("status",0);
    bool isTrigElec = trigPID == 11 && trigStatus > -4000 && trigStatus < -2000;

    //if(qa->Pass(runnum,evnum)) {
        // First, check the PID
        for(int i=0; i<PART.getRows(); i++){ //#1
	    Particle P;
	    int charge = PART.getInt("charge",i); // Gets the charge of the particle (1, -1, or 0)
            int pid = PART.getInt("pid",i);
            int status = PART.getInt("status",i);
	    bool pidCut = false; // Used for trigger electrons
	    if( PID == 11 ) pidCut = status < -2000 && status > -4000;
	    else if( abs(PID) == 1 ) pidCut = PID == charge;
	    else pidCut = pid == PID; // This just causes the program to effectively ignore this flag if it's not an electron the user's looking for
            if( /*pid == PID &&*/ pidCut /*&& isTrigElec*/ ){ //#2
                // For each particle, get its kinematic and detector geometry info
                double partP[3] = { PART.getFloat("px",i), PART.getFloat("py",i), PART.getFloat("pz",i) };
                double vx = PART.getFloat("vx",i); double vy = PART.getFloat("vy",i); double vz = PART.getFloat("vz",i); // the z vertex
                int corrHelicity = 1; // There's no helicity info for the MC data, so this is just set to 1 so I can resuse my data structures
		//int corrHelicity = HEL.getInt("helicity",i) * corr_factor; // Corrects the helicity state for HWP and target polarization for asymmetry
                double chi2pid = PART.getFloat("chi2pid",i); // How closely PID matches time-of-flight information based on particle momentum
		int charge = PART.getInt("charge",i); // Gets the charge of the particle (1, -1, or 0)
		
                // Now loop through and get the relevant calorimeter information
                double pcalXY[2] = {0,0}; double ecinXY[2] = {0,0}; double eoutXY[2] = {0,0}; // x,y positions for each section of the ECAL
                double pcalL[3] = {0,0,0}; // Lu, Lv, Lw in PCAL
                double ecinL[3] = {0,0,0};
                double eoutL[3] = {0,0,0};
                double calEnergy[3] = {0,0,0}; // energy deposited in PCAL, ECIN, EOUT
                int sectorECAL = 0; // Tracks the ECAL sector
		
                for(int j=0; j<CAL.getRows(); j++){
                    //cout <<" Check and make sure this is the correct particle for the loop"<<endl;
                    if( CAL.getInt("pindex",j) == i ){
                        int layer = CAL.getInt("layer",j);
                        if( layer == 1 ){ // PCAL
                            pcalXY[0] = CAL.getFloat("x",j); pcalXY[1] = CAL.getFloat("y",j);
                            calEnergy[0] = CAL.getFloat("energy",j); sectorECAL = CAL.getInt("sector",j);
                            pcalL[0] = CAL.getFloat("lu",j); pcalL[1] = CAL.getFloat("lv",j); pcalL[2] = CAL.getFloat("lw",j);
                        }
                        else if( layer == 4 ){ // ECIN
                            ecinXY[0] = CAL.getFloat("x",j); ecinXY[1] = CAL.getFloat("y",j);
                            calEnergy[1] = CAL.getFloat("energy",j); sectorECAL = CAL.getInt("sector",j);
                            ecinL[0] = CAL.getFloat("lu",j); ecinL[1] = CAL.getFloat("lv",j); ecinL[2]= CAL.getFloat("lw",j);
                        }
                        else if( layer == 7 ){ // EOUT
                            eoutXY[0] = CAL.getFloat("x",j); eoutXY[1] = CAL.getFloat("y",j);
                            calEnergy[2] = CAL.getFloat("energy",j); sectorECAL = CAL.getInt("sector",j);
                            eoutL[0] = CAL.getFloat("lu",j); eoutL[1] = CAL.getFloat("lv",j); eoutL[2] = CAL.getFloat("lw",j);
                        }
                    }
                } // end of CAL bank loop
		
                // Now get the DC information
                double chi2 = 0; double ndf = 0;
                int sectorDC = 0;
                //cout <<" Loop through TRACK and get the reduced chi2 for the DC track"<<endl;
                for(int j=0; j<TRACK.getRows(); j++){
                    if( TRACK.getInt("pindex",j) == i ){
                        chi2 = TRACK.getFloat("chi2",j); ndf = 1.0*TRACK.getInt("NDF",j);
                        sectorDC = TRACK.getInt("sector",j);
                    }
                }
                //cout <<" Now go throught the TRAJ bank and get the edge variables for each region of DC\n";
                double edgeRegion[3] = {-1,-1,-1}; // Tracks edge for R1, R2, R3
		double X_POS_DC[3] = {0,0,0}; double Y_POS_DC[3] = {0,0,0}; double Z_POS_DC[3] = {0,0,0};
                for(int j=0; j<TRAJ.getRows(); j++){
                    if( TRAJ.getInt("pindex",j) == i ){
                        int detector = TRAJ.getInt("detector",j); // 6 corresponds to DC
                        int layer = TRAJ.getInt("layer",j); // 6 is R1, 18 is R2, 36 is R3
                        if( detector == 6 && layer == 6 ){
                            edgeRegion[0] = TRAJ.getFloat("edge",j);
			    X_POS_DC[0] = TRAJ.getFloat("x",j);
			    Y_POS_DC[0] = TRAJ.getFloat("y",j);
			    Z_POS_DC[0] = TRAJ.getFloat("z",j);
                        }
                        else if( detector == 6 && layer == 18 ){
                            edgeRegion[1] = TRAJ.getFloat("edge",j);
			    X_POS_DC[1] = TRAJ.getFloat("x",j);
			    Y_POS_DC[1] = TRAJ.getFloat("y",j);
			    Z_POS_DC[1] = TRAJ.getFloat("z",j);
                        }
                        else if( detector == 6 && layer == 36 ){
                            edgeRegion[2] = TRAJ.getFloat("edge",j);
			    X_POS_DC[2] = TRAJ.getFloat("x",j);
			    Y_POS_DC[2] = TRAJ.getFloat("y",j);
			    Z_POS_DC[2] = TRAJ.getFloat("z",j);
                        }
                    }
                }
                //cout <<" Now get the number of photoelectrons in the HTCC\n";
                double nphe = 0; int HTCCsector = 0;
                for(int j=0; j<CHKV.getRows(); j++){
                    if( CHKV.getInt("pindex",j) == i && CHKV.getInt("detector",j) == 15 ){
                        nphe = CHKV.getFloat("nphe",j);
			HTCCsector = CHKV.getInt("sector",j);
                    }
                }
		stringstream sout;
                if( ndf != 0 && abs(charge) == 1 ){
		    //evCount++;
                    sout << pid<<" "<< status<<" "<< chi2pid<<" "<< corrHelicity<<" "<< charge <<" "<< partP[0]<<" "<< partP[1]<<" "<< partP[2]<<" "<< vz<<" "<< chi2/ndf<<" "<< edgeRegion[0]<<" "<< edgeRegion[1]<<" ";
                    sout<<edgeRegion[2]<<" "<<sectorDC<<" "<<sectorECAL<<" "<< HTCCsector <<" "<< pcalXY[0]<<" "<< pcalXY[1]<<" "<< calEnergy[0]<<" "<< pcalL[0]<<" "<< pcalL[1]<<" "<< pcalL[2]<<" ";
                    sout << ecinXY[0]<<" "<< ecinXY[1]<<" "<< calEnergy[1]<<" "<< ecinL[0]<<" "<< ecinL[1]<<" "<< ecinL[2]<<" ";
                    sout << eoutXY[0]<<" "<< eoutXY[1]<<" "<< calEnergy[2]<<" "<< eoutL[0]<<" "<< eoutL[1]<<" "<< eoutL[2]<<" "<< nphe <<" ";// << endl;
		    sout << X_POS_DC[0]<<" "<< Y_POS_DC[0]<<" "<< Z_POS_DC[0] <<" ";
		    sout << X_POS_DC[1]<<" "<< Y_POS_DC[1]<<" "<< Z_POS_DC[1] <<" ";
		    sout << X_POS_DC[2]<<" "<< Y_POS_DC[2]<<" "<< Z_POS_DC[2] <<" ";
		    sout << vx <<" "<< vy <<" "<< 0 <<" "<< 0 <<" "<< runnum; // Passing in the zeros just keeps the raster position turned off
		    P.SetData( sout );
		    Particles.push_back( P );
		} // end of if checking ndf and charge
/*
		else if( PID == 22 || PID == 2112 ){
                    sout << pid<<" "<< status<<" "<< chi2pid<<" "<< corrHelicity<<" "<< charge <<" "<< partP[0]<<" "<< partP[1]<<" "<< partP[2]<<" "<< vz<<" "<< 0 <<" "<< edgeRegion[0]<<" "<< edgeRegion[1]<<" ";
                    sout<<edgeRegion[2]<<" "<<sectorDC<<" "<<sectorECAL<<" "<< HTCCsector <<" "<< pcalXY[0]<<" "<< pcalXY[1]<<" "<< calEnergy[0]<<" "<< pcalL[0]<<" "<< pcalL[1]<<" "<< pcalL[2]<<" ";
                    sout << ecinXY[0]<<" "<< ecinXY[1]<<" "<< calEnergy[1]<<" "<< ecinL[0]<<" "<< ecinL[1]<<" "<< ecinL[2]<<" ";
                    sout << eoutXY[0]<<" "<< eoutXY[1]<<" "<< calEnergy[2]<<" "<< eoutL[0]<<" "<< eoutL[1]<<" "<< eoutL[2]<<" "<< nphe <<" ";// << endl;
		    sout << X_POS_DC[0]<<" "<< Y_POS_DC[0]<<" "<< Z_POS_DC[0] <<" ";
		    sout << X_POS_DC[1]<<" "<< Y_POS_DC[1]<<" "<< Z_POS_DC[1] <<" ";
		    sout << X_POS_DC[2]<<" "<< Y_POS_DC[2]<<" "<< Z_POS_DC[2] <<" ";
		    sout << vx <<" "<< vy <<" "<< 0 <<" "<< 0 <<" "<< runnum; // Passing in the zeros just keeps the raster position turned off
		    P.SetData( sout );
		    Particles.push_back( P );
		}
*/
		//else P->Reset(); // That way, the rest of the input loop knows it's an invalid read in
	    } //end of #2
	    //else P->Reset();
	} //end of #1
	
	// Before returning the particle, make sure it's a valid hit
	//if( !P->IsGoodECAL() || !P->IsGoodDC() ) P->Reset();

    //}
	
	return Particles;

}


#endif
