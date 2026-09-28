
// This program calculates values for D and Eta using the new A1, R values

#include "../Header_Files/Binning_Classes.h"

using namespace std;

vector<int> palette = { 632, 800, 416, 600, 840, 880, 900, 616, 820, 432, 920, 860 };

void CalculateNewEtaA2(){

  TMultiGraph* mgR = new TMultiGraph();
  TLegend* legR = new TLegend(0.1,0.7,0.6,0.9); legR->SetNColumns(3);
  legR->SetHeader("Q^{2} Bin Midpoints [GeV^{2}]","C");

  TMultiGraph* mgEtaA2 = new TMultiGraph();
  TLegend* legEtaA2 = new TLegend(0.1,0.7,0.6,0.9); legEtaA2->SetNColumns(3);
  legEtaA2->SetHeader("Q^{2} Bin Midpoints [GeV^{2}]","C");

  TMultiGraph* mgDepol = new TMultiGraph();
  TLegend* legDepol = new TLegend(0.1,0.7,0.6,0.9); legDepol->SetNColumns(3);
  legDepol->SetHeader("Q^{2} Bin Midpoints [GeV^{2}]","C");

  // This is kinda a flimsy way of doing it, but I'll open and close the file for each Q2 bin
  for(size_t i=0; i<Q2_Bin_Bounds.size()-1; i++){

    double qmin = Q2_Bin_Bounds[i];
    double qmax = Q2_Bin_Bounds[i+1];

    ifstream fin("../Input_Text_Files/Asymmetry_Parameterizations.txt");

    string line;
    getline(fin,line); // Throw away header row

    vector<double> XVals, EtaA2Vals, DepolVals, RVals;

    // Loop through the file and calculate new values of D, Eta, etc.
    while( getline(fin, line) ){
	stringstream sin(line);
	double Q2, W2, x, F1, F2, R, A1, A2, g1, g2;
	sin >> Q2 >> W2 >> x >> F1 >> F2 >> R >> A1 >> A2 >> g1 >> g2;
	
	double Ep = avgBeamEnergy - ( (W2 + Q2 - (nucleon_mass*nucleon_mass)) / (2*nucleon_mass)   ); // Mean scattered energy for this Q2, W2
	double theta = 2*asin( sqrt( Q2/(4*avgBeamEnergy*Ep) ) ); // Average scattered angle in radians
	double tau = Q2 / (4*nucleon_mass*nucleon_mass*x*x); // tau kinematic variable
	double eps = 1.0 / ( 1.0 + 2*(1.0+tau)*pow( tan(theta/2.0),2)  ); // virtual photon polarization epsilon, with theta passed in with radians
	
	double depol = ( 1.0 - (eps*Ep/avgBeamEnergy) ) / ( 1.0 + (eps*R) ); // The depolarization factor
	double eta = (eps * sqrt(Q2) ) / ( avgBeamEnergy - (eps*Ep) ); // The eta kinematic variable
	
	// Loop thru x bins to find the correct x vals
	double xmin = 0;double xmax = 0;
	for(size_t j=0; j< X_Bin_Bounds.size()-1; j++){
	    if( x >= X_Bin_Bounds[j] && x < X_Bin_Bounds[j+1] ){
		xmin = X_Bin_Bounds[j];
		xmax = X_Bin_Bounds[j+1];
	    }
	}

	if( Q2 >= qmin && Q2 < qmax ){
	    double qmid = (qmin + qmax)/2.0;
	    double xmid = (xmin + xmax)/2.0;

	    XVals.push_back( x );
	    EtaA2Vals.push_back( eta*A2 );
	    DepolVals.push_back( depol );
	    RVals.push_back( R );
	    //fout << qmin <<"  "<< qmax <<"  "<< qmid <<"  "<< Q2 <<"  "<< xmin <<"  "<< xmax <<"  "<< xmid <<"  "<< x <<
	    //"  "<< W2 <<"  "<< F1 <<"  "<< F2 <<"  "<< R <<"  "<< A1 <<"  "<< A2 <<"  "<< g1 <<"  "<< g2 <<"  "<< depol <<"  "<< eta << endl;
	}
    }
    fin.close();

    stringstream bl; bl << "Q^{2} = " << setprecision(4) << (qmin+qmax)/2.0;
    string binLabel = bl.str();

    // Make the plots for this Q2 bin
    TGraph* grR = new TGraph( XVals.size(), XVals.data(), RVals.data() );
    grR->SetMarkerColor( palette[i] ); grR->SetMarkerStyle( kFullCircle );
    mgR->Add( grR, "p" );
    legR->AddEntry( grR, binLabel.c_str(), "p" );

    TGraph* grEtaA2 = new TGraph( XVals.size(), XVals.data(), EtaA2Vals.data() );
    grEtaA2->SetMarkerColor( palette[i] ); grEtaA2->SetMarkerStyle( kFullTriangleUp );
    mgEtaA2->Add( grEtaA2, "p" );
    legEtaA2->AddEntry( grEtaA2, binLabel.c_str(), "p" );

    TGraph* grDepol = new TGraph( XVals.size(), XVals.data(), DepolVals.data() );
    grDepol->SetMarkerColor( palette[i] ); grDepol->SetMarkerStyle( kFullSquare );
    mgDepol->Add( grDepol, "p" );
    legDepol->AddEntry( grDepol, binLabel.c_str(), "p" );

  }

  // Now make the plots
  TCanvas* c = new TCanvas("c","c",2000,600); c->Divide(3,1);
  c->cd(3);
  mgR->SetTitle("Model Values of R; x; R");
  mgR->GetXaxis()->SetLimits(0.1,0.8);
  mgR->GetYaxis()->SetRangeUser(0.0,0.4);
  mgR->Draw("ap");
  legR->Draw("same");

  c->cd(1);
  mgEtaA2->SetTitle("Model Values of #eta A_{2}; x; #eta A_{2}");
  mgEtaA2->GetXaxis()->SetLimits(0.1,0.8);
  mgEtaA2->GetYaxis()->SetRangeUser(0.0,0.05);
  mgEtaA2->Draw("ap");
  legEtaA2->Draw("same");

  c->cd(2);
  mgDepol->SetTitle("Model Values of Depolarization Factor; x; D");
  mgDepol->GetXaxis()->SetLimits(0.1,0.8);
  mgDepol->GetYaxis()->SetRangeUser(0.0,1.2);
  mgDepol->Draw("ap");
  legDepol->Draw("same");

}
