
#include "../Header_Files/Binning_Classes.h"

using namespace std;

const int nBits = 64;

void PrintSetBits( bitset<21>& mask, TH1D* hist ){

    for(int i=0; i<nBits; i++){
	if( mask[i] == 1 ){
	    //cout <<"---> Bit "<< i <<" is set!\n";
	    hist->Fill( i );
	    if( i >= 19 ) cout <<"---> Bit "<< i <<" is set!\n";
	}
    }

}

//void test(int run_start, int run_end){
void test(){


    vector<double> Angles = {-1, 0, 7.5, 8, 9, 11, 13.5, 17 };
    vector<int> Sectors = {-1, 7, 1, 2, 3, 4, 5, 6};
    vector<int> Runs = {16700, 17100, 17400, 17600};

    int it = Angles.size();
    int rt = Runs.size();

    for(int r=0; r<rt; r++){
	int Run = Runs[r];
        for(int i=0; i<it; i++){
	    cout <<"Before P_elec = "<< 5 <<", After P_corr = "<< MomentumCorrection( 5, Angles[i], Sectors[i], Run ) << endl;
        }
    }

/* Old stuff for testing error bits
  RunPeriod RP; // Holds run information

  //vector<int> Runs = RP.getDataSet( period, "Neg", "Neg" );
  vector<int> Runs = {16128, 16134, 16137, 16138, 16170, 16178};

  // Histogram to track which bits have been set...
  TH1D* BitTracker = new TH1D("BitTracker", string("Error Bits for Runs in "+ period +"; Error Bits;").c_str(), 100, -0.625, 24.375 );

  //while( run <= run_end ){
  for( int run : Runs ){

    string target = RP.getTargetType( run );

    cout << "Starting run "<< run << endl;

    TFile* fin = TFile::Open( string("/volatile/clas12/holmberg/Summer_22_Data/Electron_"+target+"_"+to_string(run)+"_Data.root").c_str(), "READ");
    if( fin->IsZombie() ){ run++; continue; }
    //if( !fin ){ run++; continue; }

    TTreeReader reader("Particle", fin);

    TTreeReaderValue<int> bitMask(reader, "qadb_mask_s1");
    TTreeReaderValue<int> bitMask2(reader, "qadb_mask_s2");

    bitset<21> SetBits(0b011110000110000011111);
    int valBits = SetBits.to_ulong();

    int nEvents = 0;

    int prevBit_s1 = 0;
    int prevBit_s2 = 0;

    while( reader.Next() ){
	if( *bitMask != 0 ){

	    bitset<21> binaryMask( *bitMask );
	    bitset<21> binaryMask2( *bitMask2 );
	    if( *bitMask != prevBit_s1 or *bitMask2 != prevBit_s2 ){
		//if( (binaryMask & SetBits) != 0 ) cout <<"Reject bin with errors: "<< binaryMask << endl;
		if( (*bitMask & valBits) != 0 ) cout <<"Reject bin with errors: "<< binaryMask << endl;

	        //cout << "Sector 1 bit mask (int): "<< *bitMask << endl;
	        //cout << "Sector 1 bit mask (bit): "<< binaryMask << endl;
	        //
		PrintSetBits( binaryMask, BitTracker );

	        //cout << "Sector 2 bit mask (int): "<< *bitMask2 << endl;
	        //cout << "Sector 2 bit mask (bit): "<< binaryMask2 << endl;
		//PrintSetBits( binaryMask2 );

		prevBit_s1 = *bitMask;
		prevBit_s2 = *bitMask2;
	    }
	    
	    nEvents++;
	}
    }
    fin->Close();

    cout << "Finished run "<< run <<"!\n";
  } // End of while loop


  TCanvas* c = new TCanvas("c","c",800,600);
  c->cd();
  BitTracker->SetFillColor(kRed);
  BitTracker->Draw("hist");
*/
}

