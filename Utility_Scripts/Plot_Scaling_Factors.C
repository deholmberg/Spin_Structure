#include "../Header_Files/Binning_Classes.h"
#include "../Header_Files/Background_Runs.h"

// This script just makes a plot of the scaling factors in the file "SummerAndFall22_CH2_Scaling_Factors.txt"

using namespace std;

TGraphErrors* Plot(string filePath, double& avgVal, double& avgValErr, string title="Scaling Factors", bool useRadLenCorr=false){

	// Holds the scaling factors, their statistical errors, and the corresponding x values.
	// Multiple scaling factors for each x because there's a Q2 dependence as well.
	vector<double> XBins, ScaleF, ScaleFerr, Zeros;

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
			Zeros.push_back( 0 );
			ScaleF.push_back( sf );
			ScaleFerr.push_back( sferr );
		}
	}
	fin.close();

	TGraphErrors* SF_Plt = new TGraphErrors( XBins.size(), XBins.data(), ScaleF.data(), Zeros.data(), ScaleFerr.data() );
	SF_Plt->SetTitle( title.c_str() );
	SF_Plt->SetMarkerStyle(kFullCircle);
	SF_Plt->SetMarkerColor(kViolet);
	// Calculate the weighted average, if applicable
	cout << "Results for "<< filePath <<":\n";
	TFitResultPtr r;
	r = SF_Plt->Fit("pol0","S","Q"); // Fit to a constant
	avgVal = r->Value(0); avgValErr = r->Error(0);
	//SF_Plt->GetYaxis()->SetRangeUser(0.8,1.2);

	return SF_Plt;

	//TCanvas* c = new TCanvas("c","c",800,600);
	//c->cd();
	//SF_Plt->Draw("ap");
	
}

// Plots the average of two scaling factor plots taken together
TGraphErrors* PlotAverage(string filePath1, string filePath2, double& avgVal, double& avgValErr, string title="Scaling Factors", bool useRadLenCorr=false){

	// Holds the scaling factors, their statistical errors, and the corresponding x values.
	// Multiple scaling factors for each x because there's a Q2 dependence as well.
	vector<double> XBins, ScaleF, ScaleFerr, Zeros;

	string line, line2; // Used for reading in data
	ifstream fin( filePath1 );
	ifstream fin2( filePath2 );
	ofstream fout( "../Input_Text_Files/Solenoid_Scaling_Factors.txt");
	// Kids, this is sloppy. Don't do this lol
	while( getline( fin, line ) && getline( fin2, line2 ) ){
		stringstream sin(line);
		double qmin, qmax, xmin, xmax, sf, sferr, qavg, xavg;
		sin >> qmin >> qmax >> xmin >> xmax >> sf >> sferr >> qavg >> xavg;
		stringstream sin2(line2);
		double qmin2, qmax2, xmin2, xmax2, sf2, sferr2, qavg2, xavg2;
		sin2 >> qmin2 >> qmax2 >> xmin2 >> xmax2 >> sf2 >> sferr2 >> qavg2 >> xavg2;

		if( useRadLenCorr ){
		    double corr = ET_Scale_Factor( xavg, qavg );
		    sf = sf*corr;
		    sferr = sferr*corr;
		    double corr2 = ET_Scale_Factor( xavg2, qavg2 );
		    sf2 = sf2*corr2;
		    sferr2 = sferr2*corr2;
		}

		double sfAvg = 0; double sfAvgErr = 0;
		if( sf != 0 && sf2 != 0 ){ // If it's zero, there's no data there
			//cout << xavg <<" "<< sf << endl;
			sfAvg = (sf + sf2)/2.0;
			sfAvgErr = 0.5*sqrt( sferr*sferr + sferr2*sferr2 );
			XBins.push_back( (xavg+xavg2)/2.0 ); // Not true stat-weighted value, but close enough lol
			Zeros.push_back( 0 );
			ScaleF.push_back( sfAvg );
			ScaleFerr.push_back( sfAvgErr );

			fout << qmin <<"     "<< qmax <<"     "<< xmin <<"     "<< xmax <<"     "<< sfAvg <<"     "<< sfAvgErr <<"     "<< (qavg+qavg2)/2.0 <<"     "<< (xavg+xavg2)/2.0 <<"     "<< endl;
		}
		else fout << qmin <<"     "<< qmax <<"     "<< xmin <<"     "<< xmax <<"     "<< sfAvg <<"     "<< sfAvgErr <<"     "<< 0 <<"     "<< 0 <<"     "<< endl;

	}
	fin.close();
	fin2.close();
	fout.close();

	TGraphErrors* SF_Plt = new TGraphErrors( XBins.size(), XBins.data(), ScaleF.data(), Zeros.data(), ScaleFerr.data() );
	SF_Plt->SetTitle( title.c_str() );
	SF_Plt->SetMarkerStyle(kFullCircle);
	SF_Plt->SetMarkerColor(kViolet);
	// Calculate the weighted average, if applicable
/*
	TFitResultPtr r;
	r = SF_Plt->Fit("pol0","S","Q"); // Fit to a constant
	avgVal = r->Value(0); avgValErr = r->Error(0);
*/
	//SF_Plt->GetYaxis()->SetRangeUser(0.8,1.2);

	return SF_Plt;

	//TCanvas* c = new TCanvas("c","c",800,600);
	//c->cd();
	//SF_Plt->Draw("ap");
	
}

void Plot_Scaling_Factors(){

    double avg, avgErr; // Used for plotting

    TGraphErrors* grCHCD = Plot("../Input_Text_Files/SuFa22_CH2_Scaling_Factors.txt", avg, avgErr,"Spring CD2/CH2 Scaling Factors; Bjorken X; n_{CD2}/n_{CH2}");
    grCHCD->SetMarkerColor(kOrange);
    grCHCD->GetYaxis()->SetRangeUser(0.8,1.2);
    
    TGraphErrors* grETC = Plot("../Input_Text_Files/Fall_ET_Scaling_Factors.txt", avg, avgErr,"Fall ET/C Scaling Factors; Bjorken X; n_{ET}/n_{C}");
    grETC->SetMarkerColor(kGreen);
    grETC->GetYaxis()->SetRangeUser(0.1,0.35);
/*
    TGraphErrors* grCH2 = Plot("../Input_Text_Files/Fall_CH2_Ratios.txt", avg, avgErr,"Fall CH2_{Pos}/CH2_{Neg} Ratio; Bjorken X; n_{CH2,pos}/n_{CH2,neg}");
    grCH2->SetMarkerColor(kRed);
    grCH2->GetYaxis()->SetRangeUser(0.8,1.2);

    string ch2line = to_string(avg);
    TF1* CH2Line = new TF1("CH2Line",ch2line.c_str(),0,1);
    CH2Line->SetLineColor(kRed);
*/
    TGraphErrors* grFC = Plot("../Input_Text_Files/Fall_F_Scaling_Factors.txt", avg, avgErr,"Fall F/C Scaling Factors; Bjorken X; n_{F}/n_{C}");
    grFC->GetYaxis()->SetRangeUser(0,0.03);
    grFC->SetMarkerColor(kCyan);

/*
    TGraphErrors* grCETC = Plot("../Input_Text_Files/Fall_ET_Scaling_Factors.txt", avg, avgErr,"Fall ET/C Corrected Scaling Factors; Bjorken X; C*n_{ET}/n_{C}", true);
    grCETC->SetMarkerColor(kGreen);
    grCETC->GetYaxis()->SetRangeUser(0.1,0.35);
*/
/*
    TGraphErrors* grC = Plot("../Input_Text_Files/Fall_C_Ratios.txt", avg, avgErr,"Fall C_{Pos}/C_{Neg} Ratio; Bjorken X; n_{C,pos}/n_{C,neg}");
    grC->SetMarkerColor(kBlue);
    grC->GetYaxis()->SetRangeUser(0.8,1.2);
    string cline = to_string(avg);
    TF1* CLine = new TF1("CLine",cline.c_str(),0,1);
    CLine->SetLineColor(kBlue);
*/
    TGraphErrors* avgScale = PlotAverage("../Input_Text_Files/Fall_CH2_Ratios.txt", "Fall_C_Ratios.txt", avg, avgErr,"Average Solenoid Scaling Ratio; Bjorken X; n_{Pos}/n_{Neg}");
    avgScale->SetMarkerColor(kViolet);
    avgScale->GetYaxis()->SetRangeUser(0.8,1.2);


    TCanvas* c = new TCanvas("c","c",1400,1000); c->Divide(3,2);
    c->cd(1); grCHCD->Draw("AP");
    c->cd(2); grETC->Draw("AP");
    c->cd(3); grFC->Draw("AP");
    c->cd(4); grCH2->Draw("AP"); CH2Line->Draw("same");
    c->cd(5); grC->Draw("AP"); CLine->Draw("same");
    c->cd(6); avgScale->Draw("AP");

    // Plot only the solenoid scaling factors...
    TCanvas* c2 = new TCanvas("c2","c2",1800,1000); c2->Divide(3,1);
    c2->cd(1); grC->Draw("AP");
    c2->cd(2); grCH2->Draw("AP");
    c2->cd(3); avgScale->Draw("AP");
//    c->cd(5); grCETC->Draw("AP");

}
