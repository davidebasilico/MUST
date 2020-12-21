#include <fstream>
#include <string>
#include <iostream>
#include <cstdlib>
#include <sstream>
#include <iomanip>

#include <TApplication.h>
#include <TVector3.h>
#include <TFile.h>
#include <TChain.h>
#include <TH1.h>
#include "TMath.h"
#include <TF1.h>
#include "TH1D.h"
#include "TH1.h"
#include "TH2F.h"
#include "TH3F.h"
#include <TROOT.h>
#include "TStyle.h"
#include <TCanvas.h>
#include <TRandom3.h>
#include <TSpectrum.h>
#include <TMultiGraph.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TLegend.h>
#include <TVirtualFitter.h>

#include <TRandom3.h>
#include <TFractionFitter.h>
#include <vector>

#include <TFitResult.h>
#include <TMatrixDSym.h>
#include "TMinuit.h"

#include <TString.h>

#ifndef ROOT_TH1D
#endif

#define PI 3.141592653589793

using namespace std;

int N = 13;
TH1D** hPDF = new TH1D*[N];
TH1D** hPDF_Tag = new TH1D*[N];

TH1D *PseudoDataset_Tag;
TH1D *PseudoDataset_Sub;

int Beg, End;

double density = 896.0; // [kg/m3]
double radius_FV = 14.0;
double Mass_FV = 4./3.*PI*pow(radius_FV,3)*density; // [m3]
double *Threshold = new double[N];
double *Events = new double[N];
double *Events_Tag = new double[N];
double *Events_Sub = new double[N];

double *InjRate = new double[N];
double *InjRate_Tag = new double[N];
double *InjRate_Sub = new double[N];
double *Fix = new double[N]();
double *Sigma = new double[N]();

double days;
double Exposure_LT;
double Efficiency;
double Exposure_LT_Tag;
double Exposure_LT_Sub;

double GetFloatPrecision(double value, double precision){
    return (floor((value * pow(10, precision) + 0.5)) / pow(10, precision));
}

double fitfunction(double *x, double *par){

	int bin = hPDF[0]->FindBin(x[0]);

	double *y = new double[N]();
	double sum=0;

	for(int i=0; i<N; i++){
		y[i] = par[i]* hPDF[i]->GetBinContent(bin);
		sum = sum+y[i];
		}

	return sum;
	}


double LogPoisson( double x, double par) {
    if (x < 0) return 0;
    else if (x == 0.) return -par;
    else return x*log(par)-par-TMath::LnGamma(x+1.); //altrimenti questo, LnGamma=ln(gamma), da geogebra torna che è sempre negativa se x>0
//    else return x*log(x/par)+par-x; //altrimenti questo, LnGamma=ln(gamma), da geogebra torna che è sempre negativa se x>0
}

void fcn(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag){

    Int_t i;
    Double_t Model_Tag=0, Model_Sub=0, Data_Tag=0, Data_Sub=0, Chi2=0, Penalty=0;

    for(int bin=Beg; bin<End; bin++){

    	Model_Tag = 0; Data_Tag=0;
    	Model_Sub = 0; Data_Sub=0;

	for(int i=0; i<N; i++){
		
		// Common species for Tag and Sub spectra (no C11, no C11_2, no C10_2, no He6_2)
		if(i!=9 && i!=10 && i!=11 && i!=12){
			Model_Tag += par[i]* hPDF[i]->GetBinContent(bin)*Exposure_LT_Tag/Exposure_LT;
			Model_Sub += par[i]* hPDF[i]->GetBinContent(bin)*Exposure_LT_Sub/Exposure_LT;
			}

                 // C11 for Sub spectrum
	         if(i==11){Model_Sub += par[i]* hPDF[i]->GetBinContent(bin)*Exposure_LT_Sub/Exposure_LT;}

                 // C10_2 or He_6  or C11_2
	         if(i==12 || i==9 || i==10){Model_Tag += par[i]* hPDF[i]->GetBinContent(bin)*Exposure_LT_Tag/Exposure_LT;}

	}

	Data_Tag += PseudoDataset_Tag->GetBinContent(bin);
	Data_Sub += PseudoDataset_Sub->GetBinContent(bin);

	if(Data_Tag!=0 && Data_Sub!=0) Chi2 += -2.*LogPoisson(Data_Tag,Model_Tag) -2.*LogPoisson(Data_Sub,Model_Sub);
	
	}

    for (int i=0; i<N; i++) {                                                                        // SETTING PENALTY
    	if (Sigma[i]==0.) Penalty+=0.;
	else Penalty+=pow((((par[i]/Exposure_LT)-InjRate[i])/Sigma[i]),2);
    }

    f = Chi2 + Penalty;
}


double MUST(string Configuration_Text, string Output_Rootfile, string Output_Text){

	//######################### A) READING INPUTS ############################################################################

	ifstream ReadCfgFile;
	ReadCfgFile.open(Configuration_Text.c_str());

	double TotalNumberEvents = 0;
	double TotalNumberEventsSub = 0;
	double TotalNumberEventsTag = 0;
	double DutyCycle;
	double TaggingPower;
	double Efficiency;
	string PDFs_Path;

	cout << "Cfg file: " << Configuration_Text.c_str() << endl;

	ReadCfgFile >> Beg;
	ReadCfgFile >> End;
	ReadCfgFile >> days;
	ReadCfgFile >> DutyCycle;
	ReadCfgFile >> TaggingPower;
	ReadCfgFile >> Efficiency;
	ReadCfgFile >> PDFs_Path;
	
	// Exposure expressed in 100 tons
	Exposure_LT = Mass_FV / pow(10,5) * days * DutyCycle;
	Exposure_LT_Tag = Exposure_LT*(1-Efficiency);
	Exposure_LT_Sub = Exposure_LT*Efficiency;

  	// Getting PDFs from PDFs rootfile (PDFs_Path)
	TFile *file_pdf = new TFile(PDFs_Path.c_str());
	hPDF[0] = (TH1D*)file_pdf->Get("be7_charge");
	hPDF[1] = (TH1D*)file_pdf->Get("pep_charge");
	hPDF[2] = (TH1D*)file_pdf->Get("bi210_charge");
	hPDF[3] = (TH1D*)file_pdf->Get("k40_charge");
	hPDF[4] = (TH1D*)file_pdf->Get("kr85_charge");
	hPDF[5] = (TH1D*)file_pdf->Get("u238_charge");
	hPDF[6] = (TH1D*)file_pdf->Get("th232_charge");
	hPDF[7] = (TH1D*)file_pdf->Get("po210_charge");
	hPDF[8] = (TH1D*)file_pdf->Get("cno_charge");
	hPDF[9] = (TH1D*)file_pdf->Get("c10_charge");
	hPDF[10] = (TH1D*)file_pdf->Get("he6_charge");
	hPDF[11] = (TH1D*)file_pdf->Get("c11_charge");
	hPDF[12] = (TH1D*)hPDF[11]->Clone();

  	int NEnergyBins = hPDF[0]->GetNbinsX();

	//Reading information from configuration file (Configuration_Text)
	for(int i=0; i<N; i++){
		ReadCfgFile >> InjRate[i];                             // Injected rate for each species
		ReadCfgFile >> Fix[i];                                 // Can be 0 (species is free or penalty'd) o 1 (species is fixed)
		ReadCfgFile >> Sigma[i];                               // If not zero, it introduces a gaussian penalty with the specified sigma
		Threshold[i]=hPDF[i]->Integral();
		InjRate_Tag[i]=InjRate[i];
		InjRate_Sub[i]=InjRate[i];
		}

	// Rates for Sub and Tag datasets: C11 and C11_2 , that is [11] and [12]
	InjRate_Sub[11] = InjRate[11]*(1.-TaggingPower);
	InjRate_Sub[12]= 0;
	InjRate_Tag[11] = 0;
	InjRate_Tag[12] = (InjRate[11]*Exposure_LT-InjRate_Sub[11]*Exposure_LT_Sub)/(Exposure_LT-Exposure_LT_Sub);

	InjRate[11]=InjRate_Sub[11];
	InjRate[12]=InjRate_Tag[12];

  	// Rates for Sub and Tag datasets: C10_2 and He6_2 , that is [9] and [10]
	InjRate_Sub[9]=0;
	InjRate_Sub[10]=0;
	InjRate_Tag[9]=InjRate[9];
	InjRate_Tag[10]=InjRate[10];

	//######################### B) GENERATING THE PSEUDO DATASETS ################################################################

 	// Definition of PseudoDatasets for Sub and Tag spectra
	PseudoDataset_Tag = new TH1D("PseudoDataset_Tag","PseudoDataset_Tag",NEnergyBins,0,NEnergyBins);
	PseudoDataset_Sub = new TH1D("PseudoDataset_Sub","PseudoDataset_Sub",NEnergyBins,0,NEnergyBins);

	gRandom = new TRandom3(0);
	gRandom->SetSeed(0);                                                             	// Sets a TRandom3 seed different for each running

	for(int i=0; i<N; i++){
		Events[i]  = InjRate[i] * Exposure_LT / Threshold[i];
		Events_Tag[i]  = InjRate_Tag[i] * Exposure_LT_Tag * Threshold[i];              // Calculating injected number of events in Tag spectrum for each species
		Events_Sub[i]  = InjRate_Sub[i] * Exposure_LT_Sub * Threshold[i];              // Calculating injected number of events in Sub spectrum for each species
		PseudoDataset_Tag->FillRandom(hPDF[i],Events_Tag[i]);                          // Randomly filling Tag spectrum
		PseudoDataset_Sub->FillRandom(hPDF[i],Events_Sub[i]);                          // Randomly filling Sub spectrum

		TotalNumberEvents += InjRate[i]/Threshold[i]*Exposure_LT;                      // Total number of events
		TotalNumberEventsSub += InjRate_Sub[i]/Threshold[i]*Exposure_LT_Sub;           // Total number of events for Sub spectrum
		TotalNumberEventsTag += InjRate_Tag[i]/Threshold[i]*Exposure_LT_Tag;           // Total number of events for Tag spectrum

	//cout << i << setw(12) << InjRate_Sub[i] << setw(12) << InjRate_Tag[i] << setw(12) << Events_Sub[i] << setw(12) << Events_Tag[i] << setw(12) << InjRate[i]/Threshold[i]*Exposure_LT << endl;
	}

	//######################### C) SPECTRAL FIT ####################################################################################


	TH1D *hPDF_Total = new TH1D("hPDF_Total","hPDF_Total",NEnergyBins,0,NEnergyBins);

	TMinuit *gMinuit = new TMinuit(N);
	gMinuit->SetFCN(fcn);

	Double_t arglist[2];
	Int_t ierflg = 0;
	TString vinit[13] = {"Be7", "pep", "Bi210", "K40", "Kr85", "U238", "Th232","Po210","CNO","C10_2","He6_2","C11","C11_2"};	// default species

 	for(int i=0; i<N; i++){
    		gMinuit->mnparm(i, vinit[i], 1.1*InjRate[i]/Threshold[i]*Exposure_LT, 1, 0,10.*InjRate[i]/Threshold[i]*Exposure_LT,ierflg);
		if (Fix[i] == 1 ) gMinuit->mnparm(i, vinit[i], 1.*InjRate[i]/Threshold[i]*Exposure_LT, 0, InjRate[i]/Threshold[i]*Exposure_LT,InjRate[i]/Threshold[i]*Exposure_LT,ierflg);
	}
	
	// Preparing  MIGRAD
	TVirtualFitter::SetMaxIterations(2e5);
	arglist[0] = 1;
	gMinuit->mnexcm("SET ERR",arglist,1,ierflg);
	arglist[0] = 5000;
	arglist[1] = 1;
	gMinuit->mnexcm("MIGRAD",arglist,2,ierflg);

 	// Print results
	Double_t amin,edm,errdef;
	Int_t nvpar,nparx,icstat;
	gMinuit->mnstat(amin,edm,errdef,nvpar,nparx,icstat);

	//######################### D) PREPARING PLOTS AND WRITING OUTPUTS ##############################################################

	double *Results = new double[N]();
	double *ErrorResults = new double[N]();

	for(int i=0; i<N; i++){
	     gMinuit->GetParameter(i,Results[i],ErrorResults[i]); //ritorna valore parametro e l'errore
	}

	TFile *foutput = new TFile (Output_Rootfile.c_str(), "RECREATE");
	foutput->cd();

	PseudoDataset_Tag->Write("PseudoDataset_Tag");
	PseudoDataset_Sub->Write("PseudoDataset_Sub");

	TH1D *hPDF_Total_Sub = new TH1D("hPDF_Total_Sub","hPDF_Total_Sub",NEnergyBins,0,NEnergyBins);
	TH1D *hPDF_Total_Tag = new TH1D("hPDF_Total_Tag","hPDF_Total_Tag",NEnergyBins,0,NEnergyBins);

	TCanvas *c = new TCanvas("c","c",1500,700);
	gPad->SetLogy();	gStyle->SetOptTitle(0); gStyle->SetOptStat(0);

	double yLimPlot = 0.3;

	TPad*	 Pad_UL = new TPad("Pad_UL", "Pad_UL", 0, yLimPlot, 0.5, 1.0);
	Pad_UL->Draw();
	Pad_UL->cd();
	PseudoDataset_Sub->Draw();

	double *model_Tag=new double [End-Beg];
	double *model_Sub=new double [End-Beg];
	double *data_Tag=new double [End-Beg];
	double *data_Sub=new double [End-Beg];
	double *phelectr=new double [End-Beg];
	double *y_Tag=new double [End-Beg];
	double *y_Sub=new double [End-Beg];

	for(int bin=0; bin<End-Beg; bin++) {
		double Residuals_model_Tag=0, Residuals_model_Sub=0, Residuals_data_Tag=0, Residuals_data_Sub=0;
		for(int i=0; i<N; i++){
	  	   if(i!=11)  Residuals_model_Tag += Results[i]* hPDF[i]->GetBinContent(bin+Beg)*Exposure_LT_Tag/Exposure_LT;
	   	  if(i!=9 && i!=10 && i!= 12) Residuals_model_Sub += Results[i]* hPDF[i]->GetBinContent(bin+Beg)*Exposure_LT_Sub/Exposure_LT;
	}

	Residuals_data_Tag += PseudoDataset_Tag->GetBinContent(bin+Beg);
	Residuals_data_Sub += PseudoDataset_Sub->GetBinContent(bin+Beg);
	model_Tag[bin]=Residuals_model_Tag;
	data_Tag[bin]=Residuals_data_Tag;
	model_Sub[bin]=Residuals_model_Sub;
	data_Sub[bin]=Residuals_data_Sub;

	phelectr[bin]=bin+Beg;
	y_Sub[bin]=(data_Sub[bin]-model_Sub[bin])/sqrt(data_Sub[bin]);
	y_Tag[bin]=(data_Tag[bin]-model_Tag[bin])/sqrt(data_Tag[bin]);

	}

	// Plotting PDFs in the main spectrum fit plot, after having been scaled properly

	int *Colors = new int [N]{632,632,842,799,600,921,870,800,632,844,425,616,616};

	for(int i=0; i<N; i++){
		for(int j=Beg; j<End; j++){
			hPDF[i]->SetBinError(j,0);
			}

		hPDF[i]->SetLineColor(Colors[i]);
		hPDF[i]->SetLineWidth(3);
		hPDF[i]->SetMarkerColor(Colors[i]);
		}

	for(int i=0; i<N; i++){
		hPDF[i]->Scale(Results[i]*Exposure_LT_Sub/Exposure_LT);

		for(int j=Beg; j<End;j++){
			if(i!=9 && i!=10 && i!= 12) hPDF_Total_Sub->SetBinContent(j,hPDF_Total_Sub->GetBinContent(j)+hPDF[i]->GetBinContent(j));
		}

       		if(i!=9 && i!=10 && i!= 12) hPDF[i]->Draw("HIST Lsame");			// Legend of Tag spectrum: not include C10_2, He_6, C11_2

	}

	hPDF_Total_Sub->SetLineColor(kRed);
	hPDF_Total_Sub->SetLineWidth(1);
	PseudoDataset_Sub->SetLineColor(kBlack);
	PseudoDataset_Sub->GetYaxis()->SetRangeUser(1,1e6);
	PseudoDataset_Sub->GetXaxis()->SetRangeUser(Beg,End);
	PseudoDataset_Sub->GetXaxis()->SetTitle("Reconstructed energy [p.e.]");
	PseudoDataset_Sub->GetYaxis()->SetTitle("Events");
	gPad->SetLogy();

	// LEGEND
	TLegend *leg1 = new TLegend(0.34,0.55,0.54,0.85,NULL,"brNDC");
	leg1->SetTextAlign(13);
	leg1->SetTextSize(0.04);
	leg1->SetBorderSize(0);
	leg1->SetFillStyle(0);
	for(int i=0; i<N;i++){
		if(i!=9 && i!=10 && i!= 12) leg1->AddEntry(hPDF[i],vinit[i],"l");		// Legend of Tag spectrum: only C10_2, He_6, C11_2
		}
	leg1->Draw("same");

	gStyle->SetOptStat(kFALSE);
	hPDF_Total_Sub->SetLineColor(kRed);
	hPDF_Total_Sub->SetLineWidth(1);
	hPDF_Total_Sub->Draw("same");

	c->cd();
	TPad*	 Pad_UR = new TPad("Pad_UR", "Pad_UR", 0.5, yLimPlot, 1.0, 1.0);
	Pad_UR->Draw();
	Pad_UR->cd();

	gPad->SetLogy();
	gStyle->SetOptTitle(0);
	gStyle->SetOptStat(0);

	PseudoDataset_Tag->SetLineColor(kBlack);
	PseudoDataset_Tag->GetYaxis()->SetRangeUser(1,1e6);
	PseudoDataset_Tag->GetXaxis()->SetRangeUser(Beg,End);
	PseudoDataset_Tag->Draw();
	PseudoDataset_Tag->GetXaxis()->SetTitle("Reconstructed energy [p.e.]");
	PseudoDataset_Tag->GetYaxis()->SetTitle("Events");

	for(int i=0; i<N; i++){
		  hPDF_Tag[i] = (TH1D*)hPDF[i]->Clone();
		  hPDF_Tag[i]->Scale(Exposure_LT_Tag/Exposure_LT_Sub);
		  hPDF_Tag[i]->SetLineColor(Colors[i]);
		  hPDF_Tag[i]->SetLineWidth(3);
		  hPDF_Tag[i]->SetMarkerColor(Colors[i]);
		  if(i!=11) hPDF_Tag[i]->Draw("HIST Lsame");
		}

	for(int i=0; i<N; i++){
		  for(int j=Beg; j<End;j++){
			     if(i!=11) hPDF_Total_Tag->SetBinContent(j,hPDF_Total_Tag->GetBinContent(j)+hPDF_Tag[i]->GetBinContent(j));
			}
	}

	// LEGEND
	TLegend *leg2 = new TLegend(0.34,0.75,0.54,0.85,NULL,"brNDC");
	leg2->SetTextAlign(13);
	leg2->SetTextSize(0.04);
	leg2->SetBorderSize(0);
	leg2->SetFillStyle(0);

  	for(int i=0; i<N;i++){
		if(i==9 || i==10 || i==12) leg2->AddEntry(hPDF_Tag[i],vinit[i],"l");
		}

  	leg2->Draw("same");
	hPDF_Total_Tag->SetLineColor(kRed);
	hPDF_Total_Tag->SetLineWidth(1);
	hPDF_Total_Tag->Draw("same");

  	c->cd();
  	TPad*    Pad_DR = new TPad("Pad_DR", "Pad_DR", 0.5, 0.0, 1.0, yLimPlot);
  	Pad_DR->Draw();
  	Pad_DR->cd();

	c->cd(4);
	TGraph *Residuals_Tag = new TGraph(End-Beg,phelectr,y_Tag);
	Residuals_Tag->SetTitle("Residuals (Sub)");
	Residuals_Tag->GetXaxis()->SetTitle("Reconstructed energy [p.e.]");
	Residuals_Tag->GetYaxis()->SetTitle("(D-M)/sqrt(D)");
	Residuals_Tag->GetYaxis()->CenterTitle(true);
	Residuals_Tag->GetYaxis()->SetTitleSize(.055);
	Residuals_Tag->GetXaxis()->SetTitleSize(.05);
	Residuals_Tag->GetXaxis()->SetRangeUser(Beg, End);
	Residuals_Tag->GetYaxis()->SetRangeUser(-4,4);
	Residuals_Tag->SetLineWidth(1);
	Residuals_Tag->Draw("AL");
	Residuals_Tag->Write("Residuals_Tag");

	c->cd();
	TPad*	 Pad_DL = new TPad("Pad_DL", "Pad_DL", 0.0, 0.0, 0.5, yLimPlot);
	Pad_DL->Draw();
	Pad_DL->cd();

	TGraph *Residuals_Sub = new TGraph(End-Beg,phelectr,y_Sub);
	Residuals_Sub->SetTitle("Residuals [Tag]");
	Residuals_Sub->GetXaxis()->SetTitle("Reconstructed energy [p.e.]");
	Residuals_Sub->GetYaxis()->SetTitle("(D-M)/sqrt(D)");
	Residuals_Sub->GetYaxis()->CenterTitle(true);
	Residuals_Sub->GetYaxis()->SetTitleSize(.055);
	Residuals_Sub->GetXaxis()->SetTitleSize(.05);
	Residuals_Sub->GetXaxis()->SetRangeUser(Beg, End);
	Residuals_Sub->GetYaxis()->SetRangeUser(-4,4);
	Residuals_Sub->SetLineWidth(1);
	Residuals_Sub->Draw("AL");
	Residuals_Sub->Write("Residuals_Sub");

	c->Write("c");

	TDirectory *PDFs = foutput->mkdir("PDFs");			  // Create a PDFs directory in Output_Rootfile
	PDFs->cd();

	for(int i=0; i<N; i++){hPDF[i]->Write(vinit[i]);}                 // Write PDFs in Output_Rootfile

	foutput->Close();

	TaggingPower = GetFloatPrecision(TaggingPower,2);                  // Prepare TP for Output_Text name
	Efficiency = GetFloatPrecision(Efficiency,2);                      // Prepare EFF for Output_Text name

  	Output_Text.erase(Output_Text.length()-4);
	Output_Text.append("_");
	Output_Text.append(to_string(int(TaggingPower*100)));
	Output_Text.append("_");
	Output_Text.append(to_string(int(Efficiency*100)));
	Output_Text.append(".txt");

	ofstream WriteOutputText;
	WriteOutputText.open(Output_Text.c_str(),ios::app);		   // APPENDING text output to Output_Text text file

	cout << endl;
	cout << "Exp = " << Exposure_LT << " days*100t , ExpSub = " << Exposure_LT_Sub << " days*100t , ExpTag = " << Exposure_LT_Tag << " days*100t " << endl;
	cout << "TaggingPower " << TaggingPower << "	 Efficiency " << Efficiency << endl;
	cout << "Rates expressed in [cpd/100t]" << endl;
	cout << setw(15) << "Species" << setw(15) << "InjRate" << setw(15) << "RecRate" << setw(15) << "ErrRate"  << setw(15) << "Bias[%]" << endl;
	
	double *RecRates = new double[N];
	double *ErrorRecRates = new double[N];
	double *BiasRecRates = new double[N];
	
	for(int i=0; i<N; i++){
	
		RecRates[i] = Results[i]*Threshold[i]/Exposure_LT;
		ErrorRecRates[i] = ErrorResults[i]*Threshold[i]/Exposure_LT;
		BiasRecRates[i] = 100.*(RecRates[i]/InjRate[i]-1);
	
		cout << setw(15) << vinit[i] << setprecision(4) <<  setw(15) << InjRate[i] << setw(15) << RecRates[i] << setw(15) << ErrorRecRates[i]  << setw(15) << BiasRecRates[i] << endl;
		WriteOutputText << setprecision(4) << RecRates[i] << '\t';
		WriteOutputText << setprecision(4) << ErrorRecRates[i]	 << '\t';
	}

	WriteOutputText << amin << endl;

	delete[] Threshold;
	delete[] Fix;
	delete[] Events;
	delete[] Events_Tag;
	delete[] Events_Sub;
	delete[] InjRate;
	delete[] InjRate_Tag;
	delete[] InjRate_Sub;
	delete[] Sigma;
	delete[] hPDF;
	delete[] hPDF_Tag;
	delete[] RecRates;
	delete[] ErrorRecRates;
	delete[] BiasRecRates;	

	return 0;

}

int main(int argc, char** argv) {
        string macro = argv[0];

        if(argc!=4) {
                cout << "\n     USAGE:  "<< macro << " Configuration_Text  Output_Rootfile  Output_Text \n" << endl;
                return 1;
        }

        string Configuration_Text = argv[1];
        string Output_Rootfile = argv[2];
	      string Output_Text = argv[3];

        MUST(Configuration_Text, Output_Rootfile, Output_Text);
        return 0;
}
