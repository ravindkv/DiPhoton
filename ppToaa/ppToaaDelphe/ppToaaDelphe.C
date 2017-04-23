#define ppToaaDelphe_cxx
#include "ppToaaDelphe.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <math.h>      

#include <iostream>
#include <fstream>

void ppToaaDelphe::Loop_Delphe()
{
  //////////////////////////////////////////////////////////// 
  //
  // Section - 1 : Evaluate the invariant mass of diphoton  //
  //
  ////////////////////////////////////////////////////////////
  
  const int binN = 20;
  int xmin = 200;
  int xmax = 1000;
  TH1F* h = new TH1F("Histo", "Delphe: p p -> a a [QCD], 4000K", binN, xmin,xmax);
  Long64_t nentries = fChain->GetEntriesFast();
  
  ofstream fileMgg;
  //C-cut: 2.5
  fileMgg.open("Mgg_C2p5_ppToaaDelphe_4000K_PY6Q.dat");
  //cut-0
  //fileMgg.open("Mgg_cut0_ppToaaDelphe_4000K_PY6Q.dat");
  //cut-1
  //fileMgg.open("Mgg_cut1_ppToaaDelphe_4000K_PY6Q.dat");
  //cut-2
  //fileMgg.open("Mgg_cut2_ppToaaDelphe_4000K_PY6Q.dat");
  //cut-3
  //fileMgg.open("Mgg_cut3_ppToaaDelphe_4000K_PY6Q.dat");
  //cut-4
  //fileMgg.open("Mgg_cut4_ppToaaDelphe_4000K_PY6Q.dat");

  TH1F* etaPhot1 = new TH1F("pp: etaphot1", "eta of first photon", 100, -3, +3);
  TH1F* etaPhot2 = new TH1F("pp: etaphot2", "eta of second photon", 100, -3, +3);
  TH1F* htOfjets = new TH1F("pp: htOfjets", "sum of transverse pt of jets", 100, 0, 500);

  //Loop over number of events (nentries)   
  for (Long64_t nentry=0; nentry<nentries; nentry++) 
  {
    if(nentry%100==0){
    cout<<"Event number = "<< nentry<<endl;
    }

    fChain->GetEntry(nentry);
    //fChain->GetEntry();
    double Ht = 0;
    if(Photon_size >=2){
      Float_t Mgg; 
      Float_t prod_PT = Photon_PT[0]*Photon_PT[1];
      Float_t diff_Eta = cosh(Photon_Eta[0] -Photon_Eta[1]);
      Float_t diff_Phi = cos(Photon_Phi[0]-Photon_Phi[1]);
      //https://en.wikipedia.org/wiki/Invariant_mass
      Mgg = sqrt(abs(2*(prod_PT*(diff_Eta - diff_Phi))));
     
      //Ht-cut
      for(int j=0; j<Jet_size; j++){
        Ht = Ht+ Jet_PT[j];
        //cout<<"Ht = "<<Ht<<endl;
      }
      //C-cut
      Float_t Px1 = Photon_PT[0]* cos(Photon_Phi[0]);
      Float_t Py1 = Photon_PT[0]* sin(Photon_Phi[0]);
      Float_t Pz1 = Photon_PT[0]* sinh(Photon_Eta[0]);
      Float_t Px2 = Photon_PT[1]* cos(Photon_Phi[1]);
      Float_t Py2 = Photon_PT[1]* sin(Photon_Phi[1]);
      Float_t Pz2 = Photon_PT[1]* sinh(Photon_Eta[1]);
      Float_t cutC = sqrt(pow(Px1+Px2, 2)+pow(Py1+Py2, 2)+pow(Pz1+Pz2, 2));
     
        //Ht-Cut
      //if(Photon_PT[0] >= 0.4*Mgg && Photon_PT[1] >= 0.3*Mgg && Mgg >= 200){    
      if(cutC <= 2.5*Mgg && Mgg >= 200){    
        //cut-0
        //if(Photon_PT[0] >= 0.4*Mgg && Photon_PT[1] >= 0.3*Mgg && Mgg >= 200){    
        //cut-1
        //if(cutC <= 2.0*Mgg && Photon_Eta[0]<=0.75 && Photon_Eta[1]<=0.75 &&Ht<=200){    
        //cut-2
        //if(cutC <= 0.5*Mgg && Photon_Eta[0]<=0.75 && Photon_Eta[1]<=0.75 &&Ht<=200){    
        //cut-3
        //if(cutC <= 2.0*Mgg && Photon_Eta[0]<=0.75 && Photon_Eta[1]<=0.75 &&Ht<=100){    
        //cut-4
        //if(cutC <= 0.5*Mgg && Photon_Eta[0]<=0.75 && Photon_Eta[1]<=0.75 &&Ht<=100){    
        h->Fill(Mgg); // Create a histogram of Mgg
        fileMgg<<Mgg<<"\n";
        etaPhot1->Fill(Photon_Eta[0]);
        etaPhot2->Fill(Photon_Eta[1]);
        htOfjets->Fill(Ht);
        //}
      }
    }
  }
  fileMgg.close();
  //draw the Eta, Ht distributions
  TCanvas* c1 = new TCanvas("c1","Eta, Ht distributions");
  c1->Divide(3,1);
  c1->cd(1);
  etaPhot1->Draw();
  c1->cd(2);
  etaPhot2->Draw();
  c1->cd(3);
  htOfjets->Draw();
  cout << "============================="<<endl;
  cout << "  Total Events =  " << nentries <<endl;
  cout << "============================="<<endl;
      
  //////////////////////////////////////////////////////////// 
  //
  // Section - 2 : Read the binCenter & binContent of histo //
  //
  ////////////////////////////////////////////////////////////

  Float_t binCenter[binN] = {};
  Float_t binContent[binN] = {};
  Float_t binCenterErr[binN] = {};
  Float_t binContentErr[binN] = {};

  //Read the X-axis of the histogram
  TAxis *xaxis = h->GetXaxis();
  cout<<endl;
  cout<<"================================"<<endl;
  cout<<setw(10)<<"bin"<<setw(10)<<"Center"<<setw(10)<<"Content"<<endl;
  for(int j = 0; j < binN; j++)
  {
  	//Note, the bin starts form 1 
  	binCenter[j] = xaxis->GetBinCenter(j+1);
  	binCenterErr[j] = 0*sqrt(binCenter[j]); 
  	
  	binContent[j] = h->GetBinContent(j+1); 
  	binContentErr[j] = sqrt(binContent[j]); 
  	cout<<setw(10)<<j<<setw(10)<<binCenter[j]<<setw(10)<<binContent[j]<<endl;	
  }

  TCanvas* c_1 = new TCanvas("c_1","diphoton spectrum ");
  c_1->Divide(1,2);
  c_1->cd(1);
  c_1->SetLogy();
  // create the TGraphErrors and draw it
  plotBinErrors = new TGraphErrors(binN,binCenter,binContent,
  					binCenterErr,binContentErr);
  plotBinErrors->SetTitle("Delphe: p p -> a a [QCD], 1000K");
  plotBinErrors->GetYaxis()->SetTitle("Events/40GeV");
  plotBinErrors->SetMarkerColor(1);
  plotBinErrors->SetMarkerStyle(20);
  plotBinErrors->Draw("AP");
  plotBinErrors->GetXaxis()->SetTitle("M_{gg} (GeV)");

  
  //////////////////////////////////////////////////////////// 
  //
  // Section - 3 : Do the chi2 square fitting of the histo  //
  //
  ////////////////////////////////////////////////////////////
  
  
  Float_t s = 13000*13000;
  TF1 *fitATLAS = new TF1("fit", "pow(1-pow(x/13000,1/3),[1])*pow(x/13000,[0])");

  cout<<endl;
  cout<<"================================"<<endl;
  plotBinErrors->Fit(fitATLAS);
  
  //Fitted parameteres are a0 = p0, b = p1
  float a0 = fitATLAS->GetParameter(0);
  float b = fitATLAS->GetParameter(1);

  cout<< "a0 = "<<a0<<","<<setw(5)<<"b ="<<b<<endl;
  cout<<endl;
  cout<<"================================"<<endl;

  Float_t binContentFit[binN] = {};
  Float_t binContentDiff[binN] = {};
  Float_t binContentDiffErr[binN] = {};
  
  cout<<"binCenter"<<setw(15)<<"binContent"<<setw(10)<<"Fit"
  <<setw(10)<<"Diff"<<setw(15)<<"DiffError"<<endl;

  //Evaluating fitted line at each binCenter	
  for(j = 0; j < binN; j++)
  {
  binContentFit[j] = pow((1- pow(binCenter[j]/13000, 1/3.0)) ,b)
  										*pow(binCenter[j]/13000, a0);
  //Difference between Data and the fitted background
  binContentDiff[j] =  binContent[j] - binContentFit[j];
  binContentDiffErr[j] = sqrt(binContent[j]);

  cout<<binCenter[j]<<setw(15)<<binContent[j]<<setw(20)<<
  binContentFit[j]<<setw(10)<<binContentDiff[j]<<
  setw(10)<<binContentDiffErr[j]<<endl;
  }	

  
  //////////////////////////////////////////////////////////// 
  //
  // Section - 4 : Plotting Data - Fitted background        //
  //
  ////////////////////////////////////////////////////////////
  
  
  //TCanvas* c_2 = new TCanvas("c_2","histo demo ", 500, 600,500,800);
  c_1->cd(2);
  // create the TGraphErrors and draw it
  plotDiffErrors = new TGraphErrors(binN,binCenter,binContentDiff,
                      binCenterErr,binContentDiffErr);
  plotDiffErrors->SetMarkerColor(1);
  plotDiffErrors->SetMarkerStyle(20);
  plotDiffErrors->Draw("AP");
  plotDiffErrors->GetXaxis()->SetTitle("M_{gg} (GeV)");
  plotDiffErrors->GetYaxis()->SetTitle("Data - fitted background");
  //plotDiffErrors->SetTitle("g g -> a a [QCD], 100K");

  TF1 *baseLine = new TF1("baseLine","0",0,2000); 
  baseLine->Draw("SAME");

  //Print the final fit/diff to the disk
  //c_1->Print("Mgg_diff.jpg");

  h->Write();
 // f->Close();

}

