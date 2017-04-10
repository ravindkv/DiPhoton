#define ggToaaDelphe_cxx
#include "ggToaaDelphe.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <math.h>      

void ggToaaDelphe::Loop_Delphe()
{
  //////////////////////////////////////////////////////////// 
  //
  // Section - 1 : Evaluate the invariant mass of diphoton  //
  //
  ////////////////////////////////////////////////////////////
  
  const int binN = 20;
  int xmin = 200;
  int xmax = 1000;
  //TFile* f = new TFile("Mgg_histo.root","recreate");
  TH1F* h = new TH1F("Histo"," g g -> a a [100K], ATLAS cuts ", binN, xmin,xmax);
  
  Long64_t nentries = fChain->GetEntriesFast();
  cout << "============================="<<endl;
  cout << "  Total Events =  " << nentries <<endl;
  cout << "============================="<<endl;
 
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
      
      for(int j=0; j<Jet_size; j++){
        Ht = Ht+ Jet_PT[j];
        //cout<<"Ht = "<<Ht<<endl;
      }
      //Apply ATLAS and Ht cuts
      if(Photon_PT[0] >= 0.4*Mgg && Photon_PT[1] >= 0.3*Mgg && Mgg >= 200){    
      //if(Photon_PT[0] >= 0.4*Mgg && Photon_PT[1] >= 0.3*Mgg && Mgg >= 200i && Ht >= 50){    
        h->Fill(Mgg); // Create a histogram of Mgg
      }
    }
  
  }
      
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
  plotBinErrors->SetTitle("Delphe: g g -> a a [QCD], 100K");
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
  //plotDiffErrors->SetTitle("Delphe: g g -> a a [QCD], 100K");

  TF1 *baseLine = new TF1("baseLine","0",0,2000); 
  baseLine->Draw("SAME");

  //Print the final fit/diff to the disk
  //c_1->Print("Mgg_diff.jpg");

  h->Write();
 // f->Close();

}
