#define ggToaaDelphe_cxx
#include "ggToaaDelphe.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <math.h>      

#include <iostream>
#include <fstream>

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
  TH1F* h = new TH1F("Histo"," g g -> a a [200K], ATLAS cuts ", binN, xmin,xmax);
  
  Long64_t nentries = fChain->GetEntriesFast();
  
  ofstream fileMgg;
  //New cuts with 200K(gg) and 8000K(pp) events 
  //Output data files of Mgg
  //fileMgg.open("Mgg_cut1_ATLAS_ggToaaDelphe_200K_PY8.dat");
  //fileMgg.open("Mgg_cut2_gjet_ggToaaDelphe_200K_PY8.dat");
  //fileMgg.open("Mgg_cut3_eta_ggToaaDelphe_200K_PY8.dat");
  //fileMgg.open("Mgg_cut4_C1_ggToaaDelphe_200K_PY8.dat");
  //fileMgg.open("Mgg_cut5_dR3p2_ggToaaDelphe_200K_PY8.dat");
  //fileMgg.open("Mgg_cut6_ggToaaDelphe_200K_PY8.dat");
  //fileMgg.open("Mgg_cut7_ggToaaDelphe_200K_PY8.dat");
  fileMgg.open("Mgg_cut8_bothGjet_ggToaaDelphe_200K_PY8.dat");
 
  //Histograms
  TH1F* ptOf1stJet = new TH1F("gg: ptOf1stJet", "ptOf1stJet", 100, 0, 500);
  TH1F* ptOf1stQJet = new TH1F("gg: ptOf1stQJet", "ptOf1stQJet", 100, 0, 500);
  TH1F* etaPhot1 = new TH1F("gg: etaphot1", "eta of first photon", 100, -3, +3);
  TH1F* etaPhot2 = new TH1F("gg: etaphot2", "eta of second photon", 100, -3, +3);
  TH1F* histC = new TH1F("gg: cutC", "cutC", 100, 0, 10);
  TH1F* htOfjets = new TH1F("gg: htOfjets", "sum of transverse pt of jets", 100, 0, 500);
  TH1F* nJet = new TH1F("gg: nJet", "number of jets in an event", 100, 0, 15);
  TH1F* qJet = new TH1F("gg: qJet", "number of events that have quark leading-jet", 100, 0, 5);
  TH1F* gJet = new TH1F("gg: gJet", "number of events that have gluon leading-jet", 100, 0, 5);
  TH1F* jet0q1q = new TH1F("gg: jet0q1q", "n.o.e. having leading and subleading qjet", 100, 0, 5);
  TH1F* deltaR = new TH1F("gg: delta R", "angular separation between photons", 100, 0, 25);
  ///fileMgg<<"Mgg"<<setw(10)<<"C-cut"<<setw(10)<<"Jet_pt"<<setw(10)<<"Phot1_pt"<<setw(10)<<"Phot2_pt"<<endl;
  
  //Loop over number of events (nentries)   
  int count_2phot = 0;
  int counter = 0;
  for (Long64_t nentry=0; nentry<nentries; nentry++) 
  {
    if(nentry%100==0){
    cout<<"Event number = "<< nentry<<endl;
    }
    
    fChain->GetEntry(nentry);
    //fChain->GetEntry();
    if(Photon_size >=2){
      count_2phot =  count_2phot +1;
      Float_t Mgg; 
      Float_t prod_PT = Photon_PT[0]*Photon_PT[1];
      Float_t diff_Eta = cosh(Photon_Eta[0] -Photon_Eta[1]);
      Float_t diff_Phi = cos(Photon_Phi[0]-Photon_Phi[1]);
      //https://en.wikipedia.org/wiki/Invariant_mass
      Mgg = sqrt(abs(2*(prod_PT*(diff_Eta - diff_Phi))));
      
      //Ht of jets
      double Ht = 0;
      for(int j=0; j<Jet_size; j++){
        Ht = Ht+ Jet_PT[j];
      }
      //Jet multiplicity distributions
      int nJet_ = Jet_size;
      //delta R
      double dR = sqrt(pow(Photon_Eta[0]-Photon_Eta[1], 2) + pow(Photon_Phi[0]-Photon_Phi[1], 2));
      //the highest Pt jet must have Pt greater than 25
      //if(Jet_PT[0]>=50) continue;
      //C-cut
      Float_t Px1 = Photon_PT[0]* cos(Photon_Phi[0]);
      Float_t Py1 = Photon_PT[0]* sin(Photon_Phi[0]);
      Float_t Pz1 = Photon_PT[0]* sinh(Photon_Eta[0]);
      Float_t Px2 = Photon_PT[1]* cos(Photon_Phi[1]);
      Float_t Py2 = Photon_PT[1]* sin(Photon_Phi[1]);
      Float_t Pz2 = Photon_PT[1]* sinh(Photon_Eta[1]);
      Float_t cutC = sqrt(pow(Px1+Px2, 2)+pow(Py1+Py2, 2)+pow(Pz1+Pz2, 2));
      ///fileMgg<<Mgg<<setw(10)<<cutC<<setw(10)<<Jet_PT[0]<<setw(10)<<Photon_PT[0]<<setw(10)<<Photon_PT[1]<<endl;
      ///if(Jet_Flavor[0]==21 && Jet_PT[0]>=75) continue;
      
      //cut-1
      //if(Photon_PT[0] >= 0.4*Mgg && Photon_PT[1] >= 0.3*Mgg && Mgg >= 200){    
       ////cut-2
       //if(Jet_Flavor[0]==21){    
       ////cut-3
       //if(abs(Photon_Eta[0])<=1.0 && abs(Photon_Eta[1])<=1.0){    
       if(Photon_PT[0]>=50 && Photon_PT[1]>=50){    
       counter = counter +1;
       ////cut-4
       //if(cutC <= 1*Mgg){    
       ////cut-5
       //if(dR>3.2){    
       ////cut-6
       //if(Jet_Flavor[0]==21 && abs(Photon_Eta[0])<=1.0 && abs(Photon_Eta[1])<=1.0 && cutC <= 1*Mgg && dR>3.2){
       ////cut-7
       //if(Jet_Flavor[0]==21 && cutC <= 1*Mgg && dR>3.2){
       //if(Jet_Flavor[0]==21 && Jet_Flavor[1]==21){    
        
        h->Fill(Mgg); // Create a histogram of Mgg
        fileMgg<<Mgg<<"\n";
        ptOf1stJet->Fill(Jet_PT[0]);
        int qJet_ = 1;
        int gJet_ = 1;
        if(Jet_Flavor[0]!=21){
          qJet->Fill(qJet_);
          ptOf1stQJet->Fill(Jet_PT[0]);
        }
        else{
          gJet->Fill(gJet_);
        }
        int jet0q1q_ = 1;
        if(Jet_Flavor[0]!=21 && Jet_Flavor[1]!=21){
          jet0q1q->Fill(jet0q1q_);
        }
        etaPhot1->Fill(Photon_Eta[0]);
        etaPhot2->Fill(Photon_Eta[1]);
        htOfjets->Fill(Ht);
        histC->Fill(cutC/Mgg);
        nJet->Fill(nJet_);
        deltaR->Fill(dR);
        }
      //}
    }
  }
  //draw the Eta, Ht distributions
  TCanvas* c1 = new TCanvas("c1","Eta, Ht distributions");
  c1->Divide(3,2);
  c1->SetLogy();
  c1->cd(1);
  ptOf1stJet->Draw();
  c1->cd(2);
  ptOf1stQJet->Draw();
  c1->cd(3);
  etaPhot1->Draw();
  c1->cd(4);
  etaPhot2->Draw();
  c1->cd(5);
  htOfjets->Draw();
  c1->cd(6);
  histC->Draw();
  
  TCanvas* c2 = new TCanvas("c2","nJet and delta R, nJet and delta R");
  c2->Divide(2,1);
  c2->cd(1);
  nJet->Draw();
  c2->cd(2);
  deltaR->Draw();

  TCanvas* c3 = new TCanvas("c3","qJet vs gJet");
  c3->Divide(3,1);
  c3->cd(1);
  qJet->Draw();
  c3->cd(2);
  gJet->Draw();
  c3->cd (3);
  jet0q1q->Draw();
  fileMgg.close();
  cout << "============================="<<endl;
  cout << "  Total Events =  " << nentries <<endl;
  
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
  plotBinErrors->SetTitle("Delphe: g g -> a a [QCD], 200K");
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
  ///fitATLAS->Draw();
  cout<< "a0 = "<<a0<<","<<setw(5)<<"b ="<<b<<endl;
  cout<<endl;
  cout<<"================================"<<endl;

  Float_t binContentFit[binN] = {};
  Float_t binContentDiff[binN] = {};
  Float_t binContentDiffErr[binN] = {};
  
  cout<<"binCenter"<<setw(15)<<"binContent"<<setw(10)<<"Fit"
  <<setw(10)<<"Diff"<<setw(15)<<"DiffError"<<endl;

  //Evaluating fitted line at each binCenter	
  //Get chi2 in the Mgg range
  double chi2 = 0.0;
  for(j = 0; j < binN; j++)
  {
  binContentFit[j] = pow((1- pow(binCenter[j]/13000, 1/3.0)) ,b)
  										*pow(binCenter[j]/13000, a0);
  //Difference between Data and the fitted background
  binContentDiff[j] =  binContent[j] - binContentFit[j];
  //http://ipl.physics.harvard.edu/wp-uploads/2013/03/PS3_Error_Propagation_sp13.pdf
  binContentDiffErr[j] = sqrt(binContent[j]);
  //binContentDiffErr[j] = sqrt(binContent[j]+binContentFit[j]);
  chi2 = chi2 + pow((binContent[j] - binContentFit[j]),2)/binContentFit[j];
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
  //plotDiffErrors->SetTitle("Delphe: g g -> a a [QCD], 200K");

  TF1 *baseLine = new TF1("baseLine","0",0,2000); 
  baseLine->Draw("SAME");

  //Print the final fit/diff to the disk
  //c_1->Print("Mgg_diff.jpg");
  cout<<"number of entries = "<<nentries<<endl;
  cout<<"number of events having 2 or more photons ="<<count_2phot<<endl;
  cout<<"number of events after all cuts = "<<counter<<endl;
  h->Write();
 // f->Close();

}
