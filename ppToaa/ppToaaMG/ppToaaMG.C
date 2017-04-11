#define ppToaaMG_cxx
#include "ppToaaMG.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <math.h>      
#include "Math/VectorUtil.h"

void ppToaaMG::Loop_MG()
{
    //////////////////////////////////////////////////////////// 
    //
    // Section - 1 : Evaluate the invariant mass of diphoton  //
    //
    ////////////////////////////////////////////////////////////
    
    static const int binN = 20;
    int xmin = 200;
    int xmax = 1000;
    std::string histoTitle("MG: p p -> a a [QCD], 200K");
    //TFile* f = new TFile("Mgg_histo.root","recreate");
    TH1F* h = new TH1F("Histo", "diphoton histos", binN, xmin,xmax);
	
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
        vector<double>diPhoton_PT;
        vector<double>diPhoton_Phi;
        vector<double>diPhoton_Eta;
        vector<double>diPhoton_E;
        
        double Ht = 0.0;
        fChain->GetEntry(nentry);
        for(int i=0; i<Particle_size; i++){
          if(Particle_PID[i]==22 && Particle_Status[i]==1){
 	        diPhoton_PT.push_back(Particle_PT[i]);
 	        diPhoton_Eta.push_back(Particle_Eta[i]);
 	        diPhoton_Phi.push_back(Particle_Phi[i]);
 	        diPhoton_E.push_back(Particle_E[i]);
          }
          if(Particle_PID[i]!=22 && Particle_Status[i]==1){
            Ht = Ht+ Particle_PT[i];
          }
        }
        
        //cout<<"Ht = "<<Ht<<endl;  
        Float_t Mgg;
        Float_t prod_diPhoton_PT = diPhoton_PT[0]*diPhoton_PT[1];
        Float_t diff_diPhoton_Eta = cosh(diPhoton_Eta[0] -diPhoton_Eta[1]);
        Float_t diff_diPhoton_Phi = cos(diPhoton_Phi[0]-diPhoton_Phi[1]);
        
        //https://en.wikipedia.org/wiki/Invariant_mass
        Mgg = sqrt(abs(2*(prod_diPhoton_PT*(diff_diPhoton_Eta - diff_diPhoton_Phi))));
        //cout<<"Mg = ="<<Mgg<<endl;
        
        //Apply ATLAS and Ht cuts
        if(diPhoton_PT[0] >= 0.4*Mgg && diPhoton_PT[1] >= 0.3*Mgg && Mgg >= 200){    
        //if(diPhoton_PT[0] >= 0.4*Mgg && diPhoton_PT[1] >= 0.3*Mgg && Mgg >= 200 && Ht >= 50){    
    	
            h->Fill(Mgg); // Create a histogram of Mgg
        }

        diPhoton_PT.clear();
        diPhoton_Eta.clear();
        diPhoton_Phi.clear();
        diPhoton_E.clear();
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
    c_1->SetLogy();
    c_1->Divide(1,2);
    c_1->cd(1);
	// create the TGraphErrors and draw it
	plotBinErrors = new TGraphErrors(binN,binCenter,binContent,
						binCenterErr,binContentErr);
    plotBinErrors->SetTitle("MG: p p -> a a [QCD], 200K");
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
    //plotDiffErrors->SetTitle("MG: p p -> a a [QCD], 200K");

	TF1 *baseLine = new TF1("baseLine","0",0,2000); 
	baseLine->Draw("SAME");

    //Print the final fit/diff to the disk
	//c_1->Print("Mgg_diff.jpg");

    h->Write();
   // f->Close();

}
