#define ppggToaaDelphe_cxx
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <math.h>      

#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>

void ppggToaaDelphe()
{
  //////////////////////////////////////////////////////////// 
  //                                                        //
  // Section - 1 : create histograms for pp and gg process  //
  //                                                        // 
  ////////////////////////////////////////////////////////////
/* 

  const char *title = "Delphe, gg>aa[200K]+ pp>aa[8000K], cut-1: Only ATLAS cuts";
  ifstream ggToaaAllCut("Mgg_cut1_ATLAS_ggToaaDelphe_200K_PY8.dat");
  ifstream ppToaaAllCut("Mgg_cut1_ATLAS_ppToaaDelphe_8000K_PY8.dat");
*/
  
 /* 
  const char *title = "Delphe, gg>aa[200K]+ pp>aa[8000K], cut-2: ATLAS cuts+ Events with a g-leading-jet";
  ifstream ggToaaAllCut("Mgg_cut2_gjet_ggToaaDelphe_200K_PY8.dat");
  ifstream ppToaaAllCut("Mgg_cut2_gjet_ppToaaDelphe_8000K_PY8.dat");
  
  const char *title = "Delphe, gg>aa[200K]+ pp>aa[8000K], cut-3: ATLAS cuts+ ( |eta_1| <= 1 && |eta_2| <= 1 ) ";
  ifstream ggToaaAllCut("Mgg_cut3_eta_ggToaaDelphe_200K_PY8.dat");
  ifstream ppToaaAllCut("Mgg_cut3_eta_ppToaaDelphe_8000K_PY8.dat");
  
  const char *title = "Delphe, gg>aa[200K]+ pp>aa[8000K], cut-4: ATLAS cuts+ C < 1 ";
  ifstream ggToaaAllCut("Mgg_cut4_C1_ggToaaDelphe_200K_PY8.dat");
  ifstream ppToaaAllCut("Mgg_cut4_C1_ppToaaDelphe_8000K_PY8.dat");
  
  const char *title = "Delphe, gg>aa[200K]+ pp>aa[8000K], cut-5: ATLAS cuts+ delta_R(photons) > 3.2 ";
  ifstream ggToaaAllCut("Mgg_cut5_dR3p2_ggToaaDelphe_200K_PY8.dat");
  ifstream ppToaaAllCut("Mgg_cut5_dR3p2_ppToaaDelphe_8000K_PY8.dat");
  
  const char *title = "Delphe, gg>aa[200K]+ pp>aa[8000K], cut-6: ATLAS cuts+ Events with a g-leading-jet + ( |eta_1| <= 1 && |eta_2| <= 1 )+  C < 1 + delta_R(photons) > 3.2 ";
  ifstream ggToaaAllCut("Mgg_cut6_ggToaaDelphe_200K_PY8.dat");
  ifstream ppToaaAllCut("Mgg_cut6_ppToaaDelphe_8000K_PY8.dat");
  
  const char *title = "Delphe, gg>aa[200K]+ pp>aa[8000K], cut-7: ATLAS cuts+ Events with a g-leading-jet + C < 1 + delta_R(photons) > 3.2";
  ifstream ggToaaAllCut("Mgg_cut7_ggToaaDelphe_200K_PY8.dat");
  ifstream ppToaaAllCut("Mgg_cut7_ppToaaDelphe_8000K_PY8.dat");
  */
  
  const char *title = "Delphe, gg>aa[200K]+ pp>aa[8000K], cut-8: ATLAS cuts+ leading and subleading g jets";
  ifstream ggToaaAllCut("Mgg_cut8_2gjet_ggToaaDelphe_200K_PY8.dat");
  ifstream ppToaaAllCut("Mgg_cut8_2gjet_ppToaaDelphe_8000K_PY8.dat");
  
  // Create histo of Mgg, for gg
  const int binN = 20;
  int xmin = 200;
  int xmax = 1000;
  TH1F* hgg = new TH1F("hgg", "Delphe: g g -> a a [QCD], 100K", binN, xmin,xmax);
  int NggToaaAllCut = 0;
  float ggMgg;
  while(true){
    ggToaaAllCut >> ggMgg; 
    //if(ggMgg<=280 || ggMgg>=400){
    hgg->Fill(ggMgg);
    //}
    //cout<<ggPt1<<"\t"<<ggMgg<<"\t"<<ggPt2<<endl;
    if( ggToaaAllCut.eof() ) break;
    NggToaaAllCut ++;
  }
  // Create histo of Mgg, for pp
  const int binN = 20;
  int xmin = 200;
  int xmax = 1000;
  TH1F* hpp = new TH1F("hpp", "Delphe: p p -> a a [QCD], 400K", binN, xmin,xmax);
  int NppToaaAllCut = 0;
  float ppMgg;
  while(true){
    ppToaaAllCut >>ppMgg; 
    ///if(ppMgg<=280 || ppMgg>=400){
      hpp->Fill(ppMgg);
    ///}
    //cout<<ppPt1<<"\t"<<ppMgg<<"\t"<<ppPt2<<endl;
    if( ppToaaAllCut.eof() ) break;
    NppToaaAllCut ++;
  }
  
  //Draw the histograms
  TCanvas* c1 = new TCanvas("c1","diphoton spectrum ");
  c1->Divide(2,1);
  c1->cd(1);
  hgg->Draw();
  //hgg->GetYaxis()->SetTitle("Events/40GeV");
  hgg->GetXaxis()->SetTitle("M_{gg} (GeV)");
  c1->cd(2);
  hpp->Draw();
  //hpp->GetYaxis()->SetTitle("Events/40GeV");
  hpp->GetXaxis()->SetTitle("M_{gg} (GeV)");
  
  //////////////////////////////////////////////////////////// 
  //                                                        //
  // Section - 2 : Combine events from the both histograms //
  //                                                        // 
  ////////////////////////////////////////////////////////////
  
  //Get the statistics from gg and pp
  int NppToaa = 8000000;
  float SppToaa = 1.53*102.96;
  int NppToaaPythia =  8000000;
  float SppToaaPythia = SppToaa*NppToaaPythia/NppToaa;
  int NggToaa = 200000;
  float SggToaa = 1.15*0.271;
  cout <<"Total gg events = "<<NggToaaAllCut<<endl;
  cout <<"Total pp events = "<<NppToaaAllCut<<endl;
  float SppToaaAllCut = SppToaaPythia * NppToaaAllCut/NppToaaPythia;
  float SggToaaAllCut = SggToaa * NggToaaAllCut/NggToaa;
  cout <<"Sigma of gg, after all cuts = "<<SggToaaAllCut<<endl;
  cout <<"Sigma of pp, after all cuts = "<<SppToaaAllCut<<endl;
  
  //Read the X-axis of the histogram
  Float_t binCenter[binN] = {};
  Float_t binContent[binN] = {};
  Float_t binCenterErr[binN] = {};
  Float_t binContentErr[binN] = {};
  
  TAxis * ppXaxis = hpp->GetXaxis();
  cout<<"gg"<<setw(20)<<"pp"<<setw(20)<<"pp+gg"<<endl;
  for(int j = 0; j < binN; j++)
  {
  	//Note, the bin starts form 1 
  	binCenter[j] = ppXaxis->GetBinCenter(j+1);
  	binCenterErr[j] = 0*sqrt(binCenter[j]);
    //binContent[j] = (float)hgg->GetBinContent(j+1)*((float)SggToaaAllCut/(float)SppToaaAllCut)*((float)NppToaaAllCut/(float)NggToaaAllCut);
    binContent[j] = (float)hpp->GetBinContent(j+1) + (float)hgg->GetBinContent(j+1)*((float)SggToaaAllCut/(float)SppToaaAllCut)*((float)NppToaaAllCut/(float)NggToaaAllCut);
    cout<<hgg->GetBinContent(j+1)<<setw(20)<<hpp->GetBinContent(j+1)<<setw(20)<<binContent[j]<<endl;
    binContentErr[j] = sqrt(binContent[j]); 
    binContentErr[j] = sqrt(binContent[j]); 
  }

  // create the TGraphErrors and draw it
  TCanvas* c_1 = new TCanvas("c_1","diphoton spectrum ");
  c_1->Divide(1,2);
  c_1->cd(1);
  c_1->SetLogy();
  plotBinErrors = new TGraphErrors(binN,binCenter,binContent,
  					binCenterErr,binContentErr);
 //plotBinErrors->SetTitle("Delphe: gg>aa[QCD] + pp>aa[QCD]");
  plotBinErrors->SetTitle(title);
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
  
  // ATLAS Fit 
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
  for(j = 0; j < binN; j++){
  binContentFit[j] = pow((1- pow(binCenter[j]/13000, 1/3.0)) ,b)*pow(binCenter[j]/13000, a0);
  //Difference between Data and the fitted background
  binContentDiff[j] =  binContent[j] - binContentFit[j];
  binContentDiffErr[j] = sqrt(binContent[j]);
  //binContentDiffErr[j] = sqrt(binContent[j] + binContentFit[j]);
  cout<<binCenter[j]<<setw(15)<<binContent[j]<<setw(20)<<
  binContentFit[j]<<setw(10)<<binContentDiff[j]<<
  setw(10)<<binContentDiffErr[j]<<endl;
  }	
  double chi2 = 0.0;
  for(j = 0; j < 5; j++){
    chi2 = chi2 + pow((binContent[j] - binContentFit[j]),2)/binContent[j];
  }
  cout<<"chi2 = "<<chi2<<endl;
/*
  // Gaussian Fit 
  Float_t s = 13000*13000;
  //fit with Gaussian --------- 
  TF1 *fitATLAS = new TF1("fit", "[0]*exp(-pow((x-200), 2)/pow([1], 2))");// - [2]*exp(-pow((x-300), 2)/pow([3], 2))" );
  //TF1 *fitATLAS = new TF1("fit", "exp(-pow((x-200)/[1], 2))");// - [2]*exp(-pow((x-300), 2)/pow([3], 2))" );
  cout<<endl;
  cout<<"================================"<<endl;
  plotBinErrors->Fit(fitATLAS);
  float a = fitATLAS->GetParameter(0);
  float b = fitATLAS->GetParameter(1);
  float c = fitATLAS->GetParameter(2);
  float d = fitATLAS->GetParameter(3);
  cout<< "a= "<<a<<setw(5)<<"b= "<<b<<setw(5)<<"c= "<<c<<setw(5)<<"d= "<<d<<endl;
  cout<<endl;
  cout<<"================================"<<endl;
  Float_t binContentFit[binN] = {};
  Float_t binContentDiff[binN] = {};
  Float_t binContentDiffErr[binN] = {};
  
  cout<<"binCenter"<<setw(15)<<"binContent"<<setw(10)<<"Fit"
  <<setw(10)<<"Diff"<<setw(15)<<"DiffError"<<endl;
  double chi2 = 0.0;
  //Evaluating fitted line at each binCenter	
  for(j = 0; j < binN; j++)
  {
  double x = binCenter[j]/13000.0;
  binContentFit[j] = a*exp(-pow((x-200)/b, 2) - c*exp(-pow((x-300)/d, 2));
  //binContentFit[j] = a*exp(-pow((x-200)/b, 2));
  //Difference between Data and the fitted background
  binContentDiff[j] =  binContent[j] - binContentFit[j];
  binContentDiffErr[j] = sqrt(binContent[j]);
  chi2 = chi2 + pow((binContent[j] - binContentFit[j]),2)/binContent[j];
  cout<<binCenter[j]<<setw(15)<<binContent[j]<<setw(20)<<
  binContentFit[j]<<setw(10)<<binContentDiff[j]<<
  setw(10)<<binContentDiffErr[j]<<endl;
  }	
  cout<<"chi2 = "<<chi2<<endl;
 */

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

  //h->Write();
 // f->Close();
}
