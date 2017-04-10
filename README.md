# DiPhoton

#### MG hepmc to root, run Delphe ####

* ./DelphesHepMC cards/converter_card.tcl output_MG.root input_MG.hepmc 
* ./DelphesHepMC cards/delphes_card_ATLAS.tcl output_Delphe.root input_MG.hepmc

#### Run the MG codes ####
    
* cd ggToaa/ggToaaMG/
* root -l ggToaaMG.C
* ggToaaMG t;
* t.Loop_MG();
 
#### Run the Delphe codes ####
    
* cd ggToaa/ggToaaDelphe/
* root -l ggToaaDelphe.C
* ggToaaDelphe t;
* t.Loop_Delphe();

#### Make Class of the root file ####

* TFile *f = new TFile("file.root")
* TTree *t = f->Get("Delphes")
* t->MakeClass("myclass")



