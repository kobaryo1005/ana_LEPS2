#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1.h>
#include <TString.h>
#include "Root.h"

int ORoot::OpenFile(){
    ofile = new TFile("../bin/output.root","recreate");
    tree = new TTree("tree","/tree");
    tree->Branch("proton_num",&output_data_.proton_num,"proton_num/I");
    tree->Branch("piminus_num",&output_data_.piminus_num,"piminus_num/I");
    tree->Branch("Lambda",output_data_.Lambda,"Lambda[5]/D");
    tree->Branch("Kpp",output_data_.Kpp,"Kpp[5]/D");
    tree->Branch("Proton",output_data_.Proton,"Proton[20][5]/D");
    tree->Branch("Piminus",output_data_.Piminus,"Piminus[20][5]/D");
    tree->Branch("K0",output_data_.K0,"K0[5]/D");
    return 0;
}

int ORoot::CloseFile(){
    ofile->Write();
    ofile->Close();
    return 0;
}

int ORoot::FillToRoot(){
    tree->Fill();
    return 0;
}
