#include <iostream>
#include <fstream>
#include <sstream>
#include <stdio.h>
#include <stdlib.h>
#include <TFile.h>
#include <TTree.h>
#include <TStyle.h>
#include <TLine.h>
#include <TLorentzVector.h>
#include <TROOT.h>
#include <TF1.h>
#include <TMath.h>
#include <TH1.h>
#include "Root.h"
// #include "Vertex.h"
// #include "PID.h"
// #include "FitFunction.h"

const int maxtrack = 54*4;
const int maxvertex = 500;
const int maxPrimary = 20;
float pgGen[4];  
int nPrimary;  
int geantID[maxPrimary];  
float gpx[maxPrimary];  
float gpy[maxPrimary];  
float gpz[maxPrimary];  
double dEdx[maxtrack];
double px[maxtrack];  
double py[maxtrack];  
double pz[maxtrack];  
double pp[maxtrack];  
double prbchi2[maxtrack];  
int charge[maxtrack];
double distance[maxvertex];
int itrk1[maxvertex];
int itrk2[maxvertex];
double x[maxvertex];
double y[maxvertex];
double z[maxvertex];
int n;

int main(int argc, char* argv[]){
    std::string infile = argv[1];
    std:: cout << "Reading ..." << infile << std::endl;
    TH1F *h1 = new TH1F("h1","Lambda_NoSelect",500,1.,1.5);
    TH1F *h2 = new TH1F("h2","Lambda",500,1.,1.5);
    TH1F *h3 = new TH1F("h3","Kpp",500,2.,2.5);
    TH1F *h4 = new TH1F("h4","K0",100,0.1,0.8);
    TFile *_file0 = TFile::Open(infile.c_str());
    TTree *gene = (TTree*)_file0->Get("gene");
    TTree *trk = (TTree*)_file0->Get("trk");
    TTree *vtx = (TTree*)_file0->Get("vtx");
    gene->SetBranchAddress("nPrimary",&nPrimary);
    gene->SetBranchAddress("pgGen",pgGen);
    gene->SetBranchAddress("gpx",gpx);
    gene->SetBranchAddress("gpy",gpy);
    gene->SetBranchAddress("gpz",gpz);
    gene->SetBranchAddress("geantID",geantID);
    trk->SetBranchAddress("n",&n);
    trk->SetBranchAddress("dEdx",dEdx);
    trk->SetBranchAddress("px",px);
    trk->SetBranchAddress("py",py);
    trk->SetBranchAddress("pz",pz);
    trk->SetBranchAddress("pp",pp);
    trk->SetBranchAddress("prbchi2",prbchi2);
    trk->SetBranchAddress("charge",charge);
    vtx->SetBranchAddress("distance",distance);
    vtx->SetBranchAddress("x",x);
    vtx->SetBranchAddress("y",y);
    vtx->SetBranchAddress("z",z);
    vtx->SetBranchAddress("itrk1",itrk1);
    vtx->SetBranchAddress("itrk2",itrk2);
    TLorentzVector proton[20];
    TLorentzVector piminus[20];
    double deuteron_mass = 1.875613;
    double proton_mass = 0.9382;
    double pi_mass = 0.1395;
    double z0 = -500.;
    int runnum = 0;
    int nvertex;
    int Allev  = trk ->GetEntries();
    int piminus_number;
    int proton_number;
    int proton_id[20];
    int piminus_id[20];
    int lambda_number=0;
    ORoot outfile;
    outfile.OpenFile();
    for(int i=0;i<Allev;i++){
        // for(int i=0;i<1000;i++){
        if(i%1000==0)printf("Event number %d\n",i);
        trk->GetEntry(i);
        gene->GetEntry(i);
        vtx->GetEntry(i);
        proton_number=0;
        piminus_number=0;
        nvertex = n * (n-1);
        TLorentzVector target;
        TLorentzVector gamma;
        TLorentzVector system;
        target.SetXYZM(0.,0.,0.,deuteron_mass);
        gamma.SetPxPyPzE(0.,0.,pgGen[2],pgGen[3]);
        system = target+gamma;
        // if(n>5) printf("many tracks!!\n");
        for(int j = 0;j<n;j++){
            if(prbchi2[j]<0.02)continue;
            if(charge[j]>0&&(dEdx[j]*10000.+70.*pp[j]-24.>0.&&dEdx[j]*10000+3.*pp[j]-5.>0)){
                proton[proton_number].SetXYZM(px[j],py[j],pz[j],proton_mass);
                outfile.set_proton( j,  0, px[j]);
                outfile.set_proton( j,  1, py[j]);
                outfile.set_proton( j,  2, pz[j]);
                double pp = sqrt(px[j]*px[j]+py[j]*py[j]+pz[j]*pz[j]);
                outfile.set_proton( j,  3, pp);
                proton_number++;}//proton select 
            if(charge[j]<0){
                piminus[piminus_number].SetXYZM(px[j],py[j],pz[j],pi_mass);
                outfile.set_piminus( j, 0, px[j]);
                outfile.set_piminus( j, 1, py[j]);
                outfile.set_piminus( j, 2, pz[j]);
                double pp = sqrt(px[j]*px[j]+py[j]*py[j]+pz[j]*pz[j]);
                outfile.set_piminus( j, 3, pp);
                piminus_number++;}//piminus select
        }
        outfile.set_proton_num(proton_number);
        outfile.set_piminus_num(piminus_number);
        for(int j = 0 ;j<proton_number;j++){
            for(int k = 0 ;k<piminus_number;k++){
                TLorentzVector lambda =proton[j]+piminus[k];
                h1->Fill(lambda.M());//no select
                for(int l=0;l<nvertex;l++){
                    if((itrk1[l]==proton_id[j]&&itrk2[l]==piminus_id[k])||(itrk2[l]==proton_id[j]&&itrk1[l]==piminus_id[k])){
                        if(distance[l]<10.&& sqrt(x[l]*x[l]+y[l]*y[l])>10.){
                            h2->Fill(lambda.M());//select lambda
                            outfile.set_lambda(0,lambda.Px());
                            outfile.set_lambda(1,lambda.Py());
                            outfile.set_lambda(2,lambda.Pz());
                            outfile.set_lambda(4,lambda.M());
                            lambda_number++;
                            if(proton_number>1&&lambda.M()>1.10&&lambda.M()<1.13){
                                for(int m=0;m<proton_number;m++){
                                    if(m!=j){
                                        TLorentzVector kpp = lambda + proton[m];
                                        h3->Fill(kpp.M());
                                        outfile.set_kpp(0,kpp.Px());
                                        outfile.set_kpp(1,kpp.Py());
                                        outfile.set_kpp(2,kpp.Pz());
                                        outfile.set_kpp(4,kpp.M());
                                        TLorentzVector K0 = system - kpp;
                                        h4->Fill(K0.M());
                                        // printf("K0.M = %lf\n",K0.M());
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }

        outfile.FillToRoot();
    }
    // WriteHist();
    h1->Write();
    h2->Write();
    h3->Write();
    h4->Write();

    outfile.CloseFile();
    return 0;
    }


