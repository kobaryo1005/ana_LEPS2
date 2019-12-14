double funVoigt(double *x,double *par){
    return par[0]*TMath::Voigt(x[0]-par[1],par[2],par[3],4);
}
