int PID(){
int SetPartileID(int *charge,double *dEdx,double *pp){
if(charge>0){
        if(dEdx*10000.+70.*pp-24.>0.&&dEdx*10000+3.*pp-5.>0.){
        particleID=2212;//proton
        } else{
        particleID=211;
        }
        }else{
            particleID=-211;
        }
}


}
