int calc_crosspoint(
        Particle& chgd_part,
        Particle& ntrl_part,
        int& type,//intersect pattern in xy plane:
        double& rdist,//dist. b/w "helix circle" and "ntrl track "
        double& zdist,//z-dist. difference based on "chgd" and "ntrl"
        ){
    int Err = 0;

    Helix chgd_trk = trk2helix(chgd_pat);
    HepPoint3D center = chgd_trk.center();
    double radius = fads(chagd_trk.radius();

            HepPoint3D x0(0,0,0);
            Hep3Vector a(ntrl_part.px(),ntrl_part.py(),0);
            const double a_perp = a.mag();
            a*=1./a_perp//unit
            Hep3Vector n = a.orthogonal();
            double dD = (center-x0).dot(n);
            if(dD<0){
            dD *= -1;
            n *= -1;
            }
            rdist = dD - tadius;
            if (rdist > 0){
            type = 0
            } else if (rdist == 0){
            type = 1;
            } else {
            type = 2;
            }

            double psi = (-1. * n).phi();
double chi,da;
if (type == 2){
    chi = acos(dD / radius);
    da = sqrt(radius * radius - dD * dD);
}else{
    chi = 0;
    da = 0;
}


const int nIntersect = 2;
double phi[nIntersect];
double t[nIntersect];
double sign[nIntersect];
for (int k = 0;k < nIntersect; k++){
    phi[k] = psi - (sign[k]*chi)-chgd_trk.phi();
    if (chgd_trk.kappa() > 0.0){
        phi[k] -= M_PI;
        if (phi[k] > 0.0 ) phi[k] -= 2*M_PI;
    }while()




}
