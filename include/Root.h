#pragma once
struct OutputData{
    int proton_num;
    int piminus_num;
    double Proton[20][5];//0~2:mom,3:mom.mag,4:mass
    double Piminus[20][5];
    double Lambda[5];
    double Kpp[5];
    double K0[5];
};

class ORoot{
    public:
        int OpenFile();
        int FillToRoot();
        int CloseFile();
        void set_proton_num(int value){output_data_.proton_num =  value;}
        void set_piminus_num(int value){output_data_.piminus_num = value;}
        void set_proton(int i, int j, double value){output_data_.Proton[i][j] =  value;}
        void set_piminus(int i, int j, double value){output_data_.Piminus[i][j] = value;}
        void set_lambda(int i, double value){output_data_.Lambda[i] = value;}
        void set_kpp(int i, double value){output_data_.Kpp[i] = value;}
        void set_k0(int i, double value){output_data_.K0[i] = value;}
    private:
        OutputData output_data_;
        TFile *ofile;
        TTree *tree;
};
