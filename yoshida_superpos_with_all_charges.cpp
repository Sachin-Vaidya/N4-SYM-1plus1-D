#include<iostream>
#include<cmath>
#include<complex>
#include<cstdio>
#include<cstdlib>
#include<ctime>
#include <fstream>
#include <omp.h>
#include<random>
#include<chrono>
#include <sstream>
#include <vector>
//#include<sys/types.h>
//#include<sys/stat.h>
#include</home/vaidya2/build_dir/source_dir/Eigen/Dense>
//#include "home/vaidya2/build_dir/Eigen/Core"
#include <gmp.h>
#include<iomanip>
#include<limits.h>
#include <cfloat>

using namespace std;
typedef std::complex<double> dcomp;

int main()
{
    double T;
    double Tmin = 0.;
    double Tmax = 50.;
    double Zmin = 0.;
    double Zmax = 20.;
    double h = 0.004;
    double tau = 0.001;
    double gamma = tau/pow(h,2);
    int N = round((Zmax-Zmin)/h);
    //double tau = gamma*pow(h,2);
    int N_t = round((Tmax-Tmin)/tau);
    double pi = 3.141592653589793238;
    int T_k;
    double AMP = 0.00001;
    
    cout<<"N = "<<N<<endl;
    
    int in = round(Zmax-Zmin);
    cout<<"in = "<<in<<endl;
    
    double f[6][3][2*N];
    
    double hamiltonian_density[N], momentum_density[N];
    double charge_density_so6[15][N], charge_density_su2[3][N];
    double charge_so6[15], charge_su2[3];
    double rand_amp[3][3][in];
    
    double c_n[4]={0.6756,-0.1756,-0.1756,0.6756};
    double d_n[4]={1.3512,-1.7024,1.3512,0};
    
    dcomp Phi0[3][3][N], dPhi0[3][3][N], iPhiDot0[3][3][N], diPhiDot0[3][3][N], dPhi0_x[3][3][N];
    dcomp Phi[3][3][N], iPhiDot[3][3][N];
    dcomp Phi1[3][3][N], iPhiDot1[3][3][N];
    dcomp Del2Phi[3][3][N];
    dcomp NonlinearPhi[3][3][N];
    
    cout<<"h = "<<h<<endl;
    cout<<"tau = "<<tau<<endl;
    cout<<"N = "<<N<<endl;
    cout<<"N_t = "<<N_t<<endl;
    cout<<"gamma = "<<gamma<<endl;
    
    double cccRe = gamma/4;
    double uRe = 4.0/3;
    
    dcomp p;
    dcomp p2(1,0);
    dcomp p1(0,0);
    dcomp b(0,1);
    dcomp nonlin(2,0);
    dcomp c(-1,0);
    dcomp cc(1,0);
    dcomp ccc(cccRe,0);
    dcomp cccc(100,0);
    dcomp u(uRe,2);
    dcomp v(2,3);
    dcomp omega;
    
    double coefficientRe = 1./sqrt(2);
    dcomp coefficient(coefficientRe,0);
    
    cout<<"cccRe = "<<cccRe<<endl;
    cout<<"uRe = "<<uRe<<endl;
    cout<<"coefficientRe = "<<coefficientRe<<endl;
    
    dcomp p3(0,0);
    dcomp co(0.25,0);
    dcomp f0(1./sqrt(8),0.);
    
    cout<<setprecision(16)<<"f0 = "<<f0<<endl;
    double one_third=pow(3,-1);
    cout<<"one_third = "<<one_third<<endl;
    
    double Energy_total, Momentum_total;
    
    double Hamiltonian_Density_x_derivative_term[3][3][2*N];
    double Hamiltonian_Density_nonlin[N];
    double Hamiltonian_Density_lin_dot[N];
    double Hamiltonian_Density_lin_grad[N];
    
    double NormalModes[9][N];
    
    
    //--------------Setting up fluctuation part of the initial conditions for normal modes----------------//
    std::default_random_engine de(static_cast<unsigned int>(time(0)));
    std::normal_distribution<double> nd(7.,2.);
    for(int m=0;m<3;m++)
    {
        for(int n=0;n<3;n++)
        {
            for(int nn=1;nn<=in;nn++)
            {
                rand_amp[m][n][nn-1] = nd(de);
            }
        }
    }
    
    for(int m=0;m<3;m++)
    {
        for(int n=0;n<3;n++)
        {
            for(int j=0;j<N;j++)
            {
                dPhi0[m][n][j]=0;
                diPhiDot0[m][n][j]=0;
                for(int nn=1;nn<=in;nn++)
                {
                    dcomp phase(2*pi*h*j*nn/(in),0);
                    
                    dPhi0[m][n][j]+=(rand_amp[m][n][nn-1]*exp(b*phase));
                    //diPhiDot0[m][n][j]+=((nd(de)+(b*nd(de)))*exp(b*phase));
                    
                    dPhi0_x[m][n][j]+=(b*(2*pi*nn/(in))*nd(de)*exp(b*phase));
                }
            }
        }
    }
    
    omp_set_num_threads(10);
    
    //--------------Initial conditions for fields----------------//
    //#pragma omp parallel for num_threads(20)
    for(int m=0;m<3;m++)
    {
        for(int n=0;n<3;n++)
        {
            for(int j=0;j<N;j++)
            {
                if(m==n)
                {
                    Phi[m][n][j] = f0 + AMP*dPhi0[m][n][j];
                    iPhiDot[m][n][j] = f0 + AMP*diPhiDot0[m][n][j];
                }
                else
                {
                    Phi[m][n][j] = AMP*dPhi0[m][n][j];
                    iPhiDot[m][n][j] = AMP*diPhiDot0[m][n][j];
                }
            }
        }
    }
    
    //std::string energy_file("myfile_energy"),total_charge_so6("total_charge_so6_"),total_charge_su2("total_charge_su2_"),charge_initial_density_file_so6("charge_initial_density_file_so6_"),charge_initial_density_file_su2("charge_initial_density_file_su2_"),myfile_charge_initial_so6("myfile_charge_initial_so6_"),myfile_charge_initial_su2("myfile_charge_initial_su2_"),charge_mid_density_file_so6("charge_mid_density_file_so6_"),charge_mid_density_file_su2("charge_mid_density_file_su2_"),myfile_charge_mid_so6("myfile_charge_mid_so6_"),myfile_charge_mid_su2("myfile_charge_mid_su2_"),charge_final_density_file_so6("charge_final_density_file_so6_"),charge_final_density_file_su2("charge_final_density_file_su2_"),myfile_charge_final_so6("myfile_charge_final_so6_"),myfile_charge_final_su2("myfile_charge_final_su2_"),interaction_term_at_t("/interaction_term/interaction_term_at_t"),kinetic_dot_term_at_t("/kinetic_dot_term/kinetic_dot_term_at_t"),kinetic_grad_term_at_t("/kinetic_grad_term/kinetic_grad_term_at_t"),charge_den_so6_at_t("/charge_den_so6/charge_den_so6_at_t"),charge_den_su2_at_t("/charge_den_su2/charge_den_su2_at_t"),L_squ,fields_at_t("/fields/fields_at_t_"),fields_amp_at_t("fields_amp/fields_amp_at_t_"),ext(".txt");
    //std::string energy_file("myfile_energy"),total_charge_so6("total_charge_so6_"),total_charge_su2("total_charge_su2_"),interaction_term_at_t("interaction_term/interaction_term_at_t"),kinetic_dot_term_at_t("kinetic_dot_term/kinetic_dot_term_at_t"),kinetic_grad_term_at_t("kinetic_grad_term/kinetic_grad_term_at_t"),charge_den_so6_at_t("charge_den_so6/charge_den_so6_at_t"),charge_den_su2_at_t("charge_den_su2/charge_den_su2_at_t"),hamiltonian_density_at_t("hamiltonian_density/hamiltonian_density_at_t"),fields_phi_at_t("fields_phi/fields_phi_at_t_"),fields_chi_at_t("fields_chi/fields_chi_at_t_"),fields_phi_amp_at_t("fields_phi_amp/fields_phi_amp_at_t_"),ext(".txt");
    std::string energy_file("myfile_energy"),momentum_file("myfile_momentum"),total_charge_so6("total_charge_so6_"),total_charge_su2("total_charge_su2_"),interaction_term_at_t("interaction_term/interaction_term_at_t"),kinetic_dot_term_at_t("kinetic_dot_term/kinetic_dot_term_at_t"),kinetic_grad_term_at_t("kinetic_grad_term/kinetic_grad_term_at_t"),charge_den_so6_at_t("charge_den_so6/charge_den_so6_at_t"),charge_den_su2_at_t("charge_den_su2/charge_den_su2_at_t"),hamiltonian_density_at_t("hamiltonian_density/hamiltonian_density_at_t"),momentum_density_at_t("momentum_density/momentum_density_at_t"),fields_phi_at_t("fields_phi/fields_phi_at_t_"),fields_chi_at_t("fields_chi/fields_chi_at_t_"),fields_phi_amp_at_t("fields_phi_amp/fields_phi_amp_at_t_"),normal_mode_at_t("normal_mode/normal_mode_at_t"),ext(".txt");
    
    char num[100];
    char num2[1000];
    
    string list[20];
    
    std::ofstream *total_charge_so6_file = new ofstream[15];
    for(int i=0;i<5;i++)
    {
        for(int j=i+1;j<6;j++)
        {
            stringstream so6_total;
            sprintf(num,"%d",((10*(i+1))+j+1));
            so6_total<<(total_charge_so6+num+ext);
            list[(5*i)+j-1-((i+1)*i/2)] = so6_total.str();
        }
    }
    for(int i=0;i<5;i++)
    {
        for(int j=i+1;j<6;j++)
        {
            total_charge_so6_file[(5*i)+j-1-((i+1)*i/2)].open(list[(5*i)+j-1-((i+1)*i/2)].c_str());
        }
    }
    
    
    std::ofstream *total_charge_su2_file = new ofstream[3];
    for(int i=0;i<3;i++)
    {
        stringstream su2_total;
        sprintf(num,"%d",i+1);
        su2_total<<(total_charge_su2+num+ext);
        list[15+i] = su2_total.str();
    }
    for(int i=0;i<3;i++)
    {
        total_charge_su2_file[i].open(list[15+i].c_str());
    }
    
    std::ofstream myfile_energy;
    stringstream energy_total;
    energy_total<<(energy_file+ext);
    list[19] = energy_total.str();
    myfile_energy.open(list[19].c_str());
    
    std::ofstream myfile_momentum;
    stringstream momentum_total;
    momentum_total<<(momentum_file+ext);
    list[19] = momentum_total.str();
    myfile_momentum.open(list[19].c_str());
    
    //--------------Time loop for evolution of fields----------------//
    for(int k=1;k<2*N_t+2;k++)
    {
        for(int j=0;j<N;j++)
        {
            NormalModes[0][j]=std::norm(Phi[0][0][j]-Phi[1][1][j]);
            NormalModes[1][j]=std::norm(Phi[1][1][j]-Phi[2][2][j]);
            NormalModes[2][j]=std::norm(Phi[0][0][j]+Phi[1][1][j]+Phi[2][2][j]);
            NormalModes[3][j]=std::norm(Phi[0][1][j]+Phi[1][0][j]);
            NormalModes[4][j]=std::norm(Phi[1][2][j]+Phi[2][1][j]);
            NormalModes[5][j]=std::norm(Phi[0][2][j]+Phi[2][0][j]);
            NormalModes[6][j]=std::norm(Phi[0][1][j]-Phi[1][0][j]);
            NormalModes[7][j]=std::norm(Phi[1][2][j]-Phi[2][1][j]);
            NormalModes[8][j]=std::norm(Phi[0][2][j]-Phi[2][0][j]);
        }
        
        for(int m=0;m<3;m++)
        {
            for(int n=0;n<3;n++)
            {
                for(int i=0;i<N;i++)
                {
                    f[2*(m+1)-2][n][i]=Phi[m][n][i].real();
                    f[2*(m+1)-1][n][i]=Phi[m][n][i].imag();
                    f[2*(m+1)-2][n][N+i]=(-1.*b*iPhiDot[m][n][i]).real();
                    f[2*(m+1)-1][n][N+i]=(-1.*b*iPhiDot[m][n][i]).imag();
                }
            }
        }
    
    
        //--------------Computing Energy density and Energy----------------//
        for(int i=0;i<N;i++)
        {
            Hamiltonian_Density_nonlin[i]=0;
            for(int m=0;m<3;m++)
            {
                for(int n=0;n<3;n++)
                {
                    for(int o=0;o<3;o++)
                    {
                        for(int p=0;p<3;p++)
                        {
                            Hamiltonian_Density_nonlin[i] += (2*(pow((Phi[m][n][i].real()),2)+pow((Phi[m][n][i].imag()),2))*(pow((Phi[o][p][i].real()),2)+pow((Phi[o][p][i].imag()),2)))-(2*(((Phi[o][p][i].real())*(Phi[o][n][i].real()))+((Phi[o][p][i].imag())*(Phi[o][n][i].imag())))*(((Phi[m][p][i].real())*(Phi[m][n][i].real()))+((Phi[m][p][i].imag())*(Phi[m][n][i].imag()))));
                        }
                    }
                }
            }
        }
    
        #pragma omp parallel for num_threads(20)
        for(int i=0;i<N;i++)
        {
            if(i==0 || i==N-1)
            {
                for(int m=0;m<3;m++)
                {
                    for(int n=0;n<3;n++)
                    {
                        Hamiltonian_Density_x_derivative_term[m][n][0] = (Phi[m][n][1]-Phi[m][n][N-1]).real();
                        Hamiltonian_Density_x_derivative_term[m][n][N-1] = (Phi[m][n][0]-Phi[m][n][N-2]).real();
                        Hamiltonian_Density_x_derivative_term[m][n][N] = (Phi[m][n][1]-Phi[m][n][N-1]).imag();
                        Hamiltonian_Density_x_derivative_term[m][n][2*N-1] = (Phi[m][n][0]-Phi[m][n][N-2]).imag();
                    }
                }
            }
            else
            {
                for(int m=0;m<3;m++)
                {
                    for(int n=0;n<3;n++)
                    {
                        Hamiltonian_Density_x_derivative_term[m][n][i] = (Phi[m][n][i+1]-Phi[m][n][i-1]).real();
                        Hamiltonian_Density_x_derivative_term[m][n][i+N] = (Phi[m][n][i+1]-Phi[m][n][i-1]).imag();
                    }
                }
            }
        }
    
        for(int i=0;i<N;i++)
        {
            Hamiltonian_Density_lin_dot[i]=0;
            Hamiltonian_Density_lin_grad[i]=0;
            for(int m=0;m<3;m++)
            {
                for(int n=0;n<3;n++)
                {
                    Hamiltonian_Density_lin_dot[i] += std::norm(-1.*b*iPhiDot[m][n][i]);
                    Hamiltonian_Density_lin_grad[i] += ((pow(Hamiltonian_Density_x_derivative_term[m][n][i],2)+pow(Hamiltonian_Density_x_derivative_term[m][n][N+i],2))*pow(2.*h,-2.));
                }
            }
        }
    
    
        for(int i=0;i<N;i++)
        {
            hamiltonian_density[i]=Hamiltonian_Density_lin_dot[i]+Hamiltonian_Density_lin_grad[i]+Hamiltonian_Density_nonlin[i];
        }
    
        Energy_total=0;
        for(int i=0;i<N;i++)
        {
            Energy_total+=hamiltonian_density[i]*h;
        }
        cout.precision(16);
        //cout<<"Energy = "<<Energy_total*h<<endl;
        
        
        for(int i=0;i<N;i++)
        {
            momentum_density[i]=0.;
            for(int j=0;j<3;j++)
            {
                for(int l=0;l<3;l++)
                {
                    momentum_density[i] += (Hamiltonian_Density_x_derivative_term[j][l][i]*((-1.*b*iPhiDot[j][l][i]).real()) + Hamiltonian_Density_x_derivative_term[j][l][i+N]*((-1.*b*iPhiDot[j][l][i]).imag()))/(2.*h);
                }
            }
        }
        
        Momentum_total=0;
        for(int i=0;i<N;i++)
        {
            Momentum_total+=momentum_density[i]*h;
        }
        cout.precision(16);
    
    
    
        //--------------Computing Charges and densities---------------//
        for(int i=0;i<N;i++)
        {
            for(int m=0;m<5;m++)
            {
                for(int n=m+1;n<6;n++)
                {
                    charge_density_so6[(5*m)+n-1-((m+1)*m/2)][i]=0;
                    for(int j=0;j<3;j++)
                    {
                        charge_density_so6[(5*m)+n-1-((m+1)*m/2)][i]+=((f[m][j][N+i]*f[n][j][i])-(f[m][j][i]*f[n][j][N+i]));
                    }
                }
            }
        }
    
        for(int i=0;i<N;i++)
        {
            for(int m=0;m<2;m++)
            {
                for(int n=m+1;n<3;n++)
                {
                    charge_density_su2[m+n-1][i]=0;
                    for(int j=0;j<6;j++)
                    {
                        charge_density_su2[m+n-1][i]+=((f[j][m][N+i]*f[j][n][i])-(f[j][m][i]*f[j][n][N+i]));
                    }
                }
            }
        }
    
        for(int m=0;m<15;m++)
        {
            charge_so6[m]=0;
            for(int i=0;i<N;i++)
            {
                charge_so6[m]+=charge_density_so6[m][i]*h;
            }
        }
    
        for(int m=0;m<3;m++)
        {
            charge_su2[m]=0;
            for(int i=0;i<N;i++)
            {
                charge_su2[m]+=charge_density_su2[m][i]*h;
            }
        }
        
        //----------------------writing fields, normal modes, charge densities and hamiltonial densities at multiple time steps-------------------------//
        if(k % 2 == 1)
        {
            T_k=(k-1)*tau/2;
            sprintf(num2,"%d",T_k);
            //cout<<"T_k = "<<T_k<<endl;
            std::ofstream real_phi_fields((fields_phi_at_t+num2+ext).c_str());
            std::ofstream real_phi_fields_amp((fields_phi_amp_at_t+num2+ext).c_str());
            std::ofstream real_chi_fields((fields_chi_at_t+num2+ext).c_str());
            std::ofstream kinetic_dot_term((kinetic_dot_term_at_t+num2+ext).c_str());
            std::ofstream kinetic_grad_term((kinetic_grad_term_at_t+num2+ext).c_str());
            std::ofstream interaction_term((interaction_term_at_t+num2+ext).c_str());
            std::ofstream charge_den_so6((charge_den_so6_at_t+num2+ext).c_str());
            std::ofstream charge_den_su2((charge_den_su2_at_t+num2+ext).c_str());
            std::ofstream hamiltonian_density_file((hamiltonian_density_at_t+num2+ext).c_str());
            std::ofstream momentum_density_file((momentum_density_at_t+num2+ext).c_str());
            std::ofstream normal_mode_file((normal_mode_at_t+num2+ext).c_str());
            if(real_phi_fields.is_open() || real_phi_fields_amp.is_open() || real_chi_fields.is_open() || kinetic_dot_term.is_open() || kinetic_grad_term.is_open() || interaction_term.is_open() || charge_den_so6.is_open() || charge_den_su2.is_open() || hamiltonian_density_file.is_open() || normal_mode_file.is_open())
            {
                for(int count=0;count<N;count++)
                {
                    real_phi_fields<<(count+1)*h<<"\t";
                    real_phi_fields_amp<<(count+1)*h<<"\t";
                    real_chi_fields<<(count+1)*h<<"\t";
                    for(int m=0;m<3;m++)
                    {
                        for(int n=0;n<3;n++)
                        {
                            real_phi_fields<<std::setprecision(16)<<Phi[m][n][count].real()<<"\t"<<Phi[m][n][count].imag()<<"\t";
                            real_phi_fields_amp<<std::setprecision(16)<<std::norm(Phi[m][n][count])<<"\t";
                            real_chi_fields<<std::setprecision(16)<<iPhiDot[m][n][count].real()<<"\t"<<iPhiDot[m][n][count].imag()<<"\t";
                        }
                    }
                    real_phi_fields<<"\n";
                    real_phi_fields_amp<<"\n";
                    real_chi_fields<<"\n";
    
                    kinetic_dot_term<<(count+1)*h<<"\t"<<std::setprecision(16)<<Hamiltonian_Density_lin_dot[count]<<"\n";
                    kinetic_grad_term<<(count+1)*h<<"\t"<<std::setprecision(16)<<Hamiltonian_Density_lin_grad[count]<<"\n";
                    interaction_term<<(count+1)*h<<"\t"<<std::setprecision(16)<<Hamiltonian_Density_nonlin[count]<<"\n";
                    hamiltonian_density_file<<(count+1)*h<<"\t"<<std::setprecision(16)<<hamiltonian_density[count]<<"\n";
                    momentum_density_file<<(count+1)*h<<"\t"<<std::setprecision(16)<<momentum_density[count]<<"\n";
    
                    normal_mode_file<<(count+1)*h<<"\t";
                    for(int i=0;i<9;i++)
                    {
                        normal_mode_file<<std::setprecision(16)<<NormalModes[i][count]<<"\t";
                    }
                    normal_mode_file<<"\n";
    
                    charge_den_so6<<(count+1)*h<<"\t";
                    for(int i=0;i<5;i++)
                    {
                        for(int j=i+1;j<6;j++)
                        {
                            charge_den_so6<<std::setprecision(16)<<charge_density_so6[(5*i)+j-1-((i+1)*i/2)][count]<<"\t";
                        }
                    }
                    charge_den_so6<<"\n";
    
                    charge_den_su2<<(count+1)*h<<"\t";
                    for(int i=0;i<2;i++)
                    {
                        for(int j=i+1;j<3;j++)
                        {
                            charge_den_su2<<std::setprecision(16)<<charge_density_su2[i+j-1][count]<<"\t";
                        }
                    }
                    charge_den_su2<<"\n";
                }
                real_phi_fields.close();
                real_phi_fields_amp.close();
                real_chi_fields.close();
                kinetic_dot_term.close();
                kinetic_grad_term.close();
                interaction_term.close();
                hamiltonian_density_file.close();
                normal_mode_file.close();
                charge_den_so6.close();
                charge_den_su2.close();
    
            }
            else
            {
                cout << "Unable to open file";
            }
        }
        
        //--------------Storing total energy and total charges for all time steps----------------//
        for(int i=0;i<5;i++)
        {
            for(int j=i+1;j<6;j++)
            {
                total_charge_so6_file[(5*i)+j-1-((i+1)*i/2)]<<(k-1)*tau/2<<"\t"<<std::setprecision(16)<<charge_so6[(5*i)+j-1-((i+1)*i/2)]<<endl;
            }
        }
    
        for(int i=0;i<3;i++)
        {
            total_charge_su2_file[i]<<(k-1)*tau/2<<"\t"<<std::setprecision(16)<<charge_su2[i]<<endl;
        }
    
        myfile_energy<<(k-1)*tau/2<<"\t"<<std::setprecision(16)<< Energy_total<<endl;
        
        myfile_momentum<<(k-1)*tau/2<<"\t"<<std::setprecision(16)<< Momentum_total<<endl;
    
        //--------------Field Evolution----------------//
        for(int l=0;l<4;l++)
        {
            #pragma omp parallel for num_threads(20)
            for(int m=0;m<3;m++)
            {
                for(int n=0;n<3;n++)
                {
                    for(int j=0;j<N;j++)
                    {
                        if(j==0)
                        {
                            Del2Phi[m][n][j] = (Phi[m][n][j+1] - 2.*Phi[m][n][j] + Phi[m][n][N-1])*pow(h,-2.);
                        }
                        else if(j==N-1)
                        {
                            Del2Phi[m][n][j] = (Phi[m][n][0] - 2.*Phi[m][n][j] + Phi[m][n][j-1])*pow(h,-2.);
                        }
                        else
                        {
                            Del2Phi[m][n][j] = (Phi[m][n][j+1] - 2.*Phi[m][n][j] + Phi[m][n][j-1])*pow(h,-2.);
                        }
                    }
        
                    for(int j=0;j<N;j++)
                    {
                        NonlinearPhi[m][n][j] = 0.;
                        for(int A=0;A<3;A++)
                        {
                            for(int B=0;B<3;B++)
                            {
                                NonlinearPhi[m][n][j] += 2.*(Phi[m][A][j]*Phi[B][n][j]*(std::conj(Phi[B][A][j])) + Phi[m][A][j]*(std::conj(Phi[B][n][j]))*Phi[B][A][j] - 2.*Phi[m][n][j]*(std::norm(Phi[B][A][j])));
                            }
                        }
                    }
                }
            }
            
            if(l==0)
            {
                for(int m=0;m<3;m++)
                {
                    for(int n=0;n<3;n++)
                    {
                        for(int i=0;i<N;i++)
                        {
                            Phi[m][n][i] = Phi[m][n][i] - c_n[l]*b*tau*iPhiDot[m][n][i];
                        }
                    }
                }
            }
            else
            {
                for(int m=0;m<3;m++)
                {
                    for(int n=0;n<3;n++)
                    {
                        for(int i=0;i<N;i++)
                        {
                            iPhiDot[m][n][i] = iPhiDot[m][n][i] + d_n[l-1]*b*tau*(Del2Phi[m][n][i] + NonlinearPhi[m][n][i]);
                            Phi[m][n][i] = Phi[m][n][i] - c_n[l]*b*tau*iPhiDot[m][n][i];
                        }
                    }
                }
            }
        }
    }
    
    for(int i=0;i<5;i++)
    {
        for(int j=i+1;j<6;j++)
        {
            total_charge_so6_file[(5*i)+j-1-((i+1)*i/2)].close();
        }
    }
    for(int i=0;i<3;i++)
    {
        total_charge_su2_file[i].close();
    }
    myfile_energy.close();
    myfile_momentum.close();
}

