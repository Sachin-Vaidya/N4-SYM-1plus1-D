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
//#include<sys/types.h>
//#include<sys/stat.h>

#include<iomanip>


using namespace std;
typedef complex<double> dcomp;

int main()
{
double T;
double Tmin = 0.;
double Tmax = 200.;
double Zmin = 0.;
double Zmax = 20.;
double h = 0.002;
double tau = 0.001;
double gamma = tau/pow(h,2);
int N = (Zmax-Zmin)/h;
//double tau = gamma*pow(h,2);
int N_t = (Tmax-Tmin)/(tau);
double pi = 3.141592653589793238;
double hamiltonian_density[N];
double charge_density_so6[15][N], charge_density_su2[3][N];
double charge_so6[15], charge_su2[3];
double f[6][3][2*N];
double beta=0.4;
double c_n[4]={0.6756,-0.1756,-0.1756,0.6756};
double d_n[4]={1.3512,-1.7024,1.3512,0};
/*double c1=0.6756;
double c4=0.6756;
double c2=-0.1756;
double c3=-0.1756;
double d1=2*c1;
double d3=2*c1;
double d2=2*(c2-c1);*/


double ETA;


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

double coefficientRe = 1/sqrt(2);
dcomp coefficient(coefficientRe,0);

cout<<"cccRe = "<<cccRe<<endl;
cout<<"uRe = "<<uRe<<endl;
cout<<"coefficientRe = "<<coefficientRe<<endl;

dcomp p3(0,0);

dcomp co(0.25,0);

double f0=1/sqrt(8);
cout<<setprecision(16)<<"f0 = "<<f0<<endl;
double one_third=pow(3,-1);
cout<<"one_third = "<<one_third<<endl;

double u_p[3][3][4*N];
double u0[3][3][4*N];
double u1[3][3][4*N];
double v0[3][3][4*N];
double v1[3][3][4*N];
double u3[9][N];
double U0[3][3][4*N];
double V0[3][3][4*N];

double u00[3][3][2*N];
double u000[3][3][2*N];

double du0[3][3][2*N];

double d_du0[3][3][2*N];

double du0_nm[9][2*N];

double d_du0_nm[9][2*N];


double Energy_density;

double Hamiltonian_Density_x_derivative_term[3][3][2*N];
double Hamiltonian_Density_nonlin[N];
double Hamiltonian_Density_lin_dot[N];
double Hamiltonian_Density_lin_grad[N];

//double eta[N];

double eta;
double indie_var;
double SD=100.;

double in=(Zmax-Zmin);//(2*pi);
cout<<"in = "<<in<<endl;

int integer = in;
cout<<"integer = "<<integer<<endl;

std::normal_distribution<double> nd2(0.,SD);
std::normal_distribution<double> nd(7.,2.);
//--------------Setting up fluctuation part of the initial conditions for normal modes----------------//
unsigned seed0=std::chrono::system_clock::now().time_since_epoch().count();
std::default_random_engine da(seed0);

for(int m=0;m<9;m++)
{
    double amp_u[150],amp_v[150];
    for(int i=1;i<in;i++)
    {
        /*double x1 = (rand() % 100 + 1);
        double y1 = (rand() % 100 + 1);
        double x2 = (rand() % 100 + 1);
        double y2 = (rand() % 100 + 1);*/
        double x1 = nd(da);//exp(-pow(nd2(da)/SD,2)/2)/(SD*sqrt(2*pi));
        double y1 = nd(da);//exp(-pow(nd2(da)/SD,2)/2)/(SD*sqrt(2*pi));
        double x2 = nd(da);//(-pow(nd2(da)/SD,2)/2)/(SD*sqrt(2*pi));
        double y2 = nd(da);//exp(-pow(nd2(da)/SD,2)/2)/(SD*sqrt(2*pi));
        amp_u[i-1]=x1/100;
        amp_v[i-1]=x2/100;
        amp_u[integer+i-1]=y1/100;
        amp_v[integer+i-1]=y2/100;
    }

    for(int j=0;j<N;j++)
    {

        du0_nm[m][j]=0;
        du0_nm[m][N+j]=0;
        d_du0_nm[m][j]=0;
        d_du0_nm[m][N+j]=0;
        for(int nn=1;nn<in;nn++)
        {
            if(m==0 || m==1 || m==3 || m==4 || m==5)
            {
                if(nn<=(in/(2*pi)))
                {
                    omega=b*sqrt(abs(0.5*((2*pow((2*pi*nn/(in)),2))+3-sqrt(9+(16*pow((2*pi*nn/(in)),2))))));
                    //cout<<"omega = "<<omega<<endl;
                }
                else
                {
                    omega=sqrt((0.5*((2*pow((2*pi*nn/(in)),2))+3-sqrt(9+(16*pow((2*pi*nn/(in)),2))))));
                    //cout<<"omega = "<<omega<<endl;
                }
            }
            else if(m==2)
            {
                omega=sqrt(pow((2*pi*nn/(in)),2)+3-sqrt(9+(4*pow((2*pi*nn/(in)),2))));
                //cout<<"omega = "<<abs(omega+p2)<<endl;
            }
            else
            {
                omega=sqrt(1+pow((2*pi*nn/(in)),2))-1;
                //cout<<"omega = "<<abs(omega+p2)<<endl;
            }

            dcomp amp_uc1(amp_u[nn-1],0);
            dcomp amp_uc2(amp_u[integer+nn-1],0);
            dcomp amp_vc1(amp_v[nn-1],0);
            dcomp amp_vc2(amp_v[integer+nn-1],0);
            dcomp phase(2*pi*h*j*nn/(in),0);
            /*du0_nm[m][j]+=real(((amp_uc1+(b*amp_uc2))*cos(abs(b*((phase)))))-(conj(amp_vc1+(b*amp_vc2))*cos(abs(b*((phase))))));
            d_du0_nm[m][j]+=real(((omega+p2)*(amp_uc1+(b*amp_uc2))*cos(abs(b*((phase)))))+((omega+c)*conj(amp_vc1+(b*amp_vc2))*cos(abs(b*((phase))))));//((amp_uc1*exp(b*((phase))))+(conj(amp_vc1)*exp(c*b*((phase)))));
            du0_nm[m][N+j]+=imag(((amp_uc1+(b*amp_uc2))*cos(abs(b*((phase)))))-(conj(amp_vc1+(b*amp_vc2))*cos(abs(b*((phase))))));
            d_du0_nm[m][N+j]+=imag(((omega+p2)*(amp_uc1+(b*amp_uc2))*cos(abs(b*((phase)))))+((omega+c)*conj(amp_vc1+(b*amp_vc2))*cos(abs(b*((phase))))));*///((amp_u(nn-1)*exp(b*((phase))))+(conj(amp_v(nn-1))*exp(c*b*((phase)))));*/

            du0_nm[m][j]+=real(((amp_uc1+(b*amp_uc2))*exp(b*((phase))))-(conj(amp_vc1+(b*amp_vc2))*exp(c*b*((phase)))));
            d_du0_nm[m][j]+=real(((omega+p2)*(amp_uc1+(b*amp_uc2))*exp(b*((phase))))+((omega+c)*conj(amp_vc1+(b*amp_vc2))*exp(c*b*((phase)))));//((amp_uc1*exp(b*((phase))))+(conj(amp_vc1)*exp(c*b*((phase)))));
            du0_nm[m][N+j]+=imag(((amp_uc1+(b*amp_uc2))*exp(b*((phase))))-(conj(amp_vc1+(b*amp_vc2))*exp(c*b*((phase)))));
            d_du0_nm[m][N+j]+=imag(((omega+p2)*(amp_uc1+(b*amp_uc2))*exp(b*((phase))))+((omega+c)*conj(amp_vc1+(b*amp_vc2))*exp(c*b*((phase)))));//((amp_u(nn-1)*exp(b*((phase))))+(conj(amp_v(nn-1))*exp(c*b*((phase)))));
        }

        //cout<<setprecision(32)<<du0_nm[m][j]<<"\t"<<du0_nm[m][N+j]<<endl;
        /*for(int nn=1;nn<in;nn++)
        {

            if(m==0 || m==1 || m==3 || m==4 || m==5)
            {
                if(nn<=(in/(2*pi)))
                {
                    omega=b*sqrt(abs(0.5*((2*pow((2*pi*nn/(in)),2))+3-sqrt(9+(16*pow((2*pi*nn/(in)),2))))));
                    //cout<<"diff = "<<(pi*nn/(in))-(pi*nn*0.025)<<endl;
                }
                else
                {
                    omega=sqrt((0.5*((2*pow((2*pi*nn/(in)),2))+3-sqrt(9+(16*pow((2*pi*nn/(in)),2))))));
                    //cout<<"omega = "<<omega<<endl;
                }
            }
            else if(m==2)
            {
                omega=sqrt(pow((2*pi*nn/(in)),2)+3-sqrt(9+(4*pow((2*pi*nn/(in)),2))));
                //cout<<"omega = "<<omega<<endl;
            }
            else
            {
                omega=sqrt(1+pow((2*pi*nn/(in)),2))-1;
                //cout<<"omega = "<<omega<<endl;
            }

            if(m==4 || m==7)
            {
                dcomp amp_uc1(amp_u[nn-1],0);
                dcomp amp_uc2(amp_u[integer+nn-1],0);
                dcomp amp_vc1(amp_v[nn-1],0);
                dcomp amp_vc2(amp_v[integer+nn-1],0);
                dcomp phase(2*pi*h*j*nn/(in),0);
                du0_nm[m][j]+=real(((amp_uc1+(b*amp_uc2))*exp(b*((phase))))-(conj(amp_vc1+(b*amp_vc2))*exp(c*b*((phase)))));
                d_du0_nm[m][j]+=real(((omega+p2)*(amp_uc1+(b*amp_uc2))*exp(b*((phase))))+((omega+c)*conj(amp_vc1+(b*amp_vc2))*exp(c*b*((phase)))));//((amp_uc1*exp(b*((phase))))+(conj(amp_vc1)*exp(c*b*((phase)))));
                du0_nm[m][N+j]+=imag(((amp_uc1+(b*amp_uc2))*exp(b*((phase))))-(conj(amp_vc1+(b*amp_vc2))*exp(c*b*((phase)))));
                d_du0_nm[m][N+j]+=imag(((omega+p2)*(amp_uc1+(b*amp_uc2))*exp(b*((phase))))+((omega+c)*conj(amp_vc1+(b*amp_vc2))*exp(c*b*((phase)))));//((amp_u(nn-1)*exp(b*((phase))))+(conj(amp_v(nn-1))*exp(c*b*((phase)))));

            }
        }*/
    }

}


//--------------Figuring out fluctuation part of the initial conditions for fields from above normal modes----------------//
for(int j=0;j<2*N;j++)
{
    du0[0][0][j]=(one_third)*((2*du0_nm[0][j])+du0_nm[1][j]+du0_nm[2][j]);
    du0[1][1][j]=(one_third)*((-1*du0_nm[0][j])+du0_nm[1][j]+du0_nm[2][j]);
    du0[2][2][j]=(one_third)*((-2*du0_nm[1][j])-du0_nm[0][j]+du0_nm[2][j]);
    du0[0][1][j]=(0.5)*(du0_nm[3][j]+du0_nm[6][j]);
    du0[1][0][j]=(0.5)*(du0_nm[3][j]-du0_nm[6][j]);
    du0[1][2][j]=(0.5)*(du0_nm[4][j]+du0_nm[7][j]);
    du0[2][1][j]=(0.5)*(du0_nm[4][j]-du0_nm[7][j]);
    du0[0][2][j]=(0.5)*(du0_nm[5][j]+du0_nm[8][j]);
    du0[2][0][j]=(0.5)*(du0_nm[5][j]-du0_nm[8][j]);

    d_du0[0][0][j]=(one_third)*((2*d_du0_nm[0][j])+d_du0_nm[1][j]+d_du0_nm[2][j]);
    d_du0[1][1][j]=(one_third)*((-1*d_du0_nm[0][j])+d_du0_nm[1][j]+d_du0_nm[2][j]);
    d_du0[2][2][j]=(one_third)*((-2*d_du0_nm[1][j])-d_du0_nm[0][j]+d_du0_nm[2][j]);
    d_du0[0][1][j]=(0.5)*(d_du0_nm[3][j]+d_du0_nm[6][j]);
    d_du0[1][0][j]=(0.5)*(d_du0_nm[3][j]-d_du0_nm[6][j]);
    d_du0[1][2][j]=(0.5)*(d_du0_nm[4][j]+d_du0_nm[7][j]);
    d_du0[2][1][j]=(0.5)*(d_du0_nm[4][j]-d_du0_nm[7][j]);
    d_du0[0][2][j]=(0.5)*(d_du0_nm[5][j]+d_du0_nm[8][j]);
    d_du0[2][0][j]=(0.5)*(d_du0_nm[5][j]-d_du0_nm[8][j]);

    //cout<<setprecision(32)<<one_third<<endl;
}


omp_set_num_threads(3);


//--------------Initial conditions for fields----------------//
////#pragma omp parallel for num_threads(20)
for(int m=0;m<3;m++)
{
    for(int n=0;n<3;n++)
    {
        for(int j=0;j<N;j++)
        {
            if(m==n)
            {
                u0[m][n][j]=((f0)+ 0.0001*(du0[m][n][j]));
                u0[m][n][2*N+j]=((f0)+ 0.0001*(d_du0[m][n][j]));
            }
            else
            {
                u0[m][n][j]=0.0001*(du0[m][n][j]);
                u0[m][n][2*N+j]=0.0001*(d_du0[m][n][j]);
            }
            u0[m][n][N+j]=0.0001*(du0[m][n][N+j]);
            u0[m][n][3*N+j]=0.0001*(d_du0[m][n][N+j]);

            //cout<<"u0["<<m<<"]["<<n<<"]["<<j<<"] = "<<u0[m][n][j]<<endl;
            //cout<<"u0["<<m<<"]["<<n<<"]["<<N+j<<"] = "<<u0[m][n][N+j]<<endl;
        }
    }
}

/*for(int i=0;i<N;i++)
{
    eta[i]=nd2(de);
}*/

/*unsigned seed=std::chrono::system_clock::now().time_since_epoch().count();
std::default_random_engine de(seed);
indie_var=nd2(de);

eta=exp(-pow(indie_var/SD,2)/2)/(SD*sqrt(2*pi));

//ETA=0.0001*nd2(de);

//#pragma omp parallel for
for(int m=0;m<3;m++)
{
    for(int n=0;n<3;n++)
    {
        for(int j=0;j<N;j++)
        {
            if(j==0)
            {
                u000[m][n][0]=u0[m][n][1]-u0[m][n][0]-u0[m][n][0]+u0[m][n][N-1];
                u000[m][n][N]=u0[m][n][N+1]-u0[m][n][N]-u0[m][n][N]+u0[m][n][2*N-1];
            }
            else if(j==N-1)
            {
                u000[m][n][N-1]=u0[m][n][0]-u0[m][n][N-1]-u0[m][n][N-1]+u0[m][n][N-2];
                u000[m][n][2*N-1]=u0[m][n][N]-u0[m][n][2*N-1]-u0[m][n][2*N-1]+u0[m][n][2*N-2];
            }
            else
            {
                u000[m][n][j]=u0[m][n][j+1]-u0[m][n][j]-u0[m][n][j]+u0[m][n][j-1];
                u000[m][n][N+j]=u0[m][n][N+j+1]-u0[m][n][N+j]-u0[m][n][N+j]+u0[m][n][N+j-1];
            }
        }
    }
}

for(int m=0;m<3;m++)
{
    for(int n=0;n<3;n++)
    {
        for(int j=0;j<N;j++)
        {
            u00[m][n][j]=0;
            u00[m][n][N+j]=0;
            dcomp a_c(u0[m][n][j],u0[m][n][N+j]);
            for(int A=0;A<3;A++)
            {
                dcomp a2(u0[A][n][j],u0[A][n][N+j]);
                for(int B=0;B<3;B++)
                {
                    dcomp a1(u0[m][B][j],u0[m][B][N+j]);
                    dcomp a3(u0[A][B][j],u0[A][B][N+j]);

                    u00[m][n][j]+= real((a1*a2*conj(a3))+ (a1*a3*conj(a2))-(nonlin*(a_c*a3*conj(a3))));
                    u00[m][n][N+j]+=imag((a1*a2*conj(a3))+ (a1*a3*conj(a2))-(nonlin*(a_c*a3*conj(a3))));
                }
            }
        }
    }
}*/

/*//#pragma omp parallel for
for(int m=0;m<3;m++)
{
    for(int n=0;n<3;n++)
    {
        for(int i=0;i<N;i++)
        {
            V0[m][n][i]=u0[m][n][i];//+(u0[m][n][3*N+i]*tau);
            V0[m][n][N+i]=u0[m][n][N+i];//-(u0[m][n][2*N+i]*tau);
            V0[m][n][2*N+i]=u0[m][n][2*N+i];//-((gamma)*u000[m][n][N+i])-(2*u00[m][n][N+i]*(tau));//-(eta*u0[m][n][N+i]*tau);
            V0[m][n][3*N+i]=u0[m][n][3*N+i];//+((gamma)*u000[m][n][i])+(2*u00[m][n][i]*(tau));//+(eta*u0[m][n][i]*tau);

        }

    }
}*/


//std::string energy_file("myfile_energy"),total_charge_so6("total_charge_so6_"),total_charge_su2("total_charge_su2_"),charge_initial_density_file_so6("charge_initial_density_file_so6_"),charge_initial_density_file_su2("charge_initial_density_file_su2_"),myfile_charge_initial_so6("myfile_charge_initial_so6_"),myfile_charge_initial_su2("myfile_charge_initial_su2_"),charge_mid_density_file_so6("charge_mid_density_file_so6_"),charge_mid_density_file_su2("charge_mid_density_file_su2_"),myfile_charge_mid_so6("myfile_charge_mid_so6_"),myfile_charge_mid_su2("myfile_charge_mid_su2_"),charge_final_density_file_so6("charge_final_density_file_so6_"),charge_final_density_file_su2("charge_final_density_file_su2_"),myfile_charge_final_so6("myfile_charge_final_so6_"),myfile_charge_final_su2("myfile_charge_final_su2_"),interaction_term_at_t("/interaction_term/interaction_term_at_t"),kinetic_dot_term_at_t("/kinetic_dot_term/kinetic_dot_term_at_t"),kinetic_grad_term_at_t("/kinetic_grad_term/kinetic_grad_term_at_t"),charge_den_so6_at_t("/charge_den_so6/charge_den_so6_at_t"),charge_den_su2_at_t("/charge_den_su2/charge_den_su2_at_t"),L_squ,fields_at_t("/fields/fields_at_t_"),fields_amp_at_t("fields_amp/fields_amp_at_t_"),ext(".txt");
//std::string energy_file("myfile_energy"),total_charge_so6("total_charge_so6_"),total_charge_su2("total_charge_su2_"),interaction_term_at_t("interaction_term/interaction_term_at_t"),kinetic_dot_term_at_t("kinetic_dot_term/kinetic_dot_term_at_t"),kinetic_grad_term_at_t("kinetic_grad_term/kinetic_grad_term_at_t"),charge_den_so6_at_t("charge_den_so6/charge_den_so6_at_t"),charge_den_su2_at_t("charge_den_su2/charge_den_su2_at_t"),hamiltonian_density_at_t("hamiltonian_density/hamiltonian_density_at_t"),fields_phi_at_t("fields_phi/fields_phi_at_t_"),fields_chi_at_t("fields_chi/fields_chi_at_t_"),fields_phi_amp_at_t("fields_phi_amp/fields_phi_amp_at_t_"),ext(".txt");
std::string energy_file("myfile_energy"),total_charge_so6("total_charge_so6_"),total_charge_su2("total_charge_su2_"),interaction_term_at_t("interaction_term/interaction_term_at_t"),kinetic_dot_term_at_t("kinetic_dot_term/kinetic_dot_term_at_t"),kinetic_grad_term_at_t("kinetic_grad_term/kinetic_grad_term_at_t"),charge_den_so6_at_t("charge_den_so6/charge_den_so6_at_t"),charge_den_su2_at_t("charge_den_su2/charge_den_su2_at_t"),hamiltonian_density_at_t("hamiltonian_density/hamiltonian_density_at_t"),fields_phi_at_t("fields_phi/fields_phi_at_t_"),fields_chi_at_t("fields_chi/fields_chi_at_t_"),fields_phi_amp_at_t("fields_phi_amp/fields_phi_amp_at_t_"),normal_mode_at_t("normal_mode/normal_mode_at_t"),ext(".txt");

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




//-----------------to store fields at multiple time steps--------------------//

int tt=Tmax;
double kk[tt];
for(int i=0;i<tt;i++)
{
    kk[i]=((i+1)/tau)-1;
    //cout<<"kk["<<i<<"] = "<<kk[i]<<endl;
}



//--------------Time loop for evolution of fields----------------//
for(int k=1;k<N_t-1;k++)
{
    if(k>1)
    {
        for(int j=0;j<N;j++)
        {
            u3[0][j]=pow(abs(U0[0][0][j]+(b*U0[0][0][N+j])+(c*(U0[1][1][j]+(b*U0[1][1][N+j])))),2);
            u3[1][j]=pow(abs(U0[1][1][j]+(b*U0[1][1][N+j])+(c*(U0[2][2][j]+(b*U0[2][2][N+j])))),2);
            u3[2][j]=pow(abs(U0[0][0][j]+(b*U0[0][0][N+j])+U0[1][1][j]+(b*U0[1][1][N+j])+U0[2][2][j]+(b*U0[2][2][N+j])),2);
            u3[3][j]=pow(abs(U0[0][1][j]+(b*U0[0][1][N+j])+U0[1][0][j]+(b*U0[1][0][N+j])),2);
            u3[4][j]=pow(abs(U0[1][2][j]+(b*U0[1][2][N+j])+U0[2][1][j]+(b*U0[2][1][N+j])),2);
            u3[5][j]=pow(abs(U0[0][2][j]+(b*U0[0][2][N+j])+U0[2][0][j]+(b*U0[2][0][N+j])),2);
            u3[6][j]=pow(abs(U0[0][1][j]+(b*U0[0][1][N+j])+(c*(U0[1][0][j]+(b*U0[1][0][N+j])))),2);
            u3[7][j]=pow(abs(U0[1][2][j]+(b*U0[1][2][N+j])+(c*(U0[2][1][j]+(b*U0[2][1][N+j])))),2);
            u3[8][j]=pow(abs(U0[0][2][j]+(b*U0[0][2][N+j])+(c*(U0[2][0][j]+(b*U0[2][0][N+j])))),2);
        }
    }

    //----------------------writing fields, normal modes, charge densities and hamiltonial densities at multiple time steps-------------------------//
    for(int i=0;i<tt;i++)
    {
        if(k==kk[i])
        {
            int T_k=(k+1)*tau;
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
                            real_phi_fields<<std::setprecision(16)<<u0[m][n][count]<<"\t"<<u0[m][n][N+count]<<"\t";
                            real_phi_fields_amp<<std::setprecision(16)<<pow(u0[m][n][count],2)+pow(u0[m][n][N+count],2)<<"\t";
                            real_chi_fields<<std::setprecision(16)<<u0[m][n][2*N+count]<<"\t"<<u0[m][n][3*N+count]<<"\t";
                        }
                    }
                    real_phi_fields<<"\n";
                    real_phi_fields_amp<<"\n";
                    real_chi_fields<<"\n";

                    kinetic_dot_term<<(count+1)*h<<"\t"<<std::setprecision(16)<<Hamiltonian_Density_lin_dot[count]<<"\n";
                    kinetic_grad_term<<(count+1)*h<<"\t"<<std::setprecision(16)<<Hamiltonian_Density_lin_grad[count]<<"\n";
                    interaction_term<<(count+1)*h<<"\t"<<std::setprecision(16)<<Hamiltonian_Density_nonlin[count]<<"\n";
                    hamiltonian_density_file<<(count+1)*h<<"\t"<<std::setprecision(16)<<hamiltonian_density[count]<<"\n";

                    normal_mode_file<<(count+1)*h<<"\t";
                    for(int i=0;i<9;i++)
                    {
                        normal_mode_file<<std::setprecision(16)<<u3[i][count]<<"\t";
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
    }


    for(int m=0;m<3;m++)
    {
        for(int n=0;n<3;n++)
        {
            for(int i=0;i<N;i++)
            {
                f[2*(m+1)-2][n][i]=U0[m][n][i];
                f[2*(m+1)-1][n][i]=U0[m][n][N+i];
                f[2*(m+1)-2][n][N+i]=U0[m][n][3*N+i];
                f[2*(m+1)-1][n][N+i]=-U0[m][n][2*N+i];
                //cout<<"\n"<<f[2*(m+1)-2][n][i]<<"\t"<<U0[m][n][i]<<endl;
            }
        }
    }
    
    /*for(int m=0;m<N;m++)
    {
        eta[m]=nd2(de);
    }*/
    
    unsigned seed2=std::chrono::system_clock::now().time_since_epoch().count();
    std::default_random_engine de2(seed2);
    indie_var=nd2(de2);

    eta=exp(-pow(indie_var/SD,2)/2)/(SD*sqrt(2*pi));
    
    //ETA=0.0001*nd2(de);

    //--------------Field Evolution----------------//
    for(int l=0;l<4;l++)
    {
        //#pragma omp parallel for
        for(int m=0;m<3;m++)
        {
            for(int n=0;n<3;n++)
            {
                for(int j=0;j<N;j++)
                {
                    if(j==0)
                    {
                        u000[m][n][0]=u0[m][n][1]-u0[m][n][0]-u0[m][n][0]+u0[m][n][N-1];
                        u000[m][n][N]=u0[m][n][N+1]-u0[m][n][N]-u0[m][n][N]+u0[m][n][2*N-1];
                    }
                    else if(j==N-1)
                    {
                        u000[m][n][N-1]=u0[m][n][0]-u0[m][n][N-1]-u0[m][n][N-1]+u0[m][n][N-2];
                        u000[m][n][2*N-1]=u0[m][n][N]-u0[m][n][2*N-1]-u0[m][n][2*N-1]+u0[m][n][2*N-2];
                    }
                    else
                    {      
                        u000[m][n][j]=u0[m][n][j+1]-u0[m][n][j]-u0[m][n][j]+u0[m][n][j-1];
                        u000[m][n][N+j]=u0[m][n][N+j+1]-u0[m][n][N+j]-u0[m][n][N+j]+u0[m][n][N+j-1];
                    }
                }
            }
        }

        for(int m=0;m<3;m++)
        {
            for(int n=0;n<3;n++)
            {
                for(int j=0;j<N;j++)
                {   
                    u00[m][n][j]=0;
                    u00[m][n][N+j]=0;
                    dcomp a_c(u0[m][n][j],u0[m][n][N+j]);
                    for(int A=0;A<3;A++)
                    {
                        dcomp a2(u0[A][n][j],u0[A][n][N+j]);
                        for(int B=0;B<3;B++)
                        {
                            dcomp a1(u0[m][B][j],u0[m][B][N+j]);
                            dcomp a3(u0[A][B][j],u0[A][B][N+j]);

                            u00[m][n][j]+= real((a1*a2*conj(a3))+ (a1*a3*conj(a2))-(nonlin*(a_c*a3*conj(a3))));
                            u00[m][n][N+j]+=imag((a1*a2*conj(a3))+ (a1*a3*conj(a2))-(nonlin*(a_c*a3*conj(a3))));
                        }
                    }
                }
            }
        }

        if(l==0)
        {
            //#pragma omp parallel for
            for(int m=0;m<3;m++)
            {
                for(int n=0;n<3;n++)
                {
                    for(int i=0;i<N;i++)
                    {
                        u0[m][n][i]=u0[m][n][i]+(c_n[l]*u0[m][n][3*N+i]*tau);
                        u0[m][n][N+i]=u0[m][n][N+i]-(c_n[l]*u0[m][n][2*N+i]*tau);
                    }
                }
            }
        }
        else
        {
            //#pragma omp parallel for
            for(int m=0;m<3;m++)
            {
                for(int n=0;n<3;n++)
                {
                    for(int i=0;i<N;i++)
                    {
                        u0[m][n][2*N+i]=u0[m][n][2*N+i]-((gamma*d_n[l-1])*u000[m][n][N+i])-(u00[m][n][N+i]*2*tau*d_n[l-1]);//-(eta*u0[m][n][N+i]*tau);
                        u0[m][n][3*N+i]=u0[m][n][3*N+i]+((gamma*d_n[l-1])*u000[m][n][i])+(u00[m][n][i]*2*tau*d_n[l-1]);//+(eta*u0[m][n][i]*tau);
                        u0[m][n][i]=u0[m][n][i]+(c_n[l]*u0[m][n][3*N+i]*tau);
                        u0[m][n][N+i]=u0[m][n][N+i]-(c_n[l]*u0[m][n][2*N+i]*tau);
                    }
                }
            }
        }
    }
    
    
    for(int m=0;m<3;m++)
    {
        for(int n=0;n<3;n++)
        {       
            for(int i=0;i<4*N;i++)
            {
                U0[m][n][i]=u0[m][n][i];
                u1[m][n][i]=U0[m][n][i];
            }
        }
    }
    
    /*//#pragma omp parallel for
    for(int m=0;m<3;m++)
    {
        for(int n=0;n<3;n++)
        {       
            if(k==1)
            {
                for(int i=0;i<4*N;i++)
                {
                    U0[m][n][i]=V0[m][n][i];
                }
            }
            else
            {
                for(int i=0;i<2*N;i++)
                {
                    U0[m][n][i]=V0[m][n][i]+((beta/2)*(v1[m][n][i]-(3*V0[m][n][i])+(3*u0[m][n][i])-u_p[m][n][i]));
                    U0[m][n][2*N+i]=V0[m][n][2*N+i];
                }
            }
            
            for(int i=0;i<4*N;i++)
            {
                u1[m][n][i]=U0[m][n][i];
            }

        }
    }*/


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
                        Hamiltonian_Density_nonlin[i] += (2*(pow(u1[m][n][i],2)+pow(u1[m][n][i+N],2))*(pow(u1[o][p][i],2)+pow(u1[o][p][i+N],2)))-(2*((u1[o][p][i]*u1[o][n][i])+(u1[o][p][N+i]*u1[o][n][N+i]))*((u1[m][p][i]*u1[m][n][i])+(u1[m][p][N+i]*u1[m][n][N+i])));
                    }
                }
            }
        }
    }

    ////#pragma omp parallel for num_threads(20)
    for(int i=0;i<2*N;i++)
    {
        if(i==0 || i==N-1 || i==N || i==2*N-1)
        {
            for(int m=0;m<3;m++)
            {
                for(int n=0;n<3;n++)
                {
                    Hamiltonian_Density_x_derivative_term[m][n][0] = u1[m][n][1]-u1[m][n][N-1];
                    Hamiltonian_Density_x_derivative_term[m][n][N-1] = u1[m][n][0]-u1[m][n][N-2];
                    Hamiltonian_Density_x_derivative_term[m][n][N] = u1[m][n][N+1]-u1[m][n][2*N-1];
                    Hamiltonian_Density_x_derivative_term[m][n][2*N-1] = u1[m][n][N]-u1[m][n][2*N-2];
                }
            }
        }
        else
        {
            for(int m=0;m<3;m++)
            {
                for(int n=0;n<3;n++)
                {
                    Hamiltonian_Density_x_derivative_term[m][n][i] = u1[m][n][i+1]-u1[m][n][i-1];
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
                Hamiltonian_Density_lin_dot[i] += ((pow((u1[m][n][2*N+i]),2)+pow((u1[m][n][3*N+i]),2)));
                Hamiltonian_Density_lin_grad[i] += ((pow(Hamiltonian_Density_x_derivative_term[m][n][i],2)+pow(Hamiltonian_Density_x_derivative_term[m][n][N+i],2))/pow(2*h,2));
            }
        }
    }


    for(int i=0;i<N;i++)
    {
        hamiltonian_density[i]=Hamiltonian_Density_lin_dot[i]+Hamiltonian_Density_lin_grad[i]+Hamiltonian_Density_nonlin[i];
    }

    Energy_density=0;
    for(int i=0;i<N;i++)
    {
        Energy_density+=hamiltonian_density[i];
    }
    cout.precision(16);
    //cout<<"Energy = "<<Energy_density*h<<endl;



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
            charge_so6[m]+=charge_density_so6[m][i];
        }
    }

    for(int m=0;m<3;m++)
    {
        charge_su2[m]=0;
        for(int i=0;i<N;i++)
        {
            charge_su2[m]+=charge_density_su2[m][i];
        }
    }



    //-------------Updating Fields----------------//
    ////#pragma omp parallel for num_threads(20)
    /*for(int j=0;j<4*N;j++)
    {
        for(int m=0;m<3;m++)
        {
            for(int n=0;n<3;n++)
            {
                u_p[m][n][j]=u0[m][n][j];
                u0[m][n][j]=U0[m][n][j];
                V0[m][n][j]=v1[m][n][j];
            }
        }
    }*/

    //--------------Storing total energy and total charges for all time steps----------------//
    for(int i=0;i<5;i++)
    {
        for(int j=i+1;j<6;j++)
        {
            total_charge_so6_file[(5*i)+j-1-((i+1)*i/2)]<<k*tau<<"\t"<<std::setprecision(16)<<charge_so6[(5*i)+j-1-((i+1)*i/2)]*h<<endl;
        }
    }

    for(int i=0;i<3;i++)
    {
        total_charge_su2_file[i]<<k*tau<<"\t"<<std::setprecision(16)<<charge_su2[i]*h<<endl;
    }

    myfile_energy<<(k)*tau<<"\t"<<std::setprecision(16)<< Energy_density*h<<endl;

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

}

