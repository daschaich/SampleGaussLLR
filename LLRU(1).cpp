#include <iostream>
#include <iomanip>
#include <fstream>
#include <array>
#include <string>
#include <map>
#include <random>
#include <cmath>
#include <time.h>
#include <math.h>
#include <complex>

using namespace std;

const double PI = 2*acos(0.0);

mt19937 rng( random_device{}() ); 


const int d=3; //Dimension
const int n=6; //length of the lattice
const int V=pow(n,d); //Number of lattice sites

const int N_TH=100;
//const int nequi=500;
const int nskip=1;
//const int nskip=50;
//const int nmeas=30;
const int N_SW=50;


const double Emin = 0.4*3;
const double Emax = 0.6*3;
const double deltaE =0.1;
const int N_RM=1000;



void neibinit(int neighbour[V][2*d]);
void filllinks(double theta[V][d],int start);
void metropolisupdate(double theta[V][d],double a,double beta, int neighbour[V][2*d], int nsweeps);
void metropolisupdateconst(double theta[V][d],double a,double beta, int neighbour[V][2*d], int nsweeps,double E,double delta);
//void heatbathupdate(double theta[V][d], double beta, int neighbour[V][2*d], int nsweeps);
double calcaction(double theta[V][d], int neighbour[V][2*d],double beta);
double plaqu(double theta[V][d], int neighbour[V][2*d], int i, int mu, int nu);
double plaquchange(double theta[V][d], int neighbour[V][2*d],int i, int mu, int nu, double offer, int sitechange, int muchange);






int main()
{
	rng.seed(time(NULL));
	
	
	
	int neighbour[V][2*d];
	neibinit(neighbour);
	double theta[V][d];
	double S[N_SW];
	double Dummy[N_SW];
	bool Einterval = false;
	
	double E_i;
	double a[int((Emax-Emin)/deltaE)][N_RM];
	
	double Energysum;
	double DeltaE_i;
	
	double beta=1;
	
	for(int i=0;i<(Emax-Emin)/deltaE;i++)
	{
		E_i=Emin+i*deltaE;
		a[i][0]=1;
		
		Energysum = 0;
		
		for(int n=0;n<N_RM;n++)
		{
			//Einterval=false;
			filllinks(theta,3);
			while(Einterval==false)
			{
				metropolisupdate(theta,a[i][n],beta,neighbour,1);
				if(calcaction(theta,neighbour,beta)<=(E_i+deltaE) && calcaction(theta,neighbour,beta)>=E_i)
				{
					Einterval=true;
					cout << "found Energy interval  " << calcaction(theta,neighbour,beta) <<  endl;
				}
			}
			metropolisupdateconst(theta,a[i][n],beta,neighbour,N_TH,E_i,deltaE);
			cout << "here " << calcaction(theta,neighbour,beta) << endl;
			for(int k = 0;k<N_SW;k++)
			{
				Energysum = Energysum + calcaction(theta,neighbour,beta);
				metropolisupdateconst(theta,a[i][n],beta,neighbour,nskip,E_i,deltaE);
			}
			
			DeltaE_i = 1/(N_SW)*Energysum-E_i-0.5*deltaE;
			a[i][n+1] = a[i][n] - 12*DeltaE_i/((n+1)*deltaE*deltaE);
			cout << "i:" << i << "n: " << n << "a: " << a[i][n] << "   " << a[i][n+1]-a[i][n] <<  endl;
		}
	}
	
	
}



void  neibinit ( int  neib[V][2*d] )
{
	int  i ,  n1p , n1m,  n2p , n2m , n3p , n3m;
	
	for (  int  n1=0; n1<n ;  n1++)
	{
		for (  int  n2=0; n2<n;  n2++)
		{
			for(int n3=0;n3<n;n3++)
			{
				
			
				n1p = n1+1;
				n1m = n1-1;
				n2p = n2+1;
				n2m = n2-1;
				n3p = n3+1;
				n3m = n3-1;
			
				// periodic boundary condition
				if (n1p == n)
				{
					n1p = 0;
				}
				if (n1m ==-1)
				{
					n1m = n-1;
				}
				if (n2p == n)
				{
					n2p = 0;
				}
				if (n2m ==-1)
				{
					n2m = n-1;
				}
				if	(n3p == n)
				{
					n3p = 0;
				}
				if	(n3m ==-1)
				{
					n3m = n-1;
				}
				
				i = n1 + n*n2 + n*n*n3 ;// z e n tr a l e r  punkt  ( aufger .  Index )
				neib[i][0] = n1p + n*n2 + n*n*n3;// +1dir
				neib[i][1] = n1 + n*n2p + n*n*n3;// +2dir
				neib[i][2] = n1 + n*n2 + n*n*n3p;// +3dir
				neib[i][3] = n1m + n*n2 + n*n*n3;// -1dir
				neib[i][4] = n1 + n*n2m + n*n*n3;// -2dir
				neib[i][5] = n1 + n*n2 + n*n*n3m;// -3dir
				
			}
		}
	}
}

void filllinks(double theta[V][d],int start)
{
	uniform_real_distribution<> dist(-PI, PI);
	
	if(start == 1)
	{
		for(int i=0;i<V;i++)
		{
			for(int j=0;j<d;j++)
			{
				theta[i][j]=dist(rng);
			}
		}
	}
	
	if(start == 2)
	{
		for(int i=0;i<V;i++)
		{
			for(int j=0;j<d;j++)
			{
				theta[i][j]=-PI;
			}
		}
	}
	
	if(start == 3)
	{
		for(int i=0;i<V;i++)
		{
			for(int j=0;j<d;j++)
			{
				theta[i][j]=0.0;
			}
		}
	}
}

void metropolisupdate(double theta[V][d],double a,double beta, int neighbour[V][2*d], int nsweeps)
{
	uniform_real_distribution<> dist(-PI, PI);
	//uniform_real_distribution<> dist(0, 2*PI);
	uniform_real_distribution<> dist1(0,1);
	double offer; //new offered linkvar
	double rho, r; //Probabilities for acceptance/rejection
	
	
	for(int n=0;n<nsweeps;n++)
	{
		for(int i=0;i<V;i++)
		{
			for(int j=0;j<3;j++)
			{
				offer = dist(rng);
				
				if(j==0)
				{
					//rho = exp(beta*(-cos(offer+theta[neighbour[i][0]][1]-theta[neighbour[i][1]][0]-theta[i][1])-cos(theta[neighbour[i][3]][0]+theta[neighbour[neighbour[i][0]][3]][1]-offer-theta[neighbour[i][3]][1])+cos(theta[i][0]+theta[neighbour[i][0]][1]-theta[neighbour[i][1]][0]-theta[i][1])+cos(theta[neighbour[i][3]][0]+theta[neighbour[neighbour[i][0]][3]][1]-theta[i][0]-theta[neighbour[i][3]][1])));
					
					rho = exp(-a*beta*(-cos(plaquchange(theta,neighbour,i,0,1,offer,i,j))-cos(plaquchange(theta,neighbour,i,0,2,offer,i,j))-cos(plaquchange(theta,neighbour,neighbour[i][4],0,1,offer,i,j))-cos(plaquchange(theta,neighbour,neighbour[i][5],0,2,offer,i,j))+cos(plaqu(theta,neighbour,i,0,1))+cos(plaqu(theta,neighbour,i,0,2))+cos(plaqu(theta,neighbour,neighbour[i][4],0,1))+cos(plaqu(theta,neighbour,neighbour[i][5],0,2))));
					
					//cout << beta*(cos(offer+theta[neighbour[i][0]][1]-theta[neighbour[i][1]][0]-theta[i][1])+cos(theta[neighbour[i][3]][0]+theta[neighbour[neighbour[i][0]][3]][1]-offer-theta[neighbour[i][3]][1]) - cos(theta[i][0]+theta[neighbour[i][0]][1]-theta[neighbour[i][1]][0]-theta[i][1])-cos(theta[neighbour[i][3]][0]+theta[neighbour[neighbour[i][0]][3]][1]-theta[i][0]-theta[neighbour[i][3]][1])) << endl;	
				}
				else if(j==1)
				{
					//rho = exp(beta*(-cos(theta[i][0]+theta[neighbour[i][0]][1]-theta[neighbour[i][1]][0]-offer)-cos(theta[neighbour[i][2]][0]+offer-theta[neighbour[neighbour[i][1]][2]][0]-theta[neighbour[i][2]][1])+cos(theta[i][0]+theta[neighbour[i][0]][1]-theta[neighbour[i][1]][0]-theta[i][1])+cos(theta[neighbour[i][2]][0]+theta[i][1]-theta[neighbour[neighbour[i][1]][2]][0]-theta[neighbour[i][2]][1])));
				
					rho = exp(-a*beta*(-cos(plaquchange(theta,neighbour,i,0,1,offer,i,j))-cos(plaquchange(theta,neighbour,i,1,2,offer,i,j))-cos(plaquchange(theta,neighbour,neighbour[i][3],0,1,offer,i,j))-cos(plaquchange(theta,neighbour,neighbour[i][5],1,2,offer,i,j))+cos(plaqu(theta,neighbour,i,0,1))+cos(plaqu(theta,neighbour,i,1,2))+cos(plaqu(theta,neighbour,neighbour[i][3],0,1))+cos(plaqu(theta,neighbour,neighbour[i][5],1,2))));
					
				}
				
				else
				{
					rho = exp(-a*beta*(-cos(plaquchange(theta,neighbour,i,0,2,offer,i,j))-cos(plaquchange(theta,neighbour,i,1,2,offer,i,j))-cos(plaquchange(theta,neighbour,neighbour[i][3],0,2,offer,i,j))-cos(plaquchange(theta,neighbour,neighbour[i][4],1,2,offer,i,j))+cos(plaqu(theta,neighbour,i,0,2))+cos(plaqu(theta,neighbour,i,1,2))+cos(plaqu(theta,neighbour,neighbour[i][3],0,2))+cos(plaqu(theta,neighbour,neighbour[i][4],1,2))));
					
				}
			 	r = dist1(rng);
			 	//cout << r << " " << rho << endl;
			 	if(r<=rho)
			 	{
			 		theta[i][j]=offer;
			 		//cout << "change accepted" << endl;
				}
				
			}
			
		}
	}
	
}



void metropolisupdateconst(double theta[V][d],double a,double beta, int neighbour[V][2*d], int nsweeps,double E,double delta)
{
	uniform_real_distribution<> dist(-PI, PI);
	//uniform_real_distribution<> dist(0, 2*PI);
	uniform_real_distribution<> dist1(0,1);
	double offer; //new offered linkvar
	double rho, r; //Probabilities for acceptance/rejection
	double save;
	int counter1=0;
	int counter2=0;
	
	for(int n=0;n<nsweeps;n++)
	{
		for(int i=0;i<V;i++)
		{
			for(int j=0;j<3;j++)
			{
				offer = dist(rng);
				
				if(j==0)
				{
					//rho = exp(beta*(-cos(offer+theta[neighbour[i][0]][1]-theta[neighbour[i][1]][0]-theta[i][1])-cos(theta[neighbour[i][3]][0]+theta[neighbour[neighbour[i][0]][3]][1]-offer-theta[neighbour[i][3]][1])+cos(theta[i][0]+theta[neighbour[i][0]][1]-theta[neighbour[i][1]][0]-theta[i][1])+cos(theta[neighbour[i][3]][0]+theta[neighbour[neighbour[i][0]][3]][1]-theta[i][0]-theta[neighbour[i][3]][1])));
					
					rho = exp(-a*beta*(-cos(plaquchange(theta,neighbour,i,0,1,offer,i,j))-cos(plaquchange(theta,neighbour,i,0,2,offer,i,j))-cos(plaquchange(theta,neighbour,neighbour[i][4],0,1,offer,i,j))-cos(plaquchange(theta,neighbour,neighbour[i][5],0,2,offer,i,j))+cos(plaqu(theta,neighbour,i,0,1))+cos(plaqu(theta,neighbour,i,0,2))+cos(plaqu(theta,neighbour,neighbour[i][4],0,1))+cos(plaqu(theta,neighbour,neighbour[i][5],0,2))));
					
					//cout << beta*(cos(offer+theta[neighbour[i][0]][1]-theta[neighbour[i][1]][0]-theta[i][1])+cos(theta[neighbour[i][3]][0]+theta[neighbour[neighbour[i][0]][3]][1]-offer-theta[neighbour[i][3]][1]) - cos(theta[i][0]+theta[neighbour[i][0]][1]-theta[neighbour[i][1]][0]-theta[i][1])-cos(theta[neighbour[i][3]][0]+theta[neighbour[neighbour[i][0]][3]][1]-theta[i][0]-theta[neighbour[i][3]][1])) << endl;	
				}
				else if(j==1)
				{
					//rho = exp(beta*(-cos(theta[i][0]+theta[neighbour[i][0]][1]-theta[neighbour[i][1]][0]-offer)-cos(theta[neighbour[i][2]][0]+offer-theta[neighbour[neighbour[i][1]][2]][0]-theta[neighbour[i][2]][1])+cos(theta[i][0]+theta[neighbour[i][0]][1]-theta[neighbour[i][1]][0]-theta[i][1])+cos(theta[neighbour[i][2]][0]+theta[i][1]-theta[neighbour[neighbour[i][1]][2]][0]-theta[neighbour[i][2]][1])));
				
					rho = exp(-a*beta*(-cos(plaquchange(theta,neighbour,i,0,1,offer,i,j))-cos(plaquchange(theta,neighbour,i,1,2,offer,i,j))-cos(plaquchange(theta,neighbour,neighbour[i][3],0,1,offer,i,j))-cos(plaquchange(theta,neighbour,neighbour[i][5],1,2,offer,i,j))+cos(plaqu(theta,neighbour,i,0,1))+cos(plaqu(theta,neighbour,i,1,2))+cos(plaqu(theta,neighbour,neighbour[i][3],0,1))+cos(plaqu(theta,neighbour,neighbour[i][5],1,2))));
					
				}
				
				else
				{
					rho = exp(-a*beta*(-cos(plaquchange(theta,neighbour,i,0,2,offer,i,j))-cos(plaquchange(theta,neighbour,i,1,2,offer,i,j))-cos(plaquchange(theta,neighbour,neighbour[i][3],0,2,offer,i,j))-cos(plaquchange(theta,neighbour,neighbour[i][4],1,2,offer,i,j))+cos(plaqu(theta,neighbour,i,0,2))+cos(plaqu(theta,neighbour,i,1,2))+cos(plaqu(theta,neighbour,neighbour[i][3],0,2))+cos(plaqu(theta,neighbour,neighbour[i][4],1,2))));
					
				}
			 	r = dist1(rng);
			 	//cout << r << " " << rho << endl;
			 	counter2 = counter2 +1;
			 	if(r<=rho)
			 	{
			 		save = theta[i][j];
			 		theta[i][j]=offer;
			 		if(calcaction(theta,neighbour,beta)<(E+delta) && calcaction(theta,neighbour,beta)>=E)
			 		{
			 			//cout << "offer accepted " << i <<   endl;
						counter1 = counter1 +1;	
					}
					else
					{
						theta[i][j]=save;
						j=j-1;
						//cout << offer <<  "    " << calcaction(theta,neighbour,beta) << "i: " << i << " j: " << j+1 << endl;
					}
					
			 	
				}
				
			}
			
		}
	}
	cout << counter1 << " " << counter2 <<endl;
	
}


double calcaction(double theta[V][d],  int neighbour[V][2*d],double beta)
{
	double S = 0;
	for(int i = 0;i<V;i++)
	{
		for(int j=0;j<(d-1);j++)
		{
			for(int k =j+1;k<d;k++)
			{
				if(j<k)
				{
					S = S + cos(plaqu(theta,neighbour,i,j,k));
				}
			}
		}
		
	}
	S = beta*S/double(V);
	return S;
}

double plaqu(double theta[V][d], int neighbour[V][2*d], int i, int mu, int nu)
{
	return (theta[i][mu] + theta[neighbour[i][mu]][nu] - theta[neighbour[i][nu]][mu] - theta[i][nu]);
}

double plaquchange(double theta[V][d], int neighbour[V][2*d],int i, int mu, int nu, double offer, int sitechange, int muchange)
{
	double changedplaqu;
	double save = theta[sitechange][muchange];
	theta[sitechange][muchange] = offer;
	
	changedplaqu = plaqu(theta,neighbour,i,mu,nu);
	
	theta[sitechange][muchange] = save;
	
	return changedplaqu;
}


