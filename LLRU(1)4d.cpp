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

const double beta=1.01;

const int d=4; //Dimension
const int n=8; //length of the lattice
const int V=pow(n,d); //Number of lattice sites

const int N_TH=100; //Thermalization steps
const int nskip=50; //skip steps
const int N_SW=250; //number of sweeps


const double Emin = 0.6*6; //minimum energy that is probed
const double Emax = 1.5*6; //maximum energy that is probed
const double delta =0.0012; //size of energy interval


void neibinit(int neighbour[V][2*d]); //initiates the neighbour field
void filllinks(double theta[V][d],int start); // fills the links of the lattice
void metropolisupdate(double theta[V][d],double a,double beta, int neighbour[V][2*d], int nsweeps); //metropolisupdate unrestricted
void metropolisupdateconst(double theta[V][d],double a,double beta, int neighbour[V][2*d], int nsweeps,double E,double delta); //metropolisupdate within a certain energy interval
//void heatbathupdate(double theta[V][d], double beta, int neighbour[V][2*d], int nsweeps); 
double calcaction(double theta[V][d], int neighbour[V][2*d],double beta);
double plaqu(double theta[V][d], int neighbour[V][2*d], int i, int mu, int nu); //calculates the plaquette
double plaquchange(double theta[V][d], int neighbour[V][2*d],int i, int mu, int nu, double offer, int sitechange, int muchange); //calculates the plaquette with the proposed changed link



random_device rd;
mt19937 gen(rd());
uniform_real_distribution<> dis(0, 1);


int main()
{
	double epsilon = 0.00001; // precision of rob march algorithm
	
	
	
	rng.seed(time(NULL));
	
	int neighbour[V][2*d];
	neibinit(neighbour);
	double theta[V][d];
	double S[N_SW]; //Action
	bool Einterval = false;
	
	
	
	
	

	int Njacknife = 1;
	//int Naverage = 10000;
	
	double x0[int(Emax/delta)];  //lower end of energy interval
	double a[int(Emax/delta)]; 
	double a_i[Njacknife]; 
	double a_i_new;
	int RobMarchcount; //number of iteration of the rob march algorithm
	double Reweightexpect; //Reweighted expectationvalue of the energy
	
	double s2[int(Emax/delta)]; //variance
	double a_i_del[Njacknife];
	//stuff for calculation of density from a
	double A[int(Emax/delta)][Njacknife];
	double Astar;
	double var[int(Emax/delta)];
	double errorjack[int(Emax/delta)];
	double T;
	double rhojack[Njacknife];
	double rho[int(Emax/delta)];
	
	
	ofstream myfilemonopoledens;
	myfilemonopoledens.open("rhovalues.txt");
	
	//int Eint = 0;
	for(int Eint=0;Eint<((Emax-Emin)/delta);Eint++) //scans through the energy intervals
	{
		x0[Eint] = Eint*delta+Emin;
		
		a[Eint] = -0.995;//-1.0018;
		
		//int jcount =0;
		for(int jcount = 0;jcount<Njacknife;jcount++) //generates Njacknife values of a for Jacknife method to get error and average
		{
			//set starting value for a
			if(Eint == 0)
			{
				a_i[jcount] = -0.995+2*epsilon;//-1.0018+2*epsilon;
				
				a_i_new = -0.995;//-1.0018;
			}
			else
			{
				a_i[jcount]=a[Eint-1] - 2*epsilon;
				
				a_i_new = a[Eint-1];
			}
			
			RobMarchcount = 0;
			
			Reweightexpect  = 10.0*epsilon;
			
			filllinks(theta,3);
			
			while(abs(a_i_new-a_i[jcount]) > epsilon) //only stop rob march if a doesnt change much anymore
			{
				
				cout << "difference " << a_i_new-a_i[jcount] << endl;
				
				a_i[jcount] = a_i_new;
				
				
				
				
				
				
				//cout << Einterval << endl;
				
				//finds the right energy intervall starting from the maximum energy (all link variables set to 0)
				if(Einterval == false)
				{
				
					while(Einterval==false)
					{	
						cout << "here Einterval false " << calcaction(theta,neighbour,beta)/6 << endl;
						metropolisupdate(theta,1.0,beta,neighbour,1);
						if(calcaction(theta,neighbour,beta)<=(x0[Eint]+delta) && calcaction(theta,neighbour,beta)>=x0[Eint])
						{
							Einterval=true;
							cout << "found Energy interval  " << calcaction(theta,neighbour,beta)/6 <<  endl;
						}
					}
				}
				
				
				cout << "a " << a_i[jcount] << endl;
				
				/*if(RobMarchcount == 1)
				{
					metropolisupdate(theta,a_i[jcount],beta,neighbour,N_TH);
				}
				*/
				//restricted metropolis updates once the right energy interval is found
				//else
				{
					metropolisupdateconst(theta,a_i[jcount],beta,neighbour,N_TH,x0[Eint],delta);
				}
				
				
				
				cout << "RobMarchcount " <<  RobMarchcount << " Energy " <<  calcaction(theta,neighbour,beta)/6 << endl;
				
				//y1 = 0.5*(erf(8*a_i[jcount]+x0[Eint]/16)+1);
				//y2 = 0.5*(erf(8*a_i[jcount]+(x0[Eint]+delta)/16)+1);
				
				//cout << "y1: " << y1 << endl;
				//cout << "y2: " << y2 << endl;
 				
 				
 				//calculate the reweighted expectation value of the energy
				Reweightexpect=0;
				
				for(int k = 0;k<N_SW;k++)
				{
					Reweightexpect = Reweightexpect + calcaction(theta,neighbour,beta);
					metropolisupdateconst(theta,a_i[jcount],beta,neighbour,nskip,x0[Eint],delta);
				}
				
				
				//cout << "RobMarchcount: " << RobMarchcount << endl;
				Reweightexpect = Reweightexpect/N_SW - x0[Eint] - 0.5*delta;
				
				
				//calculate the new a
				a_i_new = a_i[jcount] + 12/(delta*delta*(RobMarchcount+1))*Reweightexpect;
				cout << "a: " << a_i_new << endl;
				RobMarchcount = RobMarchcount + 1;
			
				
			}
			
			a_i[jcount] = a_i_new;
			a[Eint] = a[Eint] + a_i_new/Njacknife;
			
		}
		//calculation for error and to get rho from a
		s2[Eint] = 0;
		for(int i=0;i<Njacknife;i++)
		{
			s2[Eint]=s2[Eint] + pow(a_i[i]-a[Eint],2);
		}
		
		s2[Eint] = s2[Eint]/(Njacknife-1);
		
		for(int i =0;i<Njacknife;i++)
		{
			for(int j=0;j<Njacknife;j++)
			{
				a_i_del[j] = a_i[j];
			}
			
			A[Eint][i] = 0;
			
			for(int j=0;j<Njacknife;j++)
			{
				if(j!=i)
				{
					A[Eint][i] = A[Eint][i] + a_i[j];
				}
			}
			
			A[Eint][i] = A[Eint][i]/(Njacknife-1);
		}
		
		Astar = 0;
		
		for(int i =0;i<Njacknife;i++)
		{
			Astar = Astar + A[Eint][i];
		}
		
		Astar = Astar/Njacknife;
		
		var[Eint] = 0;
		
		for(int i = 0;i<Njacknife;i++)
		{
			var[Eint] = var[Eint] + pow(Astar-A[Eint][i],2);
		}
		
		var[Eint] = var[Eint]*(Njacknife-1);
		
		errorjack[Eint] = sqrt(var[Eint])/sqrt(Njacknife);
		
		for(int i=0;i<Njacknife;i++)
		{
			T = 0;
			for(int j=0;j<(Eint-1);j++)
			{
				T = T + a[j];
			}
			T = T + A[Eint][i] - a[0]/2 - a[Eint]/2;
			rhojack[i] = exp(delta*T);
		}
		
		rho[Eint]=0;
		
		for(int i =0;i<Njacknife;i++)
		{
			rho[Eint] = rho[Eint] + rhojack[i];
		}
		
		rho[Eint] = rho[Eint]/Njacknife;
		cout << "rho: " << rho[Eint] << endl;
		myfilemonopoledens << rho[Eint] << " \n";
		cout << "Reweightexpect: " << Reweightexpect << endl;
		cout << "a: " << a[Eint] << endl;
		cout << "Energyinterval: " << Eint*delta/6 << endl;
		cout << "Energy: " << calcaction(theta,neighbour,beta)/6 << endl;
		cout << " " << endl;
		
	}
		myfilemonopoledens.close();
	
	
}


void  neibinit ( int  neib[V][2*d] )
{
	int  i ,  n1p , n1m,  n2p , n2m , n3p , n3m, n4p, n4m;
	
	for (  int  n1=0; n1<n ;  n1++)
	{
		for (  int  n2=0; n2<n;  n2++)
		{
			for(int n3=0;n3<n;n3++)
			{
				for(int n4=0;n4<n;n4++)
				{
				
			
					n1p = n1+1;
					n1m = n1-1;
					n2p = n2+1;
					n2m = n2-1;
					n3p = n3+1;
					n3m = n3-1;
					n4p = n4+1;
					n4m = n4-1;
			
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
					if  (n4p == n)
					{
						n4p = 0;
					}
					if  (n4m ==-1)
					{
						n4m = n-1;
					}
				
					i = n1 + n*n2 + n*n*n3 + n*n*n*n4;// z e n tr a l e r  punkt  ( aufger .  Index )
					neib[i][0] = n1p + n*n2 + n*n*n3 + n*n*n*n4;// +1dir
					neib[i][1] = n1 + n*n2p + n*n*n3 + n*n*n*n4;// +2dir
					neib[i][2] = n1 + n*n2 + n*n*n3p + n*n*n*n4;// +3dir
					neib[i][3] = n1 + n*n2 + n*n*n3 + n*n*n*n4p;// +4dir
					neib[i][4] = n1m + n*n2 + n*n*n3 + n*n*n*n4;// -1dir
					neib[i][5] = n1 + n*n2m + n*n*n3 + n*n*n*n4;// -2dir
					neib[i][6] = n1 + n*n2 + n*n*n3m + n*n*n*n4;// -3dir
					neib[i][7] = n1 + n*n2 + n*n*n3 + n*n*n*n4m;// -4dir
				
				}
				
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
	uniform_real_distribution<> dist1(0,1);
	double offer; //new offered linkvar
	double rho, r; //Probabilities for acceptance/rejection
	
	
	
	for(int n=0;n<nsweeps;n++)
	{
		for(int i=0;i<V;i++)
		{
			for(int j=0;j<d;j++)
			{
				offer = dist(rng);
				
				if(j==0)
				{
					
					rho = exp(-beta*(-cos(plaquchange(theta,neighbour,i,0,1,offer,i,j))-cos(plaquchange(theta,neighbour,i,0,2,offer,i,j))-cos(plaquchange(theta,neighbour,i,0,3,offer,i,j))-cos(plaquchange(theta,neighbour,neighbour[i][5],0,1,offer,i,j))-cos(plaquchange(theta,neighbour,neighbour[i][6],0,2,offer,i,j))-cos(plaquchange(theta,neighbour,neighbour[i][7],0,3,offer,i,j))+cos(plaqu(theta,neighbour,i,0,1))+cos(plaqu(theta,neighbour,i,0,2))+cos(plaqu(theta,neighbour,i,0,3))+cos(plaqu(theta,neighbour,neighbour[i][5],0,1))+cos(plaqu(theta,neighbour,neighbour[i][6],0,2))+cos(plaqu(theta,neighbour,neighbour[i][7],0,3))));
					
				}
				
				else if(j==1)
				{
					
					rho = exp(-beta*(-cos(plaquchange(theta,neighbour,i,0,1,offer,i,j))-cos(plaquchange(theta,neighbour,i,1,2,offer,i,j))-cos(plaquchange(theta,neighbour,i,1,3,offer,i,j))-cos(plaquchange(theta,neighbour,neighbour[i][4],0,1,offer,i,j))-cos(plaquchange(theta,neighbour,neighbour[i][6],1,2,offer,i,j))-cos(plaquchange(theta,neighbour,neighbour[i][7],1,3,offer,i,j))+cos(plaqu(theta,neighbour,i,0,1))+cos(plaqu(theta,neighbour,i,1,2))+cos(plaqu(theta,neighbour,i,1,3))+cos(plaqu(theta,neighbour,neighbour[i][4],0,1))+cos(plaqu(theta,neighbour,neighbour[i][6],1,2))+cos(plaqu(theta,neighbour,neighbour[i][7],1,3))));
					
				}
				
				else if(j==2)
				{
					
					rho = exp(-beta*(-cos(plaquchange(theta,neighbour,i,0,2,offer,i,j))-cos(plaquchange(theta,neighbour,i,1,2,offer,i,j))-cos(plaquchange(theta,neighbour,i,2,3,offer,i,j))-cos(plaquchange(theta,neighbour,neighbour[i][4],0,2,offer,i,j))-cos(plaquchange(theta,neighbour,neighbour[i][5],1,2,offer,i,j))-cos(plaquchange(theta,neighbour,neighbour[i][7],2,3,offer,i,j))+cos(plaqu(theta,neighbour,i,0,2))+cos(plaqu(theta,neighbour,i,1,2))+cos(plaqu(theta,neighbour,i,2,3))+cos(plaqu(theta,neighbour,neighbour[i][4],0,2))+cos(plaqu(theta,neighbour,neighbour[i][5],1,2))+cos(plaqu(theta,neighbour,neighbour[i][7],2,3))));
					
				}
				
				else
				{
					rho = exp(-beta*(-cos(plaquchange(theta,neighbour,i,0,3,offer,i,j))-cos(plaquchange(theta,neighbour,i,1,3,offer,i,j))-cos(plaquchange(theta,neighbour,i,2,3,offer,i,j))-cos(plaquchange(theta,neighbour,neighbour[i][4],0,3,offer,i,j))-cos(plaquchange(theta,neighbour,neighbour[i][5],1,3,offer,i,j))-cos(plaquchange(theta,neighbour,neighbour[i][6],2,3,offer,i,j))+cos(plaqu(theta,neighbour,i,0,3))+cos(plaqu(theta,neighbour,i,1,3))+cos(plaqu(theta,neighbour,i,2,3))+cos(plaqu(theta,neighbour,neighbour[i][4],0,3))+cos(plaqu(theta,neighbour,neighbour[i][5],1,3))+cos(plaqu(theta,neighbour,neighbour[i][6],2,3))));
					
				}
			 	r = dist1(rng);
			 
			 	if(r<=rho)
			 	{
			 		theta[i][j]=offer;
			 	
				}
				
			}
			
		}
	}
	
}



void metropolisupdateconst(double theta[V][d],double a,double beta, int neighbour[V][2*d], int nsweeps,double E,double delta)
{
	uniform_real_distribution<> dist(-PI, PI);
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
			for(int j=0;j<d;j++)
			{
				offer = dist(rng);
				
				if(j==0)
				{
					
					rho = exp(-beta*(-cos(plaquchange(theta,neighbour,i,0,1,offer,i,j))-cos(plaquchange(theta,neighbour,i,0,2,offer,i,j))-cos(plaquchange(theta,neighbour,i,0,3,offer,i,j))-cos(plaquchange(theta,neighbour,neighbour[i][5],0,1,offer,i,j))-cos(plaquchange(theta,neighbour,neighbour[i][6],0,2,offer,i,j))-cos(plaquchange(theta,neighbour,neighbour[i][7],0,3,offer,i,j))+cos(plaqu(theta,neighbour,i,0,1))+cos(plaqu(theta,neighbour,i,0,2))+cos(plaqu(theta,neighbour,i,0,3))+cos(plaqu(theta,neighbour,neighbour[i][5],0,1))+cos(plaqu(theta,neighbour,neighbour[i][6],0,2))+cos(plaqu(theta,neighbour,neighbour[i][7],0,3))));
					
				}
				
				else if(j==1)
				{
					
					rho = exp(-beta*(-cos(plaquchange(theta,neighbour,i,0,1,offer,i,j))-cos(plaquchange(theta,neighbour,i,1,2,offer,i,j))-cos(plaquchange(theta,neighbour,i,1,3,offer,i,j))-cos(plaquchange(theta,neighbour,neighbour[i][4],0,1,offer,i,j))-cos(plaquchange(theta,neighbour,neighbour[i][6],1,2,offer,i,j))-cos(plaquchange(theta,neighbour,neighbour[i][7],1,3,offer,i,j))+cos(plaqu(theta,neighbour,i,0,1))+cos(plaqu(theta,neighbour,i,1,2))+cos(plaqu(theta,neighbour,i,1,3))+cos(plaqu(theta,neighbour,neighbour[i][4],0,1))+cos(plaqu(theta,neighbour,neighbour[i][6],1,2))+cos(plaqu(theta,neighbour,neighbour[i][7],1,3))));
					
				}
				
				else if(j==2)
				{
					
					rho = exp(-beta*(-cos(plaquchange(theta,neighbour,i,0,2,offer,i,j))-cos(plaquchange(theta,neighbour,i,1,2,offer,i,j))-cos(plaquchange(theta,neighbour,i,2,3,offer,i,j))-cos(plaquchange(theta,neighbour,neighbour[i][4],0,2,offer,i,j))-cos(plaquchange(theta,neighbour,neighbour[i][5],1,2,offer,i,j))-cos(plaquchange(theta,neighbour,neighbour[i][7],2,3,offer,i,j))+cos(plaqu(theta,neighbour,i,0,2))+cos(plaqu(theta,neighbour,i,1,2))+cos(plaqu(theta,neighbour,i,2,3))+cos(plaqu(theta,neighbour,neighbour[i][4],0,2))+cos(plaqu(theta,neighbour,neighbour[i][5],1,2))+cos(plaqu(theta,neighbour,neighbour[i][7],2,3))));
					
				}
				
				else
				{
					rho = exp(-beta*(-cos(plaquchange(theta,neighbour,i,0,3,offer,i,j))-cos(plaquchange(theta,neighbour,i,1,3,offer,i,j))-cos(plaquchange(theta,neighbour,i,2,3,offer,i,j))-cos(plaquchange(theta,neighbour,neighbour[i][4],0,3,offer,i,j))-cos(plaquchange(theta,neighbour,neighbour[i][5],1,3,offer,i,j))-cos(plaquchange(theta,neighbour,neighbour[i][6],2,3,offer,i,j))+cos(plaqu(theta,neighbour,i,0,3))+cos(plaqu(theta,neighbour,i,1,3))+cos(plaqu(theta,neighbour,i,2,3))+cos(plaqu(theta,neighbour,neighbour[i][4],0,3))+cos(plaqu(theta,neighbour,neighbour[i][5],1,3))+cos(plaqu(theta,neighbour,neighbour[i][6],2,3))));
					
				}
			 	r = dist1(rng);
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
	if(counter1<10)
	{
		cout << "counter1: "<< counter1 << " counter2: " << counter2 << endl;
		//cout << "ratio1: " << double(counter1)/double(counter2) << " ratio2: " << double(counter1)/double(counter3) << endl;
	}
	//cout << counter1 << " " << counter2 << endl;
	
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
	//S = S/double(V);
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

