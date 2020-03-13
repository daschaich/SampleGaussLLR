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

const int d=4; //Dimension
const int n=6; //length of the lattice
const int V=pow(n,d); //Number of lattice sites

const int nequi=1000;
const int nskip=100;
const int nmeas=30;


void neibinit(int neighbour[V][2*d]);
void filllinks(double theta[V][d],int start);
void metropolisupdate(double theta[V][d],double beta, int neighbour[V][2*d], int nsweeps);
void calcaction(double theta[V][d], double S[nmeas], int neighbour[V][2*d], int count);
void calcwilsonloop(double theta[V][d], complex<double> W[nmeas], int neighbour[V][2*d], int count);
double plaqu(double theta[V][d], int neighbour[V][2*d], int i, int mu, int nu);
double plaquchange(double theta[V][d], int neighbour[V][2*d],int i, int mu, int nu, double offer, int sitechange, int muchange);
complex<double> Ulink(double theta);
complex<double> Ulinkplaqu(double theta[V][d], int neighbour[V][2*d],int i, int mu, int nu);
complex<double> Ulinkplaqu2x2(double theta[V][d], int neighbour[V][2*d],int i, int mu, int nu);
double plaqured(double theta[V][d], int neighbour[V][2*d], int i, int mu, int nu);
void calcmonopoledens(double theta[V][d], double M[nmeas], int neighbour[V][2*d], int count);

int main()
{
	
    rng.seed(time(NULL));
	
	double beta=1.02;
	int neighbour[V][2*d];
	double theta[V][d];
	double S[nmeas];
	complex<double> W[nmeas];
	double M[nmeas+1];
	
	
	neibinit(neighbour);
	filllinks(theta,1);
	//ofstream myfileaction;
  	//myfileaction.open ("Actionrandom.txt");
	//ofstream myfileconfig;
	//myfileconfig.open("Configuration.txt");
	//ofstream myfilewilson;
	//myfilewilson.open("Wilsonloop1x1.txt");
	ofstream myfilemonopoledens;
	myfilemonopoledens.open("monopoledensity4d.txt");
	
	
	//calcaction(theta,S,neighbour,0);
	
	//cout << "S start: " << S[0] << endl;
//	myfileaction << S[0] << " \n";
	
while (beta>=0.99)
{	
	
	metropolisupdate(theta,beta,neighbour,nequi); //equilibration
		
	for(int imeas=0;imeas<nmeas;imeas++)
	{
		metropolisupdate(theta,beta,neighbour,nskip);
		calcmonopoledens(theta,M,neighbour,imeas);
		//calcaction(theta,S,neighbour,imeas);
		//calcwilsonloop(theta,W,neighbour,imeas);
		//myfileaction << S[imeas] << " \n";
		//myfilewilson << real(W[imeas]) << " \n";
		myfilemonopoledens << M[imeas] << " \n";
		/*
		for(int i=0;i<V;i++)
		{
			for(int j=0;j<2;j++)
			{
				myfileconfig << theta[i][j] << " \n";
			}
			
		}
		*/
		//cout << imeas << " " << S[imeas] << endl;
			
	}
	cout << beta << endl;
	beta = beta - 0.0015;
	
}
	//myfileaction.close();
	//myfileconfig.close();
	myfilemonopoledens.close();
	
	
}




void  neibinit ( int  neib[V][2*d] )
{
	int  i ,  n1p , n1m,  n2p , n2m , n3p , n3m, n4p, n4m;
	
	for (int n1=0;n1<n;n1++)
	{
		for (int n2=0;n2<n;n2++)
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
						n3m = 0;
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

void metropolisupdate(double theta[V][d],double beta, int neighbour[V][2*d], int nsweeps)
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

void calcaction(double theta[V][d], double S[nmeas], int neighbour[V][2*d], int count)
{
	S[count] = 0;
	for(int i = 0;i<V;i++)
	{
		for(int j=0;j<(d-1);j++)
		{
			for(int k =j+1;k<d;k++)
			{
				if(j<k)
				{
					S[count] = S[count] + cos(plaqu(theta,neighbour,i,j,k));
				}
			}
		}
		
	}
	S[count] = S[count]/double(V);
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

void calcwilsonloop(double theta[V][d], complex<double> W[nmeas], int neighbour[V][2*d], int count)
{
	W[count] = 0.0;
	for(int i=0;i<V;i++)
	{
		for(int j=1;j<d;j++)
		{
			for(int k=0;k<j;k++)
			{
				if(k<j)
				{
					
					W[count] = W[count] + Ulinkplaqu(theta,neighbour,i,j,k);
					
				}
			}
		}
	}
	W[count] = W[count]/double(V)/6.0;
}

complex<double> Ulink(double theta)
{
	return polar(1.0,theta);
}

complex<double> Ulinkplaqu(double theta[V][d], int neighbour[V][2*d],int i, int mu, int nu)
{
	return (Ulink(theta[i][mu])*Ulink(theta[neighbour[i][mu]][nu])*Ulink(-theta[neighbour[i][nu]][mu])*Ulink(-theta[i][nu]));
	//return (Ulink(theta[i][mu])+Ulink(theta[neighbour[i][mu]][nu])+Ulink(-theta[neighbour[i][nu]][mu])+Ulink(-theta[i][nu]));
}

complex<double> Ulinkplaqu2x2(double theta[V][d], int neighbour[V][2*d],int i, int mu, int nu)
{
	return (Ulink(theta[i][mu])*Ulink(theta[neighbour[i][mu]][mu])*Ulink(theta[neighbour[neighbour[i][mu]][mu]][nu])*Ulink(theta[neighbour[neighbour[neighbour[i][mu]][mu]][nu]][nu])*Ulink(-theta[neighbour[neighbour[neighbour[i][mu]][nu]][nu]][mu])*Ulink(-theta[neighbour[neighbour[i][nu]][nu]][mu])*Ulink(-theta[neighbour[i][nu]][nu])*Ulink(-theta[i][nu]));
}

void calcmonopoledens(double theta[V][d], double M[nmeas], int neighbour[V][2*d], int count)
{
	M[count] = 0;
	double total_phase = 0;
	
	for(int i=0;i<V;i++)
	{
		for(int j=0;j<(d-1);j++)
		{
			for(int k=j+1;k<d;k++)
			{
				if(j<k)
				{
					total_phase =  plaqu(theta,neighbour,i,j,k);
					
					if (fabs(total_phase) < PI)
					{
						M[count] = M[count] + 0;
					}
      				    
      			  	else if (total_phase >= PI && total_phase < (3*PI))
      			  	{
      			  		M[count] = M[count] + 1;	
					}
          			
        			else if (total_phase >= (3*PI) && total_phase < (5*PI))
        			{
        				M[count] = M[count] + 2;
					}
          			
        			else if (total_phase >= (5*PI) && total_phase < (7*PI))
        			{
        				M[count] = M[count] + 3;
					}
          				
        			else if (total_phase <= -PI && total_phase > -(3*PI))
        			{
        				M[count] = M[count] + 1;
					}
          			
        			else if (total_phase <= -(3*PI) && total_phase > -(5*PI))
        			{
        				M[count] = M[count] + 2;
					}
          				
        			else if (total_phase <= -(5*PI) && total_phase > -(7*PI))
        			{
        				M[count] = M[count] + 3;
					}
          				
					
					
					
					
					
					
				}
			}
		}
		
	}
	
	//cout << M[count]/(2*PI)/2/V << endl;
	M[count] = M[count]/2/V;
}

double plaqured(double theta[V][d], int neighbour[V][2*d], int i, int mu, int nu)
{
	double plaquette = plaqu(theta,neighbour,i,mu,nu);
	double n = 2*PI*int(plaqu(theta,neighbour,i,mu,nu)/(2*PI));
	//cout << "plaquette initial: " << plaquette << " n: " << n << endl;
	if(abs(plaquette-n)>PI)
	{
		//cout << "we are here" << endl;
		if(plaquette<0)
		{
			//cout << "plaquette negative" << endl;
			plaquette = plaquette - n + 2*PI;
		}
		else if(plaquette>0)
		{
			//cout << "plaquette positive" << endl;
			plaquette = plaquette - n -2*PI;
		}
		return plaquette;
	}
	else
	{
		return (plaquette - n);
	}
	
}
