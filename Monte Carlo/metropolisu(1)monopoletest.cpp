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

using namespace std;

const double PI = 2*acos(0.0);

mt19937 rng( random_device{}() ); 

const int d=3; //Dimension
const int n=5; //length of the lattice
const int V=pow(n,d); //Number of lattice sites

const int nequi=500;
const int nskip=500;
const int nmeas=100;


void neibinit(int neighbour[V][2*d]);
void filllinks(double theta[V][d],int start);
void metropolisupdate(double theta[V][d],double beta, int neighbour[V][2*d], int nsweeps);
void calcaction(double theta[V][d], double S[nmeas], int neighbour[V][2*d], int count);
void calcmonopoledens(double theta[V][d], double M[nmeas], int neighbour[V][2*d], int count);
double plaqu(double theta[V][d], int neighbour[V][2*d], int i, int mu, int nu);
double plaquchange(double theta[V][d], int neighbour[V][2*d],int i, int mu, int nu, double offer, int sitechange, int muchange);
double plaqured(double theta[V][d], int neighbour[V][2*d], int i, int mu, int nu);

int main()
{
	
    rng.seed(time(NULL));
	
	double beta=0.05;
	
	int neighbour[V][2*d];
	double theta[V][d];
	double S[nmeas+1];
	double M[nmeas+1];
	
	
	neibinit(neighbour);
	filllinks(theta,3);
	ofstream myfileaction;
  	myfileaction.open ("Actionorder.txt");
	ofstream myfileconfig;
	myfileconfig.open("Configuration.txt");
	ofstream myfilemonopoledens;
	myfilemonopoledens.open("monopoledensity.txt");
	metropolisupdate(theta,beta,neighbour,nequi); //equilibration
	
	//calcaction(theta,S,neighbour,0);
	//calcmonopoledens(theta,M,neighbour,0);
	
	//cout << "S start: " << S[0] << endl;
	//myfileaction << S[0] << " \n";
	
//while(beta<=2)
{
	metropolisupdate(theta,beta,neighbour,nequi); //equilibration
		
	for(int imeas=0;imeas<nmeas;imeas++)
	{
		metropolisupdate(theta,beta,neighbour,nskip);
		calcaction(theta,S,neighbour,imeas);
		calcmonopoledens(theta,M,neighbour,imeas);
		myfileaction << S[imeas] << " \n";
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
		//cout << imeas << " " << M[imeas] << " " << int(M[imeas]/(2*PI)) <<" " << (M[imeas]/(2*PI)-int(M[imeas]/(2*PI)))*2*PI << endl;
		//cout << imeas << " " << S[imeas] << endl;	
	}
	 
	//beta = beta + 0.1;
	//cout << beta << endl;

}

	myfileaction.close();
	myfileconfig.close();
	myfilemonopoledens.close();
	
	/*
	for(int i=0;i<V;i++)
	{
		for(int j=0;j<2;j++)
		{
			cout << "lattice site i: " << i << " link j: " <<j<< " " <<  theta[i][j] << endl;
		}
		cout << endl;
	}
	*/
	
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
					n3m = 0;
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
	//uniform_real_distribution<> dist(0, 2*PI);
	
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
					
					rho = exp(-beta*(-cos(plaquchange(theta,neighbour,i,0,1,offer,i,j))-cos(plaquchange(theta,neighbour,i,0,2,offer,i,j))-cos(plaquchange(theta,neighbour,neighbour[i][4],0,1,offer,i,j))-cos(plaquchange(theta,neighbour,neighbour[i][5],0,2,offer,i,j))+cos(plaqu(theta,neighbour,i,0,1))+cos(plaqu(theta,neighbour,i,0,2))+cos(plaqu(theta,neighbour,neighbour[i][4],0,1))+cos(plaqu(theta,neighbour,neighbour[i][5],0,2))));
					
					//cout << beta*(cos(offer+theta[neighbour[i][0]][1]-theta[neighbour[i][1]][0]-theta[i][1])+cos(theta[neighbour[i][3]][0]+theta[neighbour[neighbour[i][0]][3]][1]-offer-theta[neighbour[i][3]][1]) - cos(theta[i][0]+theta[neighbour[i][0]][1]-theta[neighbour[i][1]][0]-theta[i][1])-cos(theta[neighbour[i][3]][0]+theta[neighbour[neighbour[i][0]][3]][1]-theta[i][0]-theta[neighbour[i][3]][1])) << endl;	
				}
				else if(j==1)
				{
					//rho = exp(beta*(-cos(theta[i][0]+theta[neighbour[i][0]][1]-theta[neighbour[i][1]][0]-offer)-cos(theta[neighbour[i][2]][0]+offer-theta[neighbour[neighbour[i][1]][2]][0]-theta[neighbour[i][2]][1])+cos(theta[i][0]+theta[neighbour[i][0]][1]-theta[neighbour[i][1]][0]-theta[i][1])+cos(theta[neighbour[i][2]][0]+theta[i][1]-theta[neighbour[neighbour[i][1]][2]][0]-theta[neighbour[i][2]][1])));
				
					rho = exp(-beta*(-cos(plaquchange(theta,neighbour,i,0,1,offer,i,j))-cos(plaquchange(theta,neighbour,i,1,2,offer,i,j))-cos(plaquchange(theta,neighbour,neighbour[i][3],0,1,offer,i,j))-cos(plaquchange(theta,neighbour,neighbour[i][5],1,2,offer,i,j))+cos(plaqu(theta,neighbour,i,0,1))+cos(plaqu(theta,neighbour,i,1,2))+cos(plaqu(theta,neighbour,neighbour[i][3],0,1))+cos(plaqu(theta,neighbour,neighbour[i][5],1,2))));
					
				}
				
				else
				{
					rho = exp(-beta*(-cos(plaquchange(theta,neighbour,i,0,2,offer,i,j))-cos(plaquchange(theta,neighbour,i,1,2,offer,i,j))-cos(plaquchange(theta,neighbour,neighbour[i][3],0,2,offer,i,j))-cos(plaquchange(theta,neighbour,neighbour[i][4],1,2,offer,i,j))+cos(plaqu(theta,neighbour,i,0,2))+cos(plaqu(theta,neighbour,i,1,2))+cos(plaqu(theta,neighbour,neighbour[i][3],0,2))+cos(plaqu(theta,neighbour,neighbour[i][4],1,2))));
					
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
/*
double diracstring(double theta[V][d], int neighbour[V][2*d], int i, int mu, int nu)
{
	double plaquette = plaqu(theta,neighbour,i,)
}
*/
double plaquchange(double theta[V][d], int neighbour[V][2*d],int i, int mu, int nu, double offer, int sitechange, int muchange)
{
	double changedplaqu;
	double save = theta[sitechange][muchange];
	theta[sitechange][muchange] = offer;
	
	changedplaqu = plaqu(theta,neighbour,i,mu,nu);
	
	theta[sitechange][muchange] = save;
	
	return changedplaqu;
}

void calcmonopoledens(double theta[V][d], double M[nmeas], int neighbour[V][2*d], int count)
{
	M[count] = 0;
	int i = 15;
	double n;
	double total_phase;
	//for(int i=0;i<V;i++)
	{
		//M[count] = M[count] + plaqu(theta,neighbour,i,0,1) + plaqu(theta,neighbour,i,2,0) + plaqu(theta,neighbour,i,1,2) + plaqu(theta,neighbour,neighbour[i][0],2,1) + plaqu(theta,neighbour,neighbour[i][1],0,2) + plaqu(theta,neighbour,neighbour[i][2],1,0);	
		//M[count] = M[count] + abs(int(plaqu(theta,neighbour,i,0,1)/(2*PI)) + int(plaqu(theta,neighbour,i,2,0)/(2*PI)) + int(plaqu(theta,neighbour,i,1,2)/(2*PI)) + int(plaqu(theta,neighbour,neighbour[i][0],2,1)/(2*PI)) + int(plaqu(theta,neighbour,neighbour[i][1],0,2)/(2*PI)) + int(plaqu(theta,neighbour,neighbour[i][2],1,0)/(2*PI)));	
		//n = int(plaqu(theta,neighbour,i,0,1)/(2*PI)) + int(plaqu(theta,neighbour,i,2,0)/(2*PI)) + int(plaqu(theta,neighbour,i,1,2)/(2*PI)) + int(plaqu(theta,neighbour,neighbour[i][0],2,1)/(2*PI)) + int(plaqu(theta,neighbour,neighbour[i][1],0,2)/(2*PI)) + int(plaqu(theta,neighbour,neighbour[i][2],1,0)/(2*PI));	
		
		//int(M[imeas]/(2*PI))
		//M[count] = M[count] + plaqu(theta,neighbour,i,0,1)-2*PI*int(plaqu(theta,neighbour,i,0,1)/(2*PI)) + plaqu(theta,neighbour,i,2,0)-2*PI*int(plaqu(theta,neighbour,i,2,0)/(2*PI)) + plaqu(theta,neighbour,i,1,2)-2*PI*int(plaqu(theta,neighbour,i,1,2)/(2*PI)) + plaqu(theta,neighbour,neighbour[i][0],2,1)-2*PI*int(plaqu(theta,neighbour,neighbour[i][0],2,1)/(2*PI)) + plaqu(theta,neighbour,neighbour[i][1],0,2)-2*PI*int(plaqu(theta,neighbour,neighbour[i][1],0,2)/(2*PI)) + plaqu(theta,neighbour,neighbour[i][2],1,0)-2*PI*int(plaqu(theta,neighbour,neighbour[i][2],1,0)/(2*PI));	
		//M[count] = plaqured(theta,neighbour,i,0,1);
		M[count] = M[count] + plaqured(theta,neighbour,i,0,1) + plaqured(theta,neighbour,i,2,0) + plaqured(theta,neighbour,i,1,2) + plaqured(theta,neighbour,neighbour[i][0],2,1) + plaqured(theta,neighbour,neighbour[i][1],0,2)+ plaqured(theta,neighbour,neighbour[i][2],1,0);	
		total_phase = plaqured(theta,neighbour,i,0,1) + plaqured(theta,neighbour,i,2,0) + plaqured(theta,neighbour,i,1,2) + plaqured(theta,neighbour,neighbour[i][0],2,1) + plaqured(theta,neighbour,neighbour[i][1],0,2)+ plaqured(theta,neighbour,neighbour[i][2],1,0);
		
		if (abs(total_phase) < PI)
		{
			n = 0;
		}
          
        else if (total_phase >= PI && total_phase < (3*PI))
        {
        	n = 1;
		}
          
        else if (total_phase >= (3*PI) && total_phase < (5*PI))
        {
        	n = 2;
		}
          
        else if (total_phase >= (5*PI) && total_phase < (7*PI))
        {
        	n = 3;
		}
          
        else if (total_phase <= (-PI) && total_phase > (-3*PI))
        {
        	n = -1;
		}
          
        else if (total_phase <= (-3*PI) && total_phase > (-5*PI))
        {
        	n = -2;
		}
          
        else if (total_phase <= (-5*PI) && total_phase > (-7*PI))
        {
        	n = -3;
		}
          
		
		//M[count] = M[count] + abs(plaqured(theta,neighbour,i,0,1) + plaqured(theta,neighbour,i,0,2) + plaqured(theta,neighbour,i,2,1) + plaqured(theta,neighbour,neighbour[i][0],2,1) + plaqured(theta,neighbour,neighbour[i][1],0,2)+ plaqured(theta,neighbour,neighbour[i][2],0,1));	
		//M[count] = M[count] + abs(plaqured(theta,neighbour,i,0,1) + plaqured(theta,neighbour,i,0,2) + plaqured(theta,neighbour,i,1,2) + plaqured(theta,neighbour,neighbour[i][0],1,2) + plaqured(theta,neighbour,neighbour[i][1],0,2)+ plaqured(theta,neighbour,neighbour[i][2],0,1));	
		//M[count] = M[count] + abs(plaqured(theta,neighbour,i,0,1) + plaqured(theta,neighbour,i,0,2) + plaqured(theta,neighbour,i,1,2) + plaqured(theta,neighbour,neighbour[i][0],1,2) + plaqured(theta,neighbour,neighbour[i][1],0,2)+ plaqured(theta,neighbour,neighbour[i][2],0,1));	

	}
	
	cout << M[count]/(2*PI) << " " << n << endl;
	M[count] = M[count]/(2*PI)/2/V;
}
