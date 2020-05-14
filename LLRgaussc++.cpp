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

double myErfInv2(double);
float my_erfinvf (float a);
float my_logf (float a);

random_device rd;
mt19937 gen(rd());
uniform_real_distribution<> dis(0, 1);


int main()
{
	double epsilon = 0.00001;
	
	double Emax = 40;
	double delta = 0.5;
	int Njacknife = 20;
	int Naverage = 1000;
	
	double x0[int(Emax/delta)];
	double a[int(Emax/delta)];
	double a_i[Njacknife];
	double a_i_new;
	int RobMarchcount;
	double Reweightexpect;
	double y1;
	double y2;
	double y;
	double x;
	double s2[int(Emax/delta)];
	double a_i_del[Njacknife];
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
	for(int Eint=0;Eint<(Emax/delta);Eint++)
	{
		x0[Eint] = Eint*delta;
		
		a[Eint] = 0.0;
		
		//int jcount =0;
		for(int jcount = 0;jcount<Njacknife;jcount++)
		{
			
			if(Eint == 0)
			{
				a_i[jcount] = 0;
				
				a_i_new = -2*epsilon;
			}
			else
			{
				a_i[jcount]=a[Eint-1] - 2*epsilon;
				
				a_i_new = a[Eint-1];
			}
			
			RobMarchcount = 0;
			
			Reweightexpect  = 10.0*epsilon;
			
			while(abs(a_i_new-a_i[jcount]) > epsilon)
			{
				
				a_i[jcount] = a_i_new;
				
				y1 = 0.5*(erf(8*a_i[jcount]+x0[Eint]/16)+1);
				y2 = 0.5*(erf(8*a_i[jcount]+(x0[Eint]+delta)/16)+1);
				
				//cout << "y1: " << y1 << endl;
				//cout << "y2: " << y2 << endl;
 				
				Reweightexpect=0;
				
				for(int i=0;i<Naverage;i++)
				{
					y=(y2-y1)*dis(gen)+y1;
					//x=16*myErfInv2(2*y-1)-128*a_i[jcount];
					x=16*my_erfinvf(2*y-1)-128*a_i[jcount];
					
					//cout << "x: " << x << endl;
					if(abs(x)<10000)
					{
						Reweightexpect = Reweightexpect + x;
					//	myfilemonopoledens << x << " \n";
					}
					else
					{
						//cout << "here " << i << endl;
						i = i-1;
					}
				}
				//cout << "RobMarchcount: " << RobMarchcount << endl;
				Reweightexpect = Reweightexpect/Naverage - x0[Eint] - 0.5*delta;
				
				a_i_new = a_i[jcount] + 12/(delta*delta*(RobMarchcount+1))*Reweightexpect;
				RobMarchcount = RobMarchcount + 1;
			
				
			}
			
			a_i[jcount] = a_i_new;
			a[Eint] = a[Eint] + a_i_new/Njacknife;
			
		}
		
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
		cout << "Energyinterval: " << Eint*delta << endl;
		cout << " " << endl;
		
	}
		myfilemonopoledens.close();
	
	
}



double myErfInv2(double x)
{
   double tt1, tt2, lnx, sgn;
   sgn = (x < 0) ? -1.0f : 1.0f;

   x = (1 - x)*(1 + x);        // x = 1 - x*x;
   lnx = logf(x);

   tt1 = 2/(PI*0.147) + 0.5f * lnx;
   tt2 = 1/(0.147) * lnx;

   return(sgn*sqrtf(-tt1 + sqrtf(tt1*tt1 - tt2)));
}

float my_erfinvf (float a)
{
    float p, r, t;
    t = fmaf (a, 0.0f - a, 1.0f);
    t = my_logf (t);
    if (fabsf(t) > 6.125f) { // maximum ulp error = 2.35793
        p =              3.03697567e-10f; //  0x1.4deb44p-32 
        p = fmaf (p, t,  2.93243101e-8f); //  0x1.f7c9aep-26 
        p = fmaf (p, t,  1.22150334e-6f); //  0x1.47e512p-20 
        p = fmaf (p, t,  2.84108955e-5f); //  0x1.dca7dep-16 
        p = fmaf (p, t,  3.93552968e-4f); //  0x1.9cab92p-12 
        p = fmaf (p, t,  3.02698812e-3f); //  0x1.8cc0dep-9 
        p = fmaf (p, t,  4.83185798e-3f); //  0x1.3ca920p-8 
        p = fmaf (p, t, -2.64646143e-1f); // -0x1.0eff66p-2 
        p = fmaf (p, t,  8.40016484e-1f); //  0x1.ae16a4p-1 
    } else { // maximum ulp error = 2.35456
        p =              5.43877832e-9f;  //  0x1.75c000p-28 
        p = fmaf (p, t,  1.43286059e-7f); //  0x1.33b458p-23 
        p = fmaf (p, t,  1.22775396e-6f); //  0x1.49929cp-20 
        p = fmaf (p, t,  1.12962631e-7f); //  0x1.e52bbap-24 
        p = fmaf (p, t, -5.61531961e-5f); // -0x1.d70c12p-15 
        p = fmaf (p, t, -1.47697705e-4f); // -0x1.35be9ap-13 
        p = fmaf (p, t,  2.31468701e-3f); //  0x1.2f6402p-9 
        p = fmaf (p, t,  1.15392562e-2f); //  0x1.7a1e4cp-7 
        p = fmaf (p, t, -2.32015476e-1f); // -0x1.db2aeep-3 
        p = fmaf (p, t,  8.86226892e-1f); //  0x1.c5bf88p-1 
    }
    r = a * p;
    return r;
}


float my_logf (float a)
{
    float i, m, r, s, t;
    int e;

    m = frexpf (a, &e);
    if (m < 0.666666667f) { // 0x1.555556p-1
        m = m + m;
        e = e - 1;
    }
    i = (float)e;
    /* m in [2/3, 4/3] */
    m = m - 1.0f;
    s = m * m;
    /* Compute log1p(m) for m in [-1/3, 1/3] */
    r =             -0.130310059f;  // -0x1.0ae000p-3
    t =              0.140869141f;  //  0x1.208000p-3
    r = fmaf (r, s, -0.121484190f); // -0x1.f19968p-4
    t = fmaf (t, s,  0.139814854f); //  0x1.1e5740p-3
    r = fmaf (r, s, -0.166846052f); // -0x1.55b362p-3
    t = fmaf (t, s,  0.200120345f); //  0x1.99d8b2p-3
    r = fmaf (r, s, -0.249996200f); // -0x1.fffe02p-3
    r = fmaf (t, m, r);
    r = fmaf (r, m,  0.333331972f); //  0x1.5554fap-2
    r = fmaf (r, m, -0.500000000f); // -0x1.000000p-1
    r = fmaf (r, s, m);
    r = fmaf (i,  0.693147182f, r); //  0x1.62e430p-1 // log(2)
    if (!((a > 0.0f) && (a <= 3.40282346e+38f))) { // 0x1.fffffep+127
        r = a + a;  // silence NaNs if necessary
        if (a  < 0.0f) r = ( 0.0f / 0.0f); //  NaN
        if (a == 0.0f) r = (-1.0f / 0.0f); // -Inf
    }
    return r;
}
