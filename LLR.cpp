// -----------------------------------------------------------------
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

random_device rd;
mt19937 gen(rd());
uniform_real_distribution<> dis(0, 1);
uniform_real_distribution<> dist(-PI, PI);
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Double precision CPU time in seconds
double dclock() {
  long fine;
  fine = clock();
  return (((double)fine) / CLOCKS_PER_SEC);
}
// -----------------------------------------------------------------


// -----------------------------------------------------------------
// Initialization routines
// Initialize neighbours with periodic BCs
void neibinit(int neib[][2*d], int n) {
  int i, n1p, n1m, n2p, n2m, n3p, n3m, n4p, n4m;
  for (int n1=0;n1<n;n1++) {
    n1p = n1+1;
    n1m = n1-1;
    if (n1p == n)
      n1p = 0;
    if (n1m ==-1)
      n1m = n-1;
    for (int n2=0;n2<n;n2++) {
      n2p = n2+1;
      n2m = n2-1;
      if (n2p == n)
        n2p = 0;
      if (n2m ==-1)
        n2m = n-1;
      for(int n3=0;n3<n;n3++) {
        n3p = n3+1;
        n3m = n3-1;
        if  (n3p == n)
          n3p = 0;
        if  (n3m ==-1)
          n3m = n-1;
        for(int n4=0;n4<n;n4++) {
          n4p = n4+1;
          n4m = n4-1;
          if  (n4p == n)
            n4p = 0;
          if  (n4m ==-1)
            n4m = n-1;

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

void filllinks(double theta[][d],int start, int V) {
  if(start == 1) {
    for(int i=0;i<V;i++)
    {
      for(int j=0;j<d;j++)
        theta[i][j]=dist(rng);
    }
  }

  if(start == 2) {
    for(int i=0;i<V;i++)
    {
      for(int j=0;j<d;j++)
        theta[i][j]=-PI;
    }
  }

  if(start == 3) {
    for(int i=0;i<V;i++)
    {
      for(int j=0;j<d;j++)
        theta[i][j]=0.0;
    }
  }
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Three basic plaquette and action routines
double plaqu(double theta[][d], int neighbour[][2*d], int i, int mu, int nu,int V) {
  return (theta[i][mu] + theta[neighbour[i][mu]][nu] - theta[neighbour[i][nu]][mu] - theta[i][nu]);
}

double plaquchange(double theta[][d], int neighbour[][2*d],int i, int mu, int nu, double offer, int sitechange, int muchange,int V) {
  double changedplaqu;
  double save = theta[sitechange][muchange];
  theta[sitechange][muchange] = offer;

  changedplaqu = plaqu(theta,neighbour,i,mu,nu,V);

  theta[sitechange][muchange] = save;

  return changedplaqu;
}

double calcaction(double theta[][d], int neighbour[][2*d], double beta, int V) {
  double S = 0;
  for(int i = 0;i<V;i++) {
    for(int j=0;j<(d-1);j++) {
      for(int k =j+1;k<d;k++)
        S += cos(plaqu(theta,neighbour,i,j,k,V));
    }
  }
  S *= beta;
  return S;
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Routine to find the targeted energy interval
void findEint(double theta[][d],int neighbour[][2*d], double beta, double Eint, double delta,int V) {
  bool Efound = false;
  double td;

  if(Efound == false) {
    for(int k=-15;k<0;k++) {
      for(int i=0;i<V;i=i+2) {
        for(int j=0;j<d;j++) {
          theta[i][j] = -k*PI/16;
          //cout << "Action/V: " << calcaction(theta,neighbour,beta,V)/V << endl;
          td = calcaction(theta,neighbour,beta,V);
          if(td>=(Eint) && td<=((Eint+delta))) {
            i=V;
            j=d;
            k=0;
            Efound = true;
          }
        }
      }
    }
  }
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Four convenience routines to allow some code re-use
double Sdiff0(double theta[][d], int neighbour[][2*d], double offer, int i, int V) {
  double td = -cos(plaquchange(theta,neighbour,i,0,1,offer,i,0,V));
  td -= cos(plaquchange(theta,neighbour,i,0,2,offer,i,0,V));
  td -= cos(plaquchange(theta,neighbour,i,0,3,offer,i,0,V));
  td -= cos(plaquchange(theta,neighbour,neighbour[i][5],0,1,offer,i,0,V));
  td -= cos(plaquchange(theta,neighbour,neighbour[i][6],0,2,offer,i,0,V));
  td -= cos(plaquchange(theta,neighbour,neighbour[i][7],0,3,offer,i,0,V));
  td += cos(plaqu(theta,neighbour,i,0,1,V));
  td += cos(plaqu(theta,neighbour,i,0,2,V));
  td += cos(plaqu(theta,neighbour,i,0,3,V));
  td += cos(plaqu(theta,neighbour,neighbour[i][5],0,1,V));
  td += cos(plaqu(theta,neighbour,neighbour[i][6],0,2,V));
  td += cos(plaqu(theta,neighbour,neighbour[i][7],0,3,V));
  return td;
}
double Sdiff1(double theta[][d], int neighbour[][2*d], double offer, int i, int V) {
  double td = -cos(plaquchange(theta,neighbour,i,0,1,offer,i,1,V));
  td -= cos(plaquchange(theta,neighbour,i,1,2,offer,i,1,V));
  td -= cos(plaquchange(theta,neighbour,i,1,3,offer,i,1,V));
  td -= cos(plaquchange(theta,neighbour,neighbour[i][4],0,1,offer,i,1,V));
  td -= cos(plaquchange(theta,neighbour,neighbour[i][6],1,2,offer,i,1,V));
  td -= cos(plaquchange(theta,neighbour,neighbour[i][7],1,3,offer,i,1,V));
  td += cos(plaqu(theta,neighbour,i,0,1,V));
  td += cos(plaqu(theta,neighbour,i,1,2,V));
  td += cos(plaqu(theta,neighbour,i,1,3,V));
  td += cos(plaqu(theta,neighbour,neighbour[i][4],0,1,V));
  td += cos(plaqu(theta,neighbour,neighbour[i][6],1,2,V));
  td += cos(plaqu(theta,neighbour,neighbour[i][7],1,3,V));
  return td;
}
double Sdiff2(double theta[][d], int neighbour[][2*d], double offer, int i, int V) {
  double td = -cos(plaquchange(theta,neighbour,i,0,2,offer,i,2,V));
  td -= cos(plaquchange(theta,neighbour,i,1,2,offer,i,2,V));
  td -= cos(plaquchange(theta,neighbour,i,2,3,offer,i,2,V));
  td -= cos(plaquchange(theta,neighbour,neighbour[i][4],0,2,offer,i,2,V));
  td -= cos(plaquchange(theta,neighbour,neighbour[i][5],1,2,offer,i,2,V));
  td -= cos(plaquchange(theta,neighbour,neighbour[i][7],2,3,offer,i,2,V));
  td += cos(plaqu(theta,neighbour,i,0,2,V));
  td += cos(plaqu(theta,neighbour,i,1,2,V));
  td += cos(plaqu(theta,neighbour,i,2,3,V));
  td += cos(plaqu(theta,neighbour,neighbour[i][4],0,2,V));
  td += cos(plaqu(theta,neighbour,neighbour[i][5],1,2,V));
  td += cos(plaqu(theta,neighbour,neighbour[i][7],2,3,V));
  return td;
}
double Sdiff3(double theta[][d], int neighbour[][2*d], double offer, int i, int V) {
  double td = -cos(plaquchange(theta,neighbour,i,0,3,offer,i,3,V));
  td -= cos(plaquchange(theta,neighbour,i,1,3,offer,i,3,V));
  td -= cos(plaquchange(theta,neighbour,i,2,3,offer,i,3,V));
  td -= cos(plaquchange(theta,neighbour,neighbour[i][4],0,3,offer,i,3,V));
  td -= cos(plaquchange(theta,neighbour,neighbour[i][5],1,3,offer,i,3,V));
  td -= cos(plaquchange(theta,neighbour,neighbour[i][6],2,3,offer,i,3,V));
  td += cos(plaqu(theta,neighbour,i,0,3,V));
  td += cos(plaqu(theta,neighbour,i,1,3,V));
  td += cos(plaqu(theta,neighbour,i,2,3,V));
  td += cos(plaqu(theta,neighbour,neighbour[i][4],0,3,V));
  td += cos(plaqu(theta,neighbour,neighbour[i][5],1,3,V));
  td += cos(plaqu(theta,neighbour,neighbour[i][6],2,3,V));
  return td;
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// MRT updates with and without constraints
void metropolisupdate(double theta[][d],double a,double beta, int neighbour[][2*d], int nsweeps, int V) {
  double offer; //new offered linkvar
  double rho, r; //Probabilities for acceptance/rejection
  double td;

  for(int n=0;n<nsweeps;n++) {
    for(int i=0;i<V;i++) {
      for(int j=0;j<d;j++) {
        offer = dist(rng);

        if(j==0)
          td = Sdiff0(theta,neighbour,offer,i,V);
        else if(j==1)
          td = Sdiff1(theta,neighbour,offer,i,V);
        else if(j==2)
          td = Sdiff2(theta,neighbour,offer,i,V);
        else
          td = Sdiff3(theta,neighbour,offer,i,V);

        rho = exp(beta * a * td);
        r = dis(rng);
        if(r<=rho)
          theta[i][j]=offer;
      }
    }
  }
}

void metropolisupdateconst(double theta[][d],double a,double beta, int neighbour[][2*d], int nsweeps,double E,double delta, int V) {
  int counter1 = 0, counter2 = 0, counter3 = 0;
  double offer; //new offered linkvar
  double rho, r; //Probabilities for acceptance/rejection
  double td;

  double ECALC = calcaction(theta,neighbour,beta,V);
  double ECALCtry = ECALC;

  for(int n=0;n<nsweeps;n++) {
    for(int i=0;i<V;i++) {
      for(int j=0;j<d;j++) {
        offer = dist(rng);

        if(j==0)
          td = Sdiff0(theta,neighbour,offer,i,V);
        else if(j==1)
          td = Sdiff1(theta,neighbour,offer,i,V);
        else if(j==2)
          td = Sdiff2(theta,neighbour,offer,i,V);
        else
          td = Sdiff3(theta,neighbour,offer,i,V);

        rho = exp(beta * a * td);
        r = dis(rng);

        counter2++;

        if(r<=rho)
        {
          ECALCtry = ECALC - beta * td;
          if(ECALCtry<(E+delta) && ECALCtry>=E)
          {
            //cout << "offer accepted " << i <<   endl;
            counter1++;
            theta[i][j]=offer;
            ECALC = ECALCtry;
          }
          else
          {
            j=j-1;
            counter3++;
            //cout << offer <<  "    " << calcaction(theta,neighbour,beta) << "i: " << i << " j: " << j+1 << endl;
          }

          //cout << "ECALC: " << ECALC/V << " calcaction: " << calcaction(theta,neighbour,beta,V) << endl;
        }
      }
    }
  }
  /*
     if(counter1<10)
     {
     cout << "counter1: "<< counter1 << " counter2: " << counter2 << endl;
  //cout << "ratio1: " << double(counter1)/double(counter2) << " ratio2: " << double(counter1)/double(counter3) << endl;
  }
  */
  //cout << counter1 << " " << counter2 << " " << counter3 << endl;
}

void metropolisupdateEfind(double theta[][d],double beta,int neighbour[][2*d],int i, int j,int V) {
  double offer; //new offered linkvar
  double rho, r; //Probabilities for acceptance/rejection
  double td;

  offer = dist(rng);

  if(j==0)
    td = Sdiff0(theta,neighbour,offer,i,V);
  else if(j==1)
    td = Sdiff1(theta,neighbour,offer,i,V);
  else if(j==2)
    td = Sdiff2(theta,neighbour,offer,i,V);
  else
    td = Sdiff3(theta,neighbour,offer,i,V);

  rho = exp(-beta * td);
  r = dis(rng);
  if(r<=rho)
    theta[i][j]=offer;
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Finally main routine
int main() {
  double dtime = -dclock();
  string line;
  string saver[12];
  int savercounter = 0;
  ifstream myfileinput("input.txt");
  if (myfileinput.is_open())
  {
    while (getline(myfileinput,line))
    {
      saver[savercounter] = line;
      cout << line << '\n';
      savercounter = savercounter + 1;
    }
    myfileinput.close();

  }

  else
  {
    cout << "Unable to open file";
    return 0;
  }

  string nl = saver[0].substr(saver[0].find("=")+1);
  stringstream nstring(nl);
  int n = 0;
  nstring >> n;
  cout << n << endl;

  int V=pow(n,d);
  cout << V << endl;

  string nequil = saver[1].substr(saver[1].find("=")+1);
  stringstream nequistring(nequil);
  int N_TH = 0;
  nequistring >> N_TH;
  cout << N_TH << endl;

  string nmeasl = saver[3].substr(saver[3].find("=")+1);
  stringstream nmeasstring(nmeasl);
  int N_SW = 0;
  nmeasstring >> N_SW;
  cout << N_SW << endl;

  string nskipl = saver[2].substr(saver[2].find("=")+1);
  stringstream nskipstring(nskipl);
  int nskip = 0;
  nskipstring >> nskip;
  cout << nskip << endl;

  string betal = saver[4].substr(saver[4].find("=")+1);
  stringstream betastring(betal);
  double beta = 0;
  betastring >> beta;
  cout << beta << endl;

  string Eminl = saver[5].substr(saver[5].find("=")+1);
  stringstream Eminstring(Eminl);
  double Emin = 0;
  Eminstring >> Emin;
  cout << Emin << endl;
  Emin = Emin*V;

  // TODO: Could switch to Emin and # of intervals
  string Emaxl = saver[6].substr(saver[6].find("=")+1);
  stringstream Emaxstring(Emaxl);
  double Emax = 0;
  Emaxstring >> Emax;
  cout << Emax << endl;
  Emax = Emax*V;

  string deltal = saver[7].substr(saver[7].find("=")+1);
  stringstream deltastring(deltal);
  double delta = 0;
  deltastring >> delta;
  cout << delta << endl;
  delta=delta*V;

  string epsilonl = saver[8].substr(saver[8].find("=")+1);
  stringstream epsilonstring(epsilonl);
  double epsilon = 0;
  epsilonstring >> epsilon;
  cout << epsilon << endl;

  string Njacknife1 = saver[9].substr(saver[9].find("=")+1);
  stringstream Njacknifestring(Njacknife1);
  int Njacknife = 0;
  Njacknifestring >> Njacknife;
  cout << Njacknife << endl;

  string ait1 = saver[10].substr(saver[10].find("=")+1);
  stringstream aitstring(ait1);
  int ait=0;
  aitstring >> ait;
  cout << ait << endl;

  string seed1 = saver[11].substr(saver[11].find("=")+1);
  stringstream seedstring(seed1);
  int seed=0;
  seedstring >> seed;
  cout << seed << endl;
  rng.seed(seed);

  int neighbour[V][2*d];
  neibinit(neighbour,n);
  double theta[V][d];
  bool Einterval = false;

  // Arrays only need int((Emax-Emin)/delta+1) elements
  int Nmax = int(Emax/delta);
  double *x0 = (double*)malloc(Nmax * sizeof(double));  //lower end of energy interval
  double *a = (double*)malloc(Nmax * sizeof(double));
  double *a_i = (double*)malloc(Njacknife * sizeof(double));
  double a_i_new;
  double Reweightexpect; //Reweighted expectationvalue of the energy
  double *measurement = (double*)malloc(N_SW * sizeof(double));
  double varianz;

  double *s2 = (double*)malloc(Nmax * sizeof(double)); //variance
  //stuff for calculation of density from a
  double A[Nmax][Njacknife];
  double Astar;
  double *var = (double*)malloc(Nmax * sizeof(double));
  double T;
  double *rhojack = (double*)malloc(Njacknife * sizeof(double));
  double *rho = (double*)malloc(Nmax * sizeof(double));
  double *errorrho = (double*)malloc(Nmax * sizeof(double));
  double td;

  ofstream myfilerhovalues;
  myfilerhovalues.open("rhovalues.txt");
  ofstream myfileoutput;
  myfileoutput.open("output.txt");
  ofstream myfileoutputinnerloop;
  myfileoutputinnerloop.open("outputinnerloop.txt");
  ofstream myfileavalues;
  myfileavalues.open("avalues.txt");
  ofstream myfileavar;
  myfileavar.open("avar.txt");
  ofstream myfilerhoerror;
  myfilerhoerror.open("rhoerror.txt");
  ofstream myfileajack;
  myfileajack.open("ajack.txt");

  for(int Eint=0;Eint<=((Emax-Emin)/delta);Eint++) //scans through the energy intervals
  {
    x0[Eint] = Eint*delta+Emin;
    a[Eint] = 0.0;//-1.0018;

    //int jcount =0;
    for(int jcount = 0;jcount<Njacknife;jcount++) //generates Njacknife values of a for Jacknife method to get error and average
    {
      //set starting value for a
      if(Eint == 0)
      {
        a_i[jcount] = -0.5+2*epsilon;//-1.0018+2*epsilon;

        a_i_new = -0.5;//-1.0018;
      }
      else
      {
        a_i[jcount]=a[Eint-1] - 2*epsilon;

        a_i_new = a[Eint-1];
      }

      filllinks(theta,2,V);
      Einterval = false;

      //while(abs(a_i_new-a_i[jcount]) > epsilon) //only stop rob march if a doesnt change much anymore
      for(int RMcount=0; RMcount<ait; RMcount++)
      {
        cout << "difference " << a_i_new-a_i[jcount] << endl;
        myfileoutputinnerloop << "difference " << a_i_new-a_i[jcount] << endl;

        a_i[jcount] = a_i_new;

        //cout << Einterval << endl;

        //finds the right energy intervall starting from the maximum energy (all link variables set to 0)
        // Only runs for RMcount=0 and could be moved outside RMcount loop
        if(Einterval == false)
        {
          findEint(theta,neighbour,beta,x0[Eint],delta,V);
          Einterval = true;
          cout << "found interval" << endl;
          metropolisupdateconst(theta,1.0,beta,neighbour,2*N_TH,x0[Eint],delta,V);
        }

        cout << "action reached: " << calcaction(theta,neighbour,beta,V)/V << endl;
        //cout << "a " << a_i[jcount] << endl;
        metropolisupdateconst(theta,a_i[jcount],beta,neighbour,N_TH,x0[Eint],delta,V);

        td = calcaction(theta,neighbour,beta,V);
        cout << "RMcount " <<  RMcount << " Energy " << td/V << endl;
        myfileoutputinnerloop << "RMcount " <<  RMcount << " Energy " << td/V << endl;

        //calculate the reweighted expectation value of the energy
        Reweightexpect=0;
        for(int k = 0;k<N_SW;k++)
        {
          measurement[k] = calcaction(theta,neighbour,beta,V);
          Reweightexpect += measurement[k];
          td = Reweightexpect / (k + 1) - x0[Eint] - 0.5*delta;
          cout << "measurement[" << k << "] = " << measurement[k] << " --> Reweightexpect = " << td;
          cout << " a = " << a_i[jcount] + 12.0 * td / (delta*delta*(RMcount+1)) << endl;
          metropolisupdateconst(theta,a_i[jcount],beta,neighbour,nskip,x0[Eint],delta,V);
        }
        Reweightexpect /= N_SW;

        varianz=0;
        for(int k=0;k<N_SW;k++)
          varianz += pow(measurement[k]-Reweightexpect,2);
        varianz /= N_SW;

        Reweightexpect = Reweightexpect - x0[Eint] - 0.5*delta;

        //calculate the new a
        a_i_new = a_i[jcount] + 12/(delta*delta*(RMcount+1))*Reweightexpect;
        cout << "Reweight expect: " << Reweightexpect << endl;
        cout << "var: " << varianz << endl;
        cout << "a: " << a_i_new << endl;
        myfileoutputinnerloop <<  "a: " << a_i_new << endl;
        myfileoutputinnerloop << "Reweight expect: " << Reweightexpect << endl;
        myfileoutputinnerloop <<  "var: " << varianz << endl;
      }

      a_i[jcount] = a_i_new;
      myfileajack << a_i[jcount] << endl;
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

    errorrho[Eint]=0;
    for(int i=0;i<Njacknife;i++)
    {
      errorrho[Eint]= errorrho[Eint] + pow(rho[Eint]-rhojack[i],2);
    }
    errorrho[Eint] = sqrt(Njacknife-1)/Njacknife*errorrho[Eint];

    cout << Eint << endl;
    cout << "rho: " << rho[Eint] << endl;
    myfilerhovalues << rho[Eint] << " \n";
    myfileavalues << a[Eint] << " \n";
    myfileavar << s2[Eint] << " \n";
    myfilerhoerror << errorrho[Eint] << endl;
    td = calcaction(theta,neighbour,beta,V);
    myfileoutput << "Reweightexpect: " << Reweightexpect << "a: " << a[Eint] << "Energyinterval: " << Eint*delta/6/V <<"Energy: " << td/6/V << " \n";
    cout << "Reweightexpect: " << Reweightexpect << endl;
    cout << "a: " << a[Eint] << endl;
    cout << "Energyinterval: " << Eint*delta/6/V << endl;
    cout << "Energy: " << td/6/V << endl;
    cout << " " << endl;

  }
  myfilerhovalues.close();
  myfileavalues.close();
  myfileavar.close();
  myfilerhoerror.close();
  myfileoutput.close();
  myfileoutputinnerloop.close();
  myfileajack.close();
  dtime += dclock();
  cout << "Total time: " << dtime << " sec" << endl;
}
// -----------------------------------------------------------------
