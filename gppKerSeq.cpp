/*
Sequential GPP code that uses std:complex<double> data type. 
*/
#include <iostream>
#include <cstdlib>
#include <memory>
#include <iomanip>
#include <cmath>
#include <complex>
#include <sys/time.h>

using namespace std;

#define nstart 0
#define nend 3

//Outputs are ssxa and scha, rest of the passed parameters are the inputs
void ssxt_scht_solver(double wxt, int igp, int my_igp, int ig, std::complex<double> wtilde, std::complex<double> wtilde2, std::complex<double> Omega2, std::complex<double> matngmatmgp, std::complex<double> matngpmatmg, std::complex<double> mygpvar1, std::complex<double> mygpvar2, std::complex<double>& ssxa, std::complex<double>& scha, std::complex<double> I_eps_array_igp_myIgp)
{
    std::complex<double> expr0( 0.0 , 0.0);
    double delw2, scha_mult, ssxcutoff;
    double to1      = 1e-6;
    double sexcut   = 4.0;
    double gamma    = 0.5;
    double limitone = 1.0/(to1*4.0);
    double limittwo = pow(0.5,2);
    std::complex<double> sch(0.00, 0.00);
    std::complex<double> ssx(0.00, 0.00);

    std::complex<double> wdiff = wxt - wtilde;

    std::complex<double> cden = wdiff;
    double rden               = 1/real(cden * conj(cden));
    std::complex<double> delw = wtilde * conj(cden) * rden;
    double delwr              = real(delw * conj(delw));
    double wdiffr             = real(wdiff * conj(wdiff));

    if((wdiffr > limittwo) && (delwr < limitone)){
      sch = delw * I_eps_array_igp_myIgp;
      cden = pow(wxt,2);
      rden = real(cden * conj(cden));
      rden = 1.00 / rden;
      ssx = Omega2 * conj(cden) * rden;
    } else if (delwr > to1) {
      sch = expr0;
      cden = (double) 4.00 * wtilde2 * (delw + (double)0.50);
      rden = real(cden * conj(cden));
      rden = 1.00/rden;
        ssx = -Omega2 * conj(cden) * rden * delw;
    } else {
      sch = expr0;
      ssx = expr0;
    }

    ssxcutoff = sexcut*abs(I_eps_array_igp_myIgp);
    if((abs(ssx) > ssxcutoff) && (abs(wxt) < 0.00))
      ssx = 0.00;

    ssxa = matngmatmgp*ssx;
    scha = matngmatmgp*sch;
}

//This function writes its results to achstemp, rest of the parameters are its inputs.
void reduce_achstemp(int n1, int* inv_igp_index, int ncouls, std::complex<double> *aqsmtemp, std::complex<double> *aqsntemp, std::complex<double> *I_eps_array, std::complex<double>& achstemp, int ngpown, double* vcoul){
  double to1 = 1e-6;
  for(int my_igp = 0; my_igp< ngpown; my_igp++){
    std::complex<double> schstemp(0.0, 0.0);
    std::complex<double> schs(0.0, 0.0);
    std::complex<double> matngmatmgp(0.0, 0.0);
    std::complex<double> matngpmatmg(0.0, 0.0);
    std::complex<double> mygpvar1(0.00, 0.00);
    std::complex<double> mygpvar2(0.00, 0.00);
    int igp = inv_igp_index[my_igp];
    if(igp >= ncouls)
      igp = ncouls-1;

    if(!(igp > ncouls || igp < 0)){
      mygpvar1    = std::conj(aqsmtemp[n1*ncouls+igp]);
      mygpvar2    = aqsntemp[n1*ncouls+igp];
      schs        = -I_eps_array[my_igp*ncouls+igp];
      matngmatmgp = aqsntemp[n1*ncouls+igp] * mygpvar1;

      if(abs(schs) > to1)
          schstemp = schstemp + matngmatmgp * schs;
    } else {
        for(int ig=1; ig<ncouls; ++ig)
            schstemp = schstemp - aqsntemp[n1*ncouls+igp] * I_eps_array[my_igp*ncouls+ig] * mygpvar1;
    }
    achstemp += schstemp * vcoul[igp] *(double) 0.5;
  }
}



//Performs the calculation for the first nvband iterations.
//Outputs are ssxt and scht, rest of the passed parameters are the inputs
void flagOCC_solver(double wxt, std::complex<double> *wtilde_array, int my_igp, int n1, std::complex<double> *aqsmtemp, std::complex<double> *aqsntemp, std::complex<double> *I_eps_array, std::complex<double> &ssxt, std::complex<double> &scht, int ncouls, int igp, std::complex<double> *ssxa, std::complex<double>* scha)
{
    std::complex<double> matngmatmgp = std::complex<double>(0.0, 0.0);
    std::complex<double> matngpmatmg = std::complex<double>(0.0, 0.0);
    for(int ig=0; ig<ncouls; ++ig){
      std::complex<double> wtilde = wtilde_array[my_igp*ncouls+ig];
      std::complex<double> wtilde2 = std::pow(wtilde,2);
      std::complex<double> Omega2 = wtilde2*I_eps_array[my_igp*ncouls+ig];
      std::complex<double> mygpvar1 = std::conj(aqsmtemp[n1*ncouls+igp]);
      std::complex<double> mygpvar2 = aqsmtemp[n1*ncouls+igp];
      std::complex<double> matngmatmgp = aqsntemp[n1*ncouls+ig] * mygpvar1;
      if(ig != igp) matngpmatmg = std::conj(aqsmtemp[n1*ncouls+ig]) * mygpvar2;

      ssxt_scht_solver(wxt, igp, my_igp, ig, wtilde, wtilde2, Omega2, matngmatmgp, matngpmatmg, mygpvar1, mygpvar2, ssxa[ig], scha[ig], I_eps_array[my_igp*ncouls+ig]); 
      ssxt += ssxa[ig];
      scht += scha[ig];
    }
}

//Outputs is scht, rest of the passed parameters are the inputs
void noflagOCC_solver(double wxt, std::complex<double> *wtilde_array, int my_igp, int n1, std::complex<double> *aqsmtemp, std::complex<double> *aqsntemp, std::complex<double> *I_eps_array, std::complex<double> &ssxt, std::complex<double> &scht, int ncouls, int igp, std::complex<double> *scha){
  double to1                    = 1e-6;
  double sexcut                 = 4.0;
  double gamma                  = 0.5;
  double limitone               = 1.0/(to1*4.0);
  double limittwo               = pow(0.5,2);
  std::complex<double> mygpvar1 = std::conj(aqsmtemp[n1*ncouls+igp]);
  std::complex<double> scht_loc(0.00, 0.00);
  
  for(int ig = 0; ig<ncouls; ++ig){
    std::complex<double> wdiff = wxt - wtilde_array[my_igp*ncouls+ig];
    double wdiffr = real(wdiff * conj(wdiff));
    double rden = 1/wdiffr;

    std::complex<double> delw = wtilde_array[my_igp*ncouls+ig] * conj(wdiff) *rden; //*rden
    double delwr = real(delw * conj(delw));

    scht_loc += mygpvar1 * aqsntemp[n1*ncouls+ig] * delw * I_eps_array[my_igp*ncouls+ig] ;
  }

  scht = scht_loc;
}

int main(int argc, char** argv){
  //The input to the executable needs 4 arguments.
  if (argc != 5)
  {
      std::cout << "The correct form of input is : " << endl;
      std::cout << " ./a.out <number_bands> <number_valence_bands> <number_plane_waves> <matrix_divider> " << endl;
      exit (0);
  }

//Input parameters stored in these variables.
  const int number_bands    = atoi(argv[1]);
  const int nvband          = atoi(argv[2]);
  const int ncouls          = atoi(argv[3]);
  const int nodes_per_group = atoi(argv[4]);

//Constants that will be used later
  const int npes        = 1; 
  const int ngpown      = ncouls / (nodes_per_group * npes); 
  const double e_lk     = 10;
  const double dw       = 1;
  const double to1      = 1e-6;
  const double gamma    = 0.5;
  const double sexcut   = 4.0;
  const double limitone = 1.0/(to1*4.0);
  const double limittwo = pow(0.5,2);
  const double e_n1kq   = 6.0; 
  const double occ      = 1.0;

  //Printing out the params passed.
  std::cout << "**************************** Sequential GPP code ************************* " << std::endl;
  std::cout << "number_bands = " << number_bands \
      << "\t nvband = " << nvband \
      << "\t ncouls = " << ncouls \
      << "\t nodes_per_group  = " << nodes_per_group \
      << "\t ngpown = " << ngpown \
      << "\t nend = " << nend \
      << "\t nstart = " << nstart \
      << "\t gamma = " << gamma \
      << "\t sexcut = " << sexcut \
      << "\t limitone = " << limitone \
      << "\t limittwo = " << limittwo << endl;

  // Memory allocation of input data structures.
  // Two dimensional arrays from theory have been initialized as a single dimension in m*n format for performance.
  std::complex<double> *acht_n1_loc = new std::complex<double> [number_bands];
  std::complex<double> *aqsmtemp = new std::complex<double> [number_bands*ncouls];
  std::complex<double> *aqsntemp = new std::complex<double> [number_bands*ncouls];
  std::complex<double> *I_eps_array = new std::complex<double> [ngpown*ncouls];
  std::complex<double> *wtilde_array = new std::complex<double> [ngpown*ncouls];
  int *inv_igp_index = new int[ngpown];
  double *vcoul = new double[ncouls];
  double wx_array[nend-nstart];

  //arrays that will be later used to store the output results
  std::complex<double> achtemp[nend-nstart]; 
  std::complex<double> asxtemp[nend-nstart];

  std::complex<double> achstemp = std::complex<double>(0.0, 0.0);

  //Data structures that store intermediete results
  std::complex<double> ssx_array[nend-nstart], \
      sch_array[nend-nstart], \
      scht, ssxt;
  std::complex<double> *ssxa = new std::complex<double> [ncouls];
  std::complex<double> *scha = new std::complex<double> [ncouls];

  //Printing the size of each of the input data structures.
  cout << "Size of wtilde_array = " << (ncouls*ngpown*2.0*8) / pow(1024,2) << " Mbytes" << endl;
  cout << "Size of aqsntemp = " << (ncouls*number_bands*2.0*8) / pow(1024,2) << " Mbytes" << endl;
  cout << "Size of I_eps_array array = " << (ncouls*ngpown*2.0*8) / pow(1024,2) << " Mbytes" << endl;

  //Some expressions declared to be used later in the initialization.
  std::complex<double> expr0( 0.0 , 0.0);
  std::complex<double> expr( 0.5 , 0.5);

  //Initializing the data structures
  for(int i=0; i<number_bands; i++)
  for(int j=0; j<ncouls; j++){
    aqsntemp[i*ncouls+j] = ((double)(i+j))*expr;
    aqsmtemp[i*ncouls+j] = ((double)(i+j))*expr;
  }


  for(int i=0; i<ngpown; i++){
    for(int j=0; j<ncouls; j++){
      I_eps_array[i*ncouls+j]  = ((double)(i+j))*expr;
      wtilde_array[i*ncouls+j] = ((double)(i+j))*expr;
    }

    inv_igp_index[i] = (i+1) * ncouls / ngpown;
  }

  for(int i=0; i<ncouls; i++){
    ssxa[i]  = expr0;
    scha[i]  = expr0;
    vcoul[i] = 1.0*i;
  }


  for(int iw=nstart; iw<nend; ++iw){
    achtemp[iw]                              = expr0;
    asxtemp[iw]                              = expr0;
    wx_array[iw]                             = e_lk - e_n1kq + dw*((iw+1)-2);
    if(abs(wx_array[iw]) < to1) wx_array[iw] = to1;
  }

  //Start the timer before the work begins.
  //Start the timer before the work begins.
  timeval startTimer, endTimer;
  gettimeofday(&startTimer, NULL);

  //The main work starts here
  for(int n1 = 0; n1<number_bands; ++n1) {
    bool flag_occ = n1 < nvband;
    reduce_achstemp(n1, inv_igp_index, ncouls, aqsmtemp, aqsntemp, I_eps_array, achstemp, ngpown, vcoul);

    for(int my_igp=0; my_igp<ngpown; ++my_igp){
      int igp = inv_igp_index[my_igp];
      if(igp >= ncouls)
          igp = ncouls-1;

      //Reinitialize the intermediete variables to initial values.
      for(int i=nstart; i<nend; i++){
        ssx_array[i] = expr0;
        sch_array[i] = expr0;
      }

      if(flag_occ){ //iterations from 0-nvband from the outer loop enter this conditional
        for(int iw=nstart; iw<nend; iw++){
          scht       = ssxt = expr0;
          double wxt = wx_array[iw];
          flagOCC_solver(wxt, wtilde_array, my_igp, n1, aqsmtemp, aqsntemp, I_eps_array, ssxt, scht, ncouls, igp, ssxa, scha);
          ssx_array[iw] += ssxt;
          sch_array[iw] +=(double) 0.5*scht;
          
          asxtemp[iw] += ssx_array[iw] * occ * vcoul[igp];//Store output of the first nvband iterations.
        }
      } else { //rest of the iterations i.e., nvband-number_bands enter this else-loop
        for(int iw=nstart; iw<nend; ++iw){
          scht       = ssxt = expr0;
          double wxt = wx_array[iw];
          noflagOCC_solver(wxt, wtilde_array, my_igp, n1, aqsmtemp, aqsntemp, I_eps_array, ssxt, scht, ncouls, igp, scha);

          sch_array[iw] +=(double) 0.5*scht;
        }
      }

      for(int iw=nstart; iw<nend; ++iw)
        achtemp[iw] += sch_array[iw] * vcoul[igp];//Store final output here

      acht_n1_loc[n1] += sch_array[2] * vcoul[igp];
    }
  }

  //Time Taken
  gettimeofday(&endTimer, NULL);
  double elapsedTimer = (endTimer.tv_sec - startTimer.tv_sec) +1e-6*(endTimer.tv_usec - startTimer.tv_usec);

  for(int iw=nstart; iw<nend; ++iw)
      cout << "achtemp[" << iw << "] = " << std::setprecision(15) << achtemp[iw] << endl;

  cout << "********** Time Taken **********= " << elapsedTimer << " secs" << endl;

  //Free the allocated memory
  free(acht_n1_loc);
  free(wtilde_array);
  free(aqsmtemp);
  free(aqsntemp);
  free(I_eps_array);
  free(inv_igp_index);
  free(vcoul);

  return 0;
}

