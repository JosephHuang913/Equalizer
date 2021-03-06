/***********************************************************/
/* Author: Chao-wang Huang                                 */
/* Date: Saturday, April 02, 2005                          */
/* An adaptive MMSE Equalizer (Kalman Filter) is simulated */
/* 3-path Rayleigh fading channel with equal strength      */
/* Multipath weighting factor: {0.577, 0.577, 0.577}       */
/* Algorithm: Recursive Least Squares (RLS) Algorithm      */
/***********************************************************/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <conio.h>
#include <iostream.h>
#include <limits.h>
#include <float.h>

void MMSE_EQ_RLS(double **, int *, double **, int, int, int, double, double, int);
void JakesFading(double, float, double, int, double *);
void AWGN_noise(float, double, double *);
int Error_count(int, int);

#define m	2											// channel memory
const int N = 1000;
const int Training = 100;							// Length of training sequence
#define num_packet 1000								// number of packets simulated
const int Num_path = 3;
const int Num_state = pow(2,m);
const float lamda = 0.95;								// Weighting factor of Kalman Filter (0<lamda<1)
const int Num_tap = 11;	   						// Tap number of adaptive MMSE EQ
const int delay = 5;									// Delay of the adaptive MMSE EQ
float CIR[3] = {0.577, 0.577, 0.577};			// Channel Weighting Factor
//float CIR[3] = {0.407, 0.815, 0.407};
//float CIR[3] = {0.2096, 0.9551, 0.2096};
//float CIR[3] = {0.0, 1.0, 0.0};
//const  float Pi = 3.14159265358979;
//float CIR[3] = {(1+cos(-2*Pi/3.5))/2.0, (1+cos(0.0))/2.0, (1+cos(2*Pi/3.5))/2.0};
const float vc = 120.0; 								/* speed of vehicle in km/hr */
const double C = 3.0E8;  							/* light speed */
const double fc = 2.0e9; 							/* carrier frequency */
const double sym_rate = 1E6; 						/* Symbol rate in symbol/sec */
const double Doppler = (vc*1000.0/3600.0)/(C/fc);  // Maximum Doppler frequency (Hz)
double **W, **U, ***P, **Kalman;

int main(void)
{
	time_t  t, start, end;
	int i, j, p, *data_bit, *ch, from_state, tmp_state, err_count, *Hk, *train_sym;
   double out_sym_I, out_sym_Q, **Output, err_rate;
   double snr, Eb_No, noise_pwr, noise[2], **Yk;
   double t1, t2, t3, fade1[2], fade2[2], fade3[2];
   FILE *ber, *records;

   start = time(NULL);
   printf("MSE Performance of adaptive MMSE EQ (Kalman Filter) in Multipath Fading Channel\n");
	cout << "3-path Rayleigh fading channel with equal strength" << endl;
	cout << "Multipath weighting factor: {0.577, 0.577, 0.577}" << endl;
   printf("Speed of the vehicle = %f (km/h)\n", vc);
   printf("Carrier Frequency = %e (Hz)\n", fc);
   printf("Maximum Doppler Frequency = %f (Hz)\n", Doppler);
   printf("Transmission bit Rate = %e (bps)\n", sym_rate*1);
   printf("f_d * t = %f\n", Doppler / sym_rate);
   printf("number of bits of simulation = %d\n\n", N*num_packet);
   printf("\nThis program is running. Don't close, please!\n\n");

   records = fopen("Records_MMSE_design.log", "a");
   fprintf(records, "MSE Performance of adaptive MMSE EQ (Kalman Filter) in Multipath Fading Channel\n");
   fprintf(records, "3-path Rayleigh fading channel with equal strength\n");
   fprintf(records, "Multipath weighting factor: {0.577, 0.577, 0.577}\n");
   fprintf(records, "Speed of the vehicle = %f (km/h)\n", vc);
   fprintf(records, "Carrier Frequency = %e (Hz)\n", fc);
   fprintf(records, "Maximum Doppler Frequency = %f (Hz)\n", Doppler);
   fprintf(records, "Transmission bit Rate = %e (bps)\n", sym_rate*1);
   fprintf(records, "f_d * t = %f\n", Doppler / sym_rate);
   fprintf(records, "number of bits of simulation = %d\n\n", N*num_packet);
   fprintf(records, "Eb/No     BER\n");
   fflush(records);

   data_bit = new int[2*N];
   ch = new int[m+1];
   Hk = new int[N];			// Hard Decision
   train_sym = new int[Training];	// Training Symbols
   Yk = new double*[2];
   for(i=0; i<2; i++)
   	Yk[i] = new double[2*N];
   Output = new double*[2];			// Equalizer output
   for(i=0; i<2; i++)
   	Output[i] = new double[2*N];
   W = new double*[2];	// Coefficients of adaptive MMSE EQ
   for(i=0; i<2; i++)
   	W[i] = new double[Num_tap];
   U = new double*[2];	// Input signal of adaptive MMSE EQ
   for(i=0; i<2; i++)
   	U[i] = new double[Num_tap];
   P = new double**[2];
   for(i=0; i<2; i++)
   	P[i] = new double*[Num_tap];
   for(i=0; i<2; i++)
   	for(j=0; j<Num_tap; j++)
      	P[i][j] = new double[Num_tap];
   Kalman = new double*[2];
   for(i=0; i<2; i++)
   	Kalman[i] = new double[Num_tap];

	srand((unsigned) time(&t));

/************************/
/* main simulation loop */
/************************/
   ber = fopen("ber_MMSE.log", "w");
   for(snr=0; snr<=30; snr+=2)
   {
   	err_count = 0;
   	// noise power calculation
      Eb_No = (double)snr;
      noise_pwr = 1.0/pow(10.0, Eb_No/10.0);	// BPSK, Nyquist filter assumption
      printf("Eb_No = %f, ", Eb_No);
      t1 = 10.0;  	// Initial time of Rayleigh fading pattern
      t2 = 20.0;  	// Initial time of the 2nd path
      t3 = 30.0;		// Innitial time of the 3rd path
      from_state = random(4);						// Initial State of Channel
      //from_state = 0;

/*************************************/
/* Training Mode (Training Sequence) */
/*************************************/
		for(i=0; i<Training; i++)
      {
      	train_sym[i] = random(2);		// Generate random information bit stream
         tmp_state = from_state;
         tmp_state = (tmp_state << 1) ^ (train_sym[i] & 0x01);  // read input bit
         ch[0] = 2*(tmp_state & 0x01) - 1;					// input symbol
         ch[1] = 2*((tmp_state >> 1) & 0x01) - 1;        // channel memory
         ch[2] = 2*((tmp_state >> 2) & 0x01) - 1;        // channel memory

         JakesFading(fc, vc*1000/3600.0, t1, 2, &fade1[0]);
         t1 += 1.0/sym_rate;
         JakesFading(fc, vc*1000/3600.0, t2, 2, &fade2[0]);
         t2 += 1.0/sym_rate;
         JakesFading(fc, vc*1000/3600.0, t3, 2, &fade3[0]);
         t3 += 1.0/sym_rate;

         // Calculate channel output symbol
         out_sym_I = ch[0]*CIR[0]*fade1[0] + ch[1]*CIR[1]*fade2[0] + ch[2]*CIR[2]*fade3[0];
         out_sym_Q = ch[0]*CIR[0]*fade1[1] + ch[1]*CIR[1]*fade2[1] + ch[2]*CIR[2]*fade3[1];
         from_state = tmp_state & (Num_state-1); 			// to_state of trellis diagram

         /* AWGN channel */
         AWGN_noise(0, noise_pwr, &noise[0]);
         Yk[0][i] = out_sym_I + noise[0];
         Yk[1][i] = out_sym_Q + noise[1];
      }

/********************************************/
/* Decision Directed Mode (Data bit stream) */
/********************************************/
      for(i=Training; i<Training+N; i++)
      {
      	data_bit[i] = random(2);		// Generate random information bit stream
         tmp_state = from_state;
         tmp_state = (tmp_state << 1) ^ (data_bit[i] & 0x01);  // read input bit
         ch[0] = 2*(tmp_state & 0x01) - 1;					// input symbol
         ch[1] = 2*((tmp_state >> 1) & 0x01) - 1;        // channel memory
         ch[2] = 2*((tmp_state >> 2) & 0x01) - 1;        // channel memory

         JakesFading(fc, vc*1000/3600.0, t1, 2, &fade1[0]);
         t1 += 1.0/sym_rate;
         JakesFading(fc, vc*1000/3600.0, t2, 2, &fade2[0]);
         t2 += 1.0/sym_rate;
         JakesFading(fc, vc*1000/3600.0, t3, 2, &fade3[0]);
         t3 += 1.0/sym_rate;

         // Calculate channel output symbol
         out_sym_I = ch[0]*CIR[0]*fade1[0] + ch[1]*CIR[1]*fade2[0] + ch[2]*CIR[2]*fade3[0];
         out_sym_Q = ch[0]*CIR[0]*fade1[1] + ch[1]*CIR[1]*fade2[1] + ch[2]*CIR[2]*fade3[1];
         from_state = tmp_state & (Num_state-1); 			// to_state of trellis diagram

         /* AWGN channel */
         AWGN_noise(0, noise_pwr, &noise[0]);
         Yk[0][i] = out_sym_I + noise[0];
         Yk[1][i] = out_sym_Q + noise[1];
      }

      // Training Mode
      MMSE_EQ_RLS(Yk, train_sym, Output, Num_tap, Training, delay, lamda, noise_pwr, 1);

      // Copy the 2nd packet as the 1st packet
      for(i=0; i<N; i++)
      {
      	data_bit[i] = data_bit[i+Training];

         for(j=0; j<2; j++)
         	Yk[j][i] = Yk[j][i+Training];
      }

      for(p=0; p<num_packet; p++)
      {
   		for(i=N; i<2*N; i++)
         {
         	data_bit[i] = random(2);		// Generate random information bit stream
         	tmp_state = from_state;
            tmp_state = (tmp_state << 1) ^ (data_bit[i] & 0x01);  // read input bit
            ch[0] = 2*(tmp_state & 0x01) - 1;					// input symbol
            ch[1] = 2*((tmp_state >> 1) & 0x01) - 1;        // channel memory
            ch[2] = 2*((tmp_state >> 2) & 0x01) - 1;        // channel memory

            JakesFading(fc, vc*1000/3600.0, t1, 2, &fade1[0]);
      		t1 += 1.0/sym_rate;
            JakesFading(fc, vc*1000/3600.0, t2, 2, &fade2[0]);
      		t2 += 1.0/sym_rate;
            JakesFading(fc, vc*1000/3600.0, t3, 2, &fade3[0]);
      		t3 += 1.0/sym_rate;

            // Calculate channel output symbol
         	out_sym_I = ch[0]*CIR[0]*fade1[0] + ch[1]*CIR[1]*fade2[0] + ch[2]*CIR[2]*fade3[0];
            out_sym_Q = ch[0]*CIR[0]*fade1[1] + ch[1]*CIR[1]*fade2[1] + ch[2]*CIR[2]*fade3[1];
            from_state = tmp_state & (Num_state-1); 			// to_state of trellis diagram

            /* AWGN channel */
            AWGN_noise(0, noise_pwr, &noise[0]);
         	Yk[0][i] = out_sym_I + noise[0];
            Yk[1][i] = out_sym_Q + noise[1];
         }

         // Decision Directed Mode
	      MMSE_EQ_RLS(Yk, train_sym, Output, Num_tap, N, delay, lamda, noise_pwr, 0);

   	   // Copy the 2nd packet as the 1st packet
      	for(i=0; i<N; i++)
	      {
         	// Bit error count
         	if(Output[0][i]>=0)
            	Hk[i] = 1;
            else
            	Hk[i] = 0;

            err_count += Error_count(data_bit[i], Hk[i]);

   	   	data_bit[i] = data_bit[i+N];
      	   for(j=0; j<2; j++)
         		Yk[j][i] = Yk[j][i+N];
	      }
    	}

      // Statistics and records
      cout << "Error Rate = ";
      err_rate = err_count / (long double)(N*num_packet);
      printf("%e\n", err_rate);

      fprintf(ber, "%f ", Eb_No);
      fprintf(records, "%f ", Eb_No);
      fprintf(ber, "%e\n", err_rate);
      fprintf(records, "%e\n", err_rate);
      fflush(records);
      fflush(ber);
   }

   fclose(ber);
   delete data_bit;
   delete ch;
   delete Hk;
   delete train_sym;
   for(i=0; i<2; i++)
   	delete Yk[i];
   delete Yk;
   for(i=0; i<2; i++)
   	delete Output[i];
   delete Output;
   for(i=0; i<2; i++)
   	delete W[i];
   delete W;
   for(i=0; i<2; i++)
   	delete U[i];
   delete U;
   for(i=0; i<2; i++)
   	for(j=0; j<Num_tap; j++)
      	delete P[i][j];
   for(i=0; i<2; i++)
   	delete P[i];
   delete P;
   for(i=0; i<2; i++)
   	delete Kalman[i];
   delete Kalman;

   end = time(NULL);
   printf("Total elapsed time: %.0f(sec)\n", difftime(end,start));
   fprintf(records, "\nTotal elapsed time: %.0f(sec)\n\n", difftime(end,start));
   fclose(records);
   printf("This program is ended. Press any key to continue.\n");
   getch();

   return 0;
}

void MMSE_EQ_RLS(double **Received, int *train_sym, double **Output, int Num_tap,
					  int Packet_length, int delay, double lamda, double noise_pwr, int mode)
{
/*****************************************************************/
/* Adaptive MMSE EQ (Kalman Filter, RLS direct form)             */
/* mode = 1 for Training Mode, otherwise, Decision Directed Mode */
/*****************************************************************/
   int i, j, l, Hk;
   double e[2], Out[2], **pi, Denominator[2];

   pi = new double*[2];
   for(i=0; i<2; i++)
   	pi[i] = new double[Num_tap];

	if(mode==1)
   {
   	// Initialize the adaptive MMSE filter
      for(j=0; j<Num_tap; j++)
      {
      	U[0][j] = U[1][j] = 0.0;
         W[0][j] = W[1][j] = 0.0;

         for(l=0; l<Num_tap; l++)
         	P[0][j][l] = P[1][j][l] = 0.0;

         P[0][j][j] = 1.0/(pow(1-lamda,1)*(1.0+noise_pwr));
      }

      // Input signal of adaptive MMSE EQ
      for(j=0; j<delay; j++)
      {
      	for(l=Num_tap-1; l>0; l--)
      	{
      		U[0][l] = U[0][l-1];
	         U[1][l] = U[1][l-1];
   	   }
      	U[0][0] = Received[0][j];
	      U[1][0] = Received[1][j];
      }
   }

   for(i=0; i<Packet_length; i++)
   {
   	// Input signal of adaptive MMSE EQ
   	for(l=Num_tap-1; l>0; l--)
      {
      	U[0][l] = U[0][l-1];
         U[1][l] = U[1][l-1];
      }
      U[0][0] = Received[0][i+delay];
      U[1][0] = Received[1][i+delay];

      // Compute output
      Out[0] = Out[1] = 0.0;
      for(j=0; j<Num_tap; j++)
      {
      	Out[0] += (W[0][j]*U[0][j] + W[1][j]*U[1][j]);
         Out[1] += (W[0][j]*U[1][j] - W[1][j]*U[0][j]);
      }

      // Compute Kalman gain vector
      for(j=0; j<Num_tap; j++)
      	pi[0][j] = pi[1][j] = 0.0;

      for(j=0; j<Num_tap; j++) 		// row index
      	for(l=0; l<Num_tap; l++)   // column index
         {
         	pi[0][j] += (P[0][j][l]*U[0][l] - P[1][j][l]*U[1][l]);
            pi[1][j] += (P[1][j][l]*U[0][l] + P[0][j][l]*U[1][l]);
         }

      Denominator[0] = Denominator[1] = 0.0;
      for(j=0; j<Num_tap; j++)
      {
      	Denominator[0] += (U[0][j]*pi[0][j] + U[1][j]*pi[1][j]);
         Denominator[1] += (U[0][j]*pi[1][j] - U[1][j]*pi[0][j]);
      }

      for(j=0; j<Num_tap; j++)
      {
      	Kalman[0][j] = (pi[0][j]*(Denominator[0]+lamda)+pi[1][j]*Denominator[1])
         					/(pow(Denominator[0]+lamda,2)+pow(Denominator[1],2));
         Kalman[1][j] = (pi[1][j]*(Denominator[0]+lamda)-pi[0][j]*Denominator[1])
         					/(pow(Denominator[0]+lamda,2)+pow(Denominator[1],2));
      }

      // Compute error signal (BPSK)
      if(mode==1)
      {
      	e[0] = (2*train_sym[i]-1) - Out[0];
      	e[1] = 0.0 - Out[1];
      }
      else
      {
      	if(Out[0]>=0.0)
         	Hk = 1;
         else
         	Hk = -1;

         e[0] = Hk - Out[0];
         e[1] = 0.0 - Out[1];

         Output[0][i] = Out[0];
         Output[1][i] = Out[1];
      }

      // Update Coefficients
      for(j=0; j<Num_tap; j++)
      {
      	W[0][j] += (Kalman[0][j]*e[0] + Kalman[1][j]*e[1]);
         W[1][j] += (Kalman[1][j]*e[0] - Kalman[0][j]*e[1]);
      }

      // Update inverse of the correlation matrix
      for(j=0; j<Num_tap; j++)
      	pi[0][j] = pi[1][j] = 0.0;

      for(l=0; l<Num_tap; l++)      // column index
      	for(j=0; j<Num_tap; j++)	// row index
         {
         	pi[0][l] += (U[0][j]*P[0][j][l] + U[1][j]*P[1][j][l]);
            pi[1][l] += (U[0][j]*P[1][j][l] - U[1][j]*P[0][j][l]);
         }

      for(j=0; j<Num_tap; j++) 		// row index
      	for(l=0; l<Num_tap; l++)   // column index
         {
         	P[0][j][l] = (P[0][j][l] - (Kalman[0][j]*pi[0][l] - Kalman[1][j]*pi[1][l]))/lamda;
            P[1][j][l] = (P[1][j][l] - (Kalman[0][j]*pi[1][l] + Kalman[1][j]*pi[0][l]))/lamda;
         }
   }

   for(i=0; i<2; i++)
   	delete pi[i];
   delete pi;
}

void JakesFading(double f_c/*Hz*/, float v/*m/s*/, double t/*s*/, int type, double *fade)
{
	const double C = 3.0e8;     // (m/s)
   const float Pi = 3.14159265358979;
   int n, N, N_o = 32;
   double lamda, w_m, beta_n, w_n, alpha, T_c2, T_s2, theta_n;

   lamda = C/f_c;     // wave length (meter)
   w_m = 2.0*Pi*v/lamda;    // maximum Doppler frequency
   N = 2*(2*N_o+1);

   switch(type)
   {
   	case 1:
   		alpha = 0.0;
         T_c2 = (double)N_o;
         T_s2 = (double)N_o + 1.0;
         break;
      case 2:
      	alpha = 0.0;
         T_c2 = (double)N_o + 1.0;
         T_s2 = (double)N_o;
         break;
      case 3:
      	alpha = Pi/4.0;
         T_c2 = (double)N_o + 0.5;
         T_s2 = (double)N_o + 0.5;
         break;
      default:
      	printf("\nInvalid type selection for Jake's fading channel model.\n");
         break;
   }

   if(v == 0.0)
   {
   	*(fade+0) = 1.0;
      *(fade+1) = 0.0;
   }
   else
   {
   	*(fade+0) = sqrt(1.0/T_c2)*cos(alpha)*cos(w_m*t);
      *(fade+1) = sqrt(1.0/T_s2)*sin(alpha)*cos(w_m*t);

      for(n = 1; n <= N_o; n++)
      {
      	switch(type)
         {
         	case 1:
            	beta_n = (double)n*Pi/((double)N_o+1.0);
               break;
            case 2:
            	beta_n = (double)n*Pi/(double)N_o;
               break;
            case 3:
            	beta_n = (double)n*Pi/(double)N_o;
               break;
         	default:
            	break;
         }
         w_n = w_m*cos(2.0*Pi*(double)n/(double)N);
//            theta_n = 2.0*Pi*((double)rand()/(double)RAND_MAX);  // random phase
			theta_n = 0.0;
         *(fade+0) += sqrt(2.0/T_c2)*cos(beta_n)*cos(w_n*t+theta_n);
         *(fade+1) += sqrt(2.0/T_s2)*sin(beta_n)*cos(w_n*t+theta_n);
		}
	}
}

void AWGN_noise(float mu, double variance, double *noise)
{
	const  float Pi = 3.14159265358979;
   double u1, u2;
   do
   {
   	u1 = (double)rand()/(double)RAND_MAX;
      u2 = (double)rand()/(double)RAND_MAX;
   }
   while(u1 == 0.0 || u2 == 0.0);

   *(noise+0) = (sqrt(-2.0*log(u1))*cos(2*Pi*u2))*sqrt(variance/2.0)+mu/sqrt(2.0);
   *(noise+1) = (sqrt(-2.0*log(u1))*sin(2*Pi*u2))*sqrt(variance/2.0)+mu/sqrt(2.0);
}

int Error_count(int x, int y)
{
	if(x == y)
   	return 0;
   else
   	return 1;
}

