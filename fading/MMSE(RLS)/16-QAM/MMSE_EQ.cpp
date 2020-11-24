/********************************************************************************/
/* Author: Chao-wang Huang                                                      */
/* Date: Wednesday, November 16, 2005                                           */
/* An adaptive MMSE Equalizer (Kalman Filter) is simulated                      */
/* Proakis A channel: {0.04 -0.05 0.07 -0.21 -0.5 0.72 0.36 0.0 0.21 0.03 0.07} */
/* Algorithm: Recursive Least Squares (RLS) Algorithm                           */
/* Modulation Scheme: 16-QAM                                                    */
/********************************************************************************/

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

const int N = 1000;
const int Training = 100;							// Length of training sequence
#define num_packet 1000								// number of packets simulated
const int Num_path = 11;
const float lamda = 0.999;								// Weighting factor of Kalman Filter (0<lamda<1)
const int Num_tap = 11;	   						// Tap number of adaptive MMSE EQ
const int delay = 9;									// Delay of the adaptive MMSE EQ
float CIR[11] = {0.04,-0.05,0.07,-0.21,-0.5,0.72,0.36,0.0,0.21,0.03,0.07};
//float CIR[11] = {0.0,  0.0,  0.0,  0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0};
const float vc = 0.0; 								/* speed of vehicle in km/hr */
const double C = 3.0E8;  							/* light speed */
const double fc = 2.0e9; 							/* carrier frequency */
const double sym_rate = 1E6; 						/* Symbol rate in symbol/sec */
const double Doppler = (vc*1000.0/3600.0)/(C/fc);  // Maximum Doppler frequency (Hz)
double **W, **U, ***P, **Kalman;

int main(void)
{
	time_t  t, start, end;
	int i, j, p, *data_bit, err_count, *Hk, *train_sym;
   double **Output, err_rate, *sym_I, *sym_Q;
   double snr, Eb_No, noise_pwr, noise[2], **Yk, **ch_matrix_I, **ch_matrix_Q;
   FILE *ber, *records;

   start = time(NULL);
   printf("BER Performance of adaptive MMSE EQ (Kalman Filter) in Multipath Static Channel\n");
	cout << "Proakis A channel is used" << endl;
	cout << "Multipath weighting factor: {0.04 -0.05 0.07 -0.21 -0.5 0.72 0.36 0.0 0.21 0.03 0.07}" << endl;
   printf("Speed of the vehicle = %f (km/h)\n", vc);
   printf("Carrier Frequency = %e (Hz)\n", fc);
   printf("Maximum Doppler Frequency = %f (Hz)\n", Doppler);
   printf("Transmission bit Rate = %e (bps)\n", sym_rate*1);
   printf("f_d * t = %f\n", Doppler / sym_rate);
   printf("number of bits of simulation = %d\n\n", N*num_packet);
   printf("\nThis program is running. Don't close, please!\n\n");

   records = fopen("Records_MMSE_ber.log", "a");
   fprintf(records, "BER Performance of adaptive MMSE EQ (Kalman Filter) in Multipath Fading Channel\n");
   fprintf(records, "Proakis A channel is used\n");
   fprintf(records, "Multipath weighting factor: {0.04 -0.05 0.07 -0.21 -0.5 0.72 0.36 0.0 0.21 0.03 0.07}\n");
   fprintf(records, "Speed of the vehicle = %f (km/h)\n", vc);
   fprintf(records, "Carrier Frequency = %e (Hz)\n", fc);
   fprintf(records, "Maximum Doppler Frequency = %f (Hz)\n", Doppler);
   fprintf(records, "Transmission bit Rate = %e (bps)\n", sym_rate*1);
   fprintf(records, "f_d * t = %f\n", Doppler / sym_rate);
   fprintf(records, "number of bits of simulation = %d\n\n", N*num_packet);
   fprintf(records, "Eb/No     BER\n");
   fflush(records);

   data_bit = new int[2*N];
   Hk = new int[4];			// Hard Decision
   train_sym = new int[Training];	// Training Symbols
   sym_I = new double[N/2];
   sym_Q = new double[N/2];
   ch_matrix_I = new double*[Num_path];
   for(i=0; i<Num_path; i++)
   	ch_matrix_I[i] = new double[N/2+Num_path-1];
   ch_matrix_Q = new double*[Num_path];
   for(i=0; i<Num_path; i++)
   	ch_matrix_Q[i] = new double[N/2+Num_path-1];
   Yk = new double*[2];
   for(i=0; i<2; i++)
   	Yk[i] = new double[N/2];
   Output = new double*[2];			// Equalizer output
   for(i=0; i<2; i++)
   	Output[i] = new double[N/4];
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
   for(j=0; j<Num_path; j++)
   	for(i=0; i<N/2+Num_path-1; i++)
      {
        	ch_matrix_I[j][i] = 0.0;
         ch_matrix_Q[j][i] = 0.0;
      }

/************************/
/* main simulation loop */
/************************/
   ber = fopen("ber_MMSE_16QAM.log", "w");
   for(snr=0; snr<=30; snr+=5)
   {
   	err_count = 0;
   	// noise power calculation
      Eb_No = (double)snr;
      noise_pwr = 0.25/pow(10.0, Eb_No/10.0);	// 16-QAM, Nyquist filter assumption
      printf("Eb_No = %f, ", Eb_No);

/*************************************/
/* Training Mode (Training Sequence) */
/*************************************/
		for(i=0; i<Training; i++)
      	train_sym[i] = random(2);		// Generate random information bit stream

      // 16-QAM Mapping
      for(i=0; i<Training/4; i++)
      {
      	if(train_sym[4*i] == 0 && train_sym[4*i+2] == 0)
   	  		sym_I[i] = 1.0/sqrt(10.0);
	       else if(train_sym[4*i] == 0 && train_sym[4*i+2] == 1)
   	  		sym_I[i] = 3.0/sqrt(10.0);
	     	else if(train_sym[4*i] == 1 && train_sym[4*i+2] == 1)
       		sym_I[i] = -3.0/sqrt(10.0);
	  	   else
  	   		sym_I[i] = -1.0/sqrt(10.0);

       	if(train_sym[4*i+1] == 0 && train_sym[4*i+3] == 0)
  	   		sym_Q[i] = 1.0/sqrt(10.0);
         else if(train_sym[4*i+1] == 0 && train_sym[4*i+3] == 1)
   	  		sym_Q[i] = 3.0/sqrt(10.0);
	     	else if(train_sym[4*i+1] == 1 && train_sym[4*i+3] == 1)
   	  		sym_Q[i] = -3.0/sqrt(10.0);
	  	   else
       		sym_Q[i] = -1.0/sqrt(10.0);
      }

      // Multipath Channel
      for(j=0; j<Num_path; j++)
	      for(i=0; i<Training/4; i++)
         {
   	   	ch_matrix_I[j][i+j] = sym_I[i] * CIR[j];
            ch_matrix_Q[j][i+j] = sym_Q[i] * CIR[j];
         }

/********************************************/
/* Decision Directed Mode (Data bit stream) */
/********************************************/
		for(i=Training; i<Training+N; i++)
         data_bit[i] = random(2);		// Generate random information bit stream

      for(i=Training/4; i<(Training+N)/4; i++)
      {
      	if(data_bit[4*i] == 0 && data_bit[4*i+2] == 0)
   	  		sym_I[i] = 1.0/sqrt(10.0);
	       else if(data_bit[4*i] == 0 && data_bit[4*i+2] == 1)
   	  		sym_I[i] = 3.0/sqrt(10.0);
	     	else if(data_bit[4*i] == 1 && data_bit[4*i+2] == 1)
       		sym_I[i] = -3.0/sqrt(10.0);
	  	   else
  	   		sym_I[i] = -1.0/sqrt(10.0);

       	if(data_bit[4*i+1] == 0 && data_bit[4*i+3] == 0)
  	   		sym_Q[i] = 1.0/sqrt(10.0);
         else if(data_bit[4*i+1] == 0 && data_bit[4*i+3] == 1)
   	  		sym_Q[i] = 3.0/sqrt(10.0);
	     	else if(data_bit[4*i+1] == 1 && data_bit[4*i+3] == 1)
   	  		sym_Q[i] = -3.0/sqrt(10.0);
	  	   else
       		sym_Q[i] = -1.0/sqrt(10.0);
      }

      // Multipath Channel
      for(j=0; j<Num_path; j++)
	      for(i=Training/4; i<(Training+N)/4; i++)
         {
   	   	ch_matrix_I[j][i+j] = sym_I[i] * CIR[j];
            ch_matrix_Q[j][i+j] = sym_Q[i] * CIR[j];
         }

      for(i=0; i<N/2; i++)
      	for(j=0; j<2; j++)
         	Yk[j][i] = 0.0;

      for(i=0; i<(Training+N)/4; i++)
      	for(j=0; j<Num_path; j++)
         {
         	Yk[0][i] += ch_matrix_I[j][i];
            Yk[1][i] += ch_matrix_Q[j][i];
         }

      // AWGN Channel
      for(i=0; i<(Training+N)/4; i++)
      {
      	AWGN_noise(0, noise_pwr, &noise[0]);
         Yk[0][i] += noise[0];
         Yk[1][i] += noise[1];
      }

      // Training Mode
      MMSE_EQ_RLS(Yk, train_sym, Output, Num_tap, Training/4, delay, lamda, 0.0, 1);

      // Copy the 2nd packet as the 1st packet
      for(i=0; i<N; i++)
      	data_bit[i] = data_bit[i+Training];

      for(i=0; i<N/4; i++)
         for(j=0; j<2; j++)
         	Yk[j][i] = Yk[j][i+Training/4];

      for(j=0; j<Num_path; j++)
      	for(i=0; i<N/4+j; i++)
         {
            ch_matrix_I[j][i] = ch_matrix_I[j][i+Training/4];
            ch_matrix_Q[j][i] = ch_matrix_Q[j][i+Training/4];
         }

      // Main Simulation Loop
      for(p=0; p<num_packet; p++)
      {
   		for(i=N; i<2*N; i++)
         	data_bit[i] = random(2);		// Generate random information bit stream

         for(i=N/4; i<N/2; i++)
      	{
      		if(data_bit[4*i] == 0 && data_bit[4*i+2] == 0)
   	  			sym_I[i] = 1.0/sqrt(10.0);
		       else if(data_bit[4*i] == 0 && data_bit[4*i+2] == 1)
   		  		sym_I[i] = 3.0/sqrt(10.0);
	   	  	else if(data_bit[4*i] == 1 && data_bit[4*i+2] == 1)
       			sym_I[i] = -3.0/sqrt(10.0);
		  	   else
  		   		sym_I[i] = -1.0/sqrt(10.0);

      	 	if(data_bit[4*i+1] == 0 && data_bit[4*i+3] == 0)
  	   			sym_Q[i] = 1.0/sqrt(10.0);
	         else if(data_bit[4*i+1] == 0 && data_bit[4*i+3] == 1)
   		  		sym_Q[i] = 3.0/sqrt(10.0);
	   	  	else if(data_bit[4*i+1] == 1 && data_bit[4*i+3] == 1)
   	  			sym_Q[i] = -3.0/sqrt(10.0);
		  	   else
   	    		sym_Q[i] = -1.0/sqrt(10.0);
      	}

      	// Multipath Channel
	      for(j=0; j<Num_path; j++)
		      for(i=N/4; i<N/2; i++)
      	   {
   	   		ch_matrix_I[j][i+j] = sym_I[i] * CIR[j];
            	ch_matrix_Q[j][i+j] = sym_Q[i] * CIR[j];
	         }

      	for(i=N/4; i<N/2; i++)
      		for(j=0; j<2; j++)
         		Yk[j][i] = 0.0;

   	   for(i=N/4; i<N/2; i++)
      		for(j=0; j<Num_path; j++)
         	{
         		Yk[0][i] += ch_matrix_I[j][i];
	            Yk[1][i] += ch_matrix_Q[j][i];
   	      }

      	// AWGN Channel
	      for(i=N/4; i<N/2; i++)
   	   {
      		AWGN_noise(0, noise_pwr, &noise[0]);
         	Yk[0][i] += noise[0];
	         Yk[1][i] += noise[1];
   	   }

         // Decision Directed Mode
	      MMSE_EQ_RLS(Yk, train_sym, Output, Num_tap, N/4, delay, lamda, 0.0, 0);

         // 16-QAM De-mapping (Hard Decision)
      	for(i=0; i<N/4; i++)
	      {
            if(Output[0][i] < -2.0/sqrt(10.0))
      		{
      			Hk[0] = 1;
		         Hk[2] = 1;
   		   }
      		else if(Output[0][i] >= -2.0/sqrt(10.0) && Output[0][i] < 0.0)
	      	{
   	   		Hk[0] = 1;
	      	   Hk[2] = 0;
		      }
   		   else if(Output[0][i] >= 0.0 && Output[0][i] < 2.0/sqrt(10.0))
      		{
      			Hk[0] = 0;
		         Hk[2] = 0;
   		   }
      		else
	      	{
   	   		Hk[0] = 0;
	      	   Hk[2] = 1;
		      }

   	   	if(Output[1][i] < -2.0/sqrt(10.0))
	      	{
   	   		Hk[1] = 1;
	   	      Hk[3] = 1;
   	   	}
	      	else if(Output[1][i] >= -2.0/sqrt(10.0) && Output[1][i] < 0.0)
		      {
   		   	Hk[1] = 1;
      		   Hk[3] = 0;
		      }
   		   else if(Output[1][i] >= 0.0 && Output[1][i] < 2.0/sqrt(10.0))
      		{
      			Hk[1] = 0;
	         	Hk[3] = 0;
	   	   }
   	   	else
	   	   {
   	   		Hk[1] = 0;
      	   	Hk[3] = 1;
		      }

            // Bit error count
            for(j=0; j<4; j++)
            	err_count += Error_count(data_bit[4*i+j], Hk[j]);
         }

         // Copy the 2nd packet as the 1st packet
         for(i=0; i<N; i++)
   	   	data_bit[i] = data_bit[i+N];

         for(i=0; i<N/4; i++)
      	   for(j=0; j<2; j++)
         		Yk[j][i] = Yk[j][i+N/4];

         for(j=0; j<Num_path; j++)
      		for(i=0; i<N/4+j; i++)
	         {
   	         ch_matrix_I[j][i] = ch_matrix_I[j][i+N/4];
      	      ch_matrix_Q[j][i] = ch_matrix_Q[j][i+N/4];
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
   delete Hk;
   delete train_sym;
   delete sym_I;
   delete sym_Q;
   for(i=0; i<Num_path; i++)
   	delete ch_matrix_I[i];
   delete ch_matrix_I;
   for(i=0; i<Num_path; i++)
   	delete ch_matrix_Q[i];
   delete ch_matrix_Q;
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
/* Adaptive MMSE EQ (Kalman Filter, RLS direct form) for 16-QAM  */
/* mode = 1 for Training Mode, otherwise, Decision Directed Mode */
/*****************************************************************/
   int i, j, l, detected[4];
   double sym_I, sym_Q, e[2], Out[2], **pi, Denominator[2];

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

      // Compute error signal (16-QAM)
      if(mode==1)
      {
      	// 16-QAM Mapping
         if(train_sym[4*i] == 0 && train_sym[4*i+2] == 0)
         	sym_I = 1.0/sqrt(10.0);
         else if(train_sym[4*i] == 0 && train_sym[4*i+2] == 1)
         	sym_I = 3.0/sqrt(10.0);
         else if(train_sym[4*i] == 1 && train_sym[4*i+2] == 1)
         	sym_I = -3.0/sqrt(10.0);
         else
         	sym_I = -1.0/sqrt(10.0);

         if(train_sym[4*i+1] == 0 && train_sym[4*i+3] == 0)
         	sym_Q = 1.0/sqrt(10.0);
         else if(train_sym[4*i+1] == 0 && train_sym[4*i+3] == 1)
         	sym_Q = 3.0/sqrt(10.0);
         else if(train_sym[4*i+1] == 1 && train_sym[4*i+3] == 1)
         	sym_Q = -3.0/sqrt(10.0);
         else
         	sym_Q = -1.0/sqrt(10.0);

      	e[0] = sym_I - Out[0];
      	e[1] = sym_Q - Out[1];
      }
      else
      {
         // 16-QAM De-mapping (Hard Decision)
         if(Out[0] < -2.0/sqrt(10.0))
      	{
      		detected[0] = 1;
	         detected[2] = 1;
   	   }
      	else if(Out[0] >= -2.0/sqrt(10.0) && Out[0] < 0.0)
	      {
   	   	detected[0] = 1;
      	   detected[2] = 0;
	      }
   	   else if(Out[0] >= 0.0 && Out[0] < 2.0/sqrt(10.0))
      	{
      		detected[0] = 0;
	         detected[2] = 0;
   	   }
      	else
	      {
   	   	detected[0] = 0;
      	   detected[2] = 1;
	      }

   	   if(Out[1] < -2.0/sqrt(10.0))
      	{
      		detected[1] = 1;
	         detected[3] = 1;
   	   }
      	else if(Out[1] >= -2.0/sqrt(10.0) && Out[1] < 0.0)
	      {
   	   	detected[1] = 1;
      	   detected[3] = 0;
	      }
   	   else if(Out[1] >= 0.0 && Out[1] < 2.0/sqrt(10.0))
      	{
      		detected[1] = 0;
	         detected[3] = 0;
   	   }
      	else
	      {
   	   	detected[1] = 0;
      	   detected[3] = 1;
	      }

         // 16-QAM Mapping
         if(detected[0] == 0 && detected[2] == 0)
         	sym_I = 1.0/sqrt(10.0);
         else if(detected[0] == 0 && detected[2] == 1)
         	sym_I = 3.0/sqrt(10.0);
         else if(detected[0] == 1 && detected[2] == 1)
         	sym_I = -3.0/sqrt(10.0);
         else
         	sym_I = -1.0/sqrt(10.0);

         if(detected[1] == 0 && detected[3] == 0)
         	sym_Q = 1.0/sqrt(10.0);
         else if(detected[1] == 0 && detected[3] == 1)
         	sym_Q = 3.0/sqrt(10.0);
         else if(detected[1] == 1 && detected[3] == 1)
         	sym_Q = -3.0/sqrt(10.0);
         else
         	sym_Q = -1.0/sqrt(10.0);

      	e[0] = sym_I - Out[0];
      	e[1] = sym_Q - Out[1];

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

