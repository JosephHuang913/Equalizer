/*================================================================*/
/* Author: Chao-wang Huang                                        */
/* Date: Thursday, July 20, 2006                                  */
/* An adaptive RLS Equalizer (Kalman Filter) is simulated         */
/* 3-path Rayleigh fading channel with equal strength             */
/* Static ISI Channel: Proakis B Channel {0.407,0.815,0.407}      */
/* Multipath fading channel: Jakes model {0.577, 0.577, 0.577}    */
/*================================================================*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <conio.h>
#include <iostream.h>
#include <limits.h>
#include <float.h>

void MMSE_EQ_RLS(double **, int *, double **, int, int, double, int);
void JakesFading(double, float, double, int, double *);
void AWGN_noise(float, double, double *);
int Error_count(int, int);

#define num_packet 1000								// number of packets simulated
const int Training = 25;							// Length of training sequence 25 symbols (4-QAM)
const int Slot_length = 100;						// Slot length: 100 symbols (4-QAM)
const int Num_slot = 164;							// Number of slots of one packet
const float lamda = 0.965;							// Weighting factor of Kalman Filter (0<lamda<1)
const int Num_tap = 9;	   						// Tap number of adaptive MMSE EQ
const int delay = 5;									// Delay of the adaptive MMSE EQ
const int Num_path = 3;
double CIR[3] = {0.577, 0.577, 0.577};			// Channel Weighting Factor
//double CIR[3] = {0.407, 0.815, 0.407};
const float vc = 135.0; 								/* speed of vehicle in km/hr */
const double C = 3.0E8;  							/* light speed */
const double fc = 2.0e9; 							/* carrier frequency */
const double sym_rate = 5E5; 						/* Symbol rate in symbol/sec */
const double Doppler = (vc*1000.0/3600.0)/(C/fc);  // Maximum Doppler frequency (Hz)
double **W, **U, ***P, **Kalman;					// For adaptive RLS algorithm (Kalman Filter)
//double **EQ_out;
double noise_pwr;

int main(void)
{
	time_t  t, start, end;
	int i, p, l, j, x, *data_bit, err_count, *train_sym, *Hk;
   double ***fade_pattern, **ch_matrix_Q;
   double snr, Eb_No, /*noise_pwr,*/ noise[2], **Yk, err_rate, **ch_matrix_I;
   double *sym_I, *sym_Q, *Ak, t1, t2, t3, fade1[2], fade2[2], fade3[2], **EQ_out;
   FILE *ber, *records;

   start = time(NULL);
   printf("BER Performance of Adaptive RLS EQ (4-QAM) in multipath ISI channel\n");
   printf("Modulation Scheme: 4-QAM\n");
	cout << "3-path Rayleigh fading channel with equal strength" << endl;
	cout << "Multipath weighting factor: {" << CIR[0] << "," << CIR[1] << "," << CIR[2] << "}" << endl;
	printf("Speed of the vehicle = %f (km/h)\n", vc);
   printf("Carrier Frequency = %e (Hz)\n", fc);
   printf("Maximum Doppler Frequency = %f (Hz)\n", Doppler);
   printf("Transmission bit Rate = %e (bps)\n", sym_rate*1);
   printf("f_d * t = %f\n", Doppler / sym_rate);
	//printf("Maximum Number of bits of simulation = %d\n", (K-m)*num_packet);
   printf("This program is running. Don't close, please!\n\n");

   records = fopen("Records_RLS_4QAM.log", "a");
   fprintf(records, "BER Performance of Adaptive RLS EQ (4-QAM) in multipath ISI channel\n");
   fprintf(records, "Modulation Scheme: 4-QAM\n");
   fprintf(records, "3-path Rayleigh fading channel with equal strength\n");
   fprintf(records, "Multipath weighting factor: {%.3f, %.3f, %.3f}\n", CIR[0], CIR[1], CIR[2]);
   fprintf(records, "Speed of the vehicle = %f (km/h)\n", vc);
   fprintf(records, "Carrier Frequency = %e (Hz)\n", fc);
   fprintf(records, "Maximum Doppler Frequency = %f (Hz)\n", Doppler);
   fprintf(records, "Transmission bit Rate = %e (bps)\n", sym_rate*1);
   fprintf(records, "f_d * t = %f\n", Doppler / sym_rate);
   //fprintf(records, "Maximum Number of bits of simulation = %d\n", (K-m)*num_packet);
   fprintf(records, "Eb/No     BER\n");
   fflush(records);

   data_bit = new int[2*Num_slot*Slot_length*2];
   Ak = new double[Num_slot*Slot_length*2];                    // Soft Decision data bit
   Hk = new int[Num_slot*Slot_length*2];								// Hard Decision Output
   EQ_out = new double*[Num_slot*(Slot_length+Training)];		// Output of adaptive RLS EQ
   for(i=0; i<Num_slot*(Slot_length+Training); i++)
   	EQ_out[i] = new double[2];
   train_sym = new int[2*(Num_slot*Training*2)];		 			// Training Symbols
   Yk = new double*[2*Num_slot*(Slot_length+Training)];  		// Received signal
   for(i=0; i<2*Num_slot*(Slot_length+Training); i++)
   	Yk[i] = new double[2];
   fade_pattern = new double **[Num_path];
   for(i=0; i<Num_path; i++)
   	fade_pattern[i] = new double*[2*Num_slot*(Slot_length+Training)];
   for(i=0; i<Num_path; i++)
   	for(l=0; l<2*Num_slot*(Slot_length+Training); l++)
   		fade_pattern[i][l] = new double[2];
   sym_I = new double[2*Num_slot*(Slot_length+Training)];
   sym_Q = new double[2*Num_slot*(Slot_length+Training)];
   ch_matrix_I = new double*[Num_path];
   for(i=0; i<Num_path; i++)
   	ch_matrix_I[i] = new double[2*Num_slot*(Slot_length+Training)+Num_path-1];
   ch_matrix_Q = new double*[Num_path];
   for(i=0; i<Num_path; i++)
   	ch_matrix_Q[i] = new double[2*Num_slot*(Slot_length+Training)+Num_path-1];
   W = new double*[Num_tap];							// Coefficients of adaptive RLS EQ
   for(i=0; i<Num_tap; i++)
   	W[i] = new double[2];
   U = new double*[Num_tap];							// Input signal of adaptive RLS EQ
   for(i=0; i<Num_tap; i++)
   	U[i] = new double[2];
   P = new double**[Num_tap];
   for(i=0; i<Num_tap; i++)
   	P[i] = new double*[Num_tap];
   for(i=0; i<Num_tap; i++)
   	for(j=0; j<Num_tap; j++)
      	P[i][j] = new double[2];
   Kalman = new double*[Num_tap];
   for(i=0; i<Num_tap; i++)
   	Kalman[i] = new double[2];

	srand((unsigned) time(&t));

/************************/
/* main simulation loop */
/************************/
   ber=fopen("ber_RLS_EQ.log", "w");
   for(snr=0; snr<=30; snr+=2)
   {
   	err_count = 0;
   	// noise power calculation
      Eb_No = (double)snr;
      noise_pwr = 0.5/(pow(10.0, Eb_No/10.0));	// 4-QAM, Nyquist filter assumption, Non-coding
      printf("Eb_No = %f,\n", Eb_No);
      t1 = 100.0;  	// Initial time of Rayleigh fading pattern
      t2 = 200.0;  	// Initial time of the 2nd path
      t3 = 300.0;		// Innitial time of the 3rd path

      for(j=0; j<Num_path; j++)
	   	for(i=0; i<2*Num_slot*(Slot_length+Training)+Num_path-1; i++)
   	   {
      		ch_matrix_I[j][i] = 0.0;
         	ch_matrix_Q[j][i] = 0.0;
	      }

      // Initialize the adaptive RLS filter
      for(j=0; j<Num_tap; j++)
      {
      	U[j][0] = U[j][1] = 0.0;
         W[j][0] = W[j][1] = 0.0;

         for(l=0; l<Num_tap; l++)
         	P[j][l][0] = P[j][l][1] = 0.0;

         P[j][j][0] = 1.0/(pow(1-lamda,1)*(1.0+noise_pwr));
      }
      W[delay-1][0] = 1.0;		// Central coefficients is initialized to 1.0

/*=====================================================================*/
/*           F i r s t     p a c k e t											  */
/*=====================================================================*/
      // Generate random information bit stream for the 1st packet
      for(i=0; i<Num_slot*Slot_length*2; i++)
      	data_bit[i] = random(2);		// Generate random information bit stream

      // Generate Training symbols: 25 symbols per-slot
      for(i=0; i<Num_slot*Training*2; i++)
      	train_sym[i] = random(2);

/********************************************************/
/* 4-QAM mapping and ISI Channel: {0.577, 0.577, 0.577} */
/********************************************************/
		// 4-QAM Mapping
      for(i=0; i<Num_slot; i++)
      {
      	for(j=0; j<Training; j++)
         {
      		sym_I[j+i*(Training+Slot_length)] = (2*train_sym[2*(j+i*Training)]-1)/sqrt(2.0);
         	sym_Q[j+i*(Training+Slot_length)] = (2*train_sym[2*(j+i*Training)+1]-1)/sqrt(2.0);
         }

         for(j=0; j<Slot_length; j++)
         {
         	sym_I[j+i*(Training+Slot_length)+Training] = (2*data_bit[2*(j+i*Slot_length)]-1)/sqrt(2.0);
         	sym_Q[j+i*(Training+Slot_length)+Training] = (2*data_bit[2*(j+i*Slot_length)+1]-1)/sqrt(2.0);
         }
      }

      for(i=0; i<Num_slot*(Slot_length+Training); i++)
      {
			JakesFading(fc, vc*1000/3600.0, t1, 2, &fade1[0]);
      	t1 += 1.0/sym_rate;
      	JakesFading(fc, vc*1000/3600.0, t2, 2, &fade2[0]);
      	t2 += 1.0/sym_rate;
      	JakesFading(fc, vc*1000/3600.0, t3, 2, &fade3[0]);
      	t3 += 1.0/sym_rate;

         // Record the channel fading pattern for MPIC
         for(j=0; j<2; j++)
         {
         	fade_pattern[0][i][j] = fade1[j];
         	fade_pattern[1][i][j] = fade2[j];
         	fade_pattern[2][i][j] = fade3[j];
         }
      }

      // Multipath Channel
      for(j=0; j<Num_path; j++)
	      for(i=0; i<Num_slot*(Slot_length+Training); i++)
         {
   	   	ch_matrix_I[j][i+j] = CIR[j]*(sym_I[i]*fade_pattern[j][i][0] - sym_Q[i]*fade_pattern[j][i][1]);
            ch_matrix_Q[j][i+j] = CIR[j]*(sym_I[i]*fade_pattern[j][i][1] + sym_Q[i]*fade_pattern[j][i][0]);
         }

      for(i=0; i<Num_slot*(Slot_length+Training); i++)
      	for(j=0; j<2; j++)
         	Yk[i][j] = 0.0;

      for(i=0; i<Num_slot*(Slot_length+Training); i++)
      	for(j=0; j<Num_path; j++)
         {
         	Yk[i][0] += ch_matrix_I[j][i];
            Yk[i][1] += ch_matrix_Q[j][i];
         }

      // AWGN Channel
      for(i=0; i<Num_slot*(Slot_length+Training); i++)
      {
      	AWGN_noise(0, noise_pwr, &noise[0]);
         Yk[i][0] += noise[0];
         Yk[i][1] += noise[1];
      }

/*=====================================================================*/
/*      S e c o n d     p a c k e t												  */
/*=====================================================================*/
      p = 0;
      do
      //for(p=0; p<num_packet; p++)
      {
      	// Generate random information bit stream for the 2nd packet
      	for(i=Num_slot*Slot_length*2; i<2*Num_slot*Slot_length*2; i++)
   	   	data_bit[i] = random(2);		// Generate random information bit stream

      	// Generate Training symbols: 25 symbols per-slot
	      for(i=Num_slot*Training*2; i<2*Num_slot*Training*2; i++)
   	   	train_sym[i] = random(2);

/********************************************************/
/* 4-QAM mapping and ISI Channel: {0.577, 0.577, 0.577} */
/********************************************************/
         // 4-QAM Mapping
      	for(i=Num_slot; i<2*Num_slot; i++)
	      {
   	   	for(j=0; j<Training; j++)
      	   {
      			sym_I[j+i*(Training+Slot_length)] = (2*train_sym[2*(j+i*Training)]-1)/sqrt(2.0);
         		sym_Q[j+i*(Training+Slot_length)] = (2*train_sym[2*(j+i*Training)+1]-1)/sqrt(2.0);
	         }

   	      for(j=0; j<Slot_length; j++)
      	   {
         		sym_I[j+i*(Training+Slot_length)+Training] = (2*data_bit[2*(j+i*Slot_length)]-1)/sqrt(2.0);
         		sym_Q[j+i*(Training+Slot_length)+Training] = (2*data_bit[2*(j+i*Slot_length)+1]-1)/sqrt(2.0);
	         }
   	   }

      	for(i=Num_slot*(Slot_length+Training); i<2*Num_slot*(Slot_length+Training); i++)
      	{
				JakesFading(fc, vc*1000/3600.0, t1, 2, &fade1[0]);
   	   	t1 += 1.0/sym_rate;
      		JakesFading(fc, vc*1000/3600.0, t2, 2, &fade2[0]);
      		t2 += 1.0/sym_rate;
	      	JakesFading(fc, vc*1000/3600.0, t3, 2, &fade3[0]);
   	   	t3 += 1.0/sym_rate;

         	// Record the channel fading pattern for MPIC
	         for(j=0; j<2; j++)
   	      {
      	   	fade_pattern[0][i][j] = fade1[j];
         		fade_pattern[1][i][j] = fade2[j];
         		fade_pattern[2][i][j] = fade3[j];
	         }
   	   }

      	// Multipath Channel
	      for(j=0; j<Num_path; j++)
		      for(i=Num_slot*(Slot_length+Training); i<2*Num_slot*(Slot_length+Training); i++)
      	   {
   	   		ch_matrix_I[j][i+j] = CIR[j]*(sym_I[i]*fade_pattern[j][i][0] - sym_Q[i]*fade_pattern[j][i][1]);
            	ch_matrix_Q[j][i+j] = CIR[j]*(sym_I[i]*fade_pattern[j][i][1] + sym_Q[i]*fade_pattern[j][i][0]);
	         }

	      for(i=Num_slot*(Slot_length+Training); i<2*Num_slot*(Slot_length+Training); i++)
   	   	for(j=0; j<2; j++)
      	   	Yk[i][j] = 0.0;

	      for(i=Num_slot*(Slot_length+Training); i<2*Num_slot*(Slot_length+Training); i++)
   	   	for(j=0; j<Num_path; j++)
      	   {
         		Yk[i][0] += ch_matrix_I[j][i];
            	Yk[i][1] += ch_matrix_Q[j][i];
	         }

      	// AWGN Channel
	      for(i=Num_slot*(Slot_length+Training); i<2*Num_slot*(Slot_length+Training); i++)
   	   {
      		AWGN_noise(0, noise_pwr, &noise[0]);
         	Yk[i][0] += noise[0];
	         Yk[i][1] += noise[1];
   	   }

/*============================================================================*/
/* R e c e i v e r:    A d a p t i v e     R L S    E Q                       */
/*============================================================================*/
         // Training Mode for adaptive RLS EQ and CE
         MMSE_EQ_RLS(Yk, train_sym, EQ_out, Num_tap, delay, lamda, p);
         //MMSE_EQ_RLS(Yk, train_sym, Num_tap, delay, lamda, p);

         // Soft Decision and Hard Decision of 4-QAM
         for(i=0; i<Num_slot; i++)
         	for(j=0; j<Slot_length; j++)
            {
            	//Ak[2*(j+i*Slot_length)] = tanh(0.5*(sqrt(2.0)*EQ_out[j+i*(Slot_length+Training)+Training][0]/noise_pwr));
               //Ak[2*(j+i*Slot_length)+1] = tanh(0.5*(sqrt(2.0)*EQ_out[j+i*(Slot_length+Training)+Training][1]/noise_pwr));
               Ak[2*(j+i*Slot_length)] = tanh(0.5*(sqrt(2.0)*EQ_out[j+i*(Slot_length+Training)+Training][0]));
               Ak[2*(j+i*Slot_length)+1] = tanh(0.5*(sqrt(2.0)*EQ_out[j+i*(Slot_length+Training)+Training][1]));
               //Ak[2*(j+i*Slot_length)] = sqrt(2.0) * EQ_out[j+i*(Slot_length+Training)+Training][0];
               //Ak[2*(j+i*Slot_length)+1] = sqrt(2.0) * EQ_out[j+i*(Slot_length+Training)+Training][1];

               if(EQ_out[j+i*(Slot_length+Training)+Training][0] >= 0.0)
               	Hk[2*(j+i*Slot_length)] = 1;
               else
               	Hk[2*(j+i*Slot_length)] = 0;

               if(EQ_out[j+i*(Slot_length+Training)+Training][1] >= 0.0)
               	Hk[2*(j+i*Slot_length)+1] = 1;
               else
               	Hk[2*(j+i*Slot_length)+1] = 0;
            }

         for(i=0; i<Num_slot*Slot_length*2; i++)
         	err_count += Error_count(data_bit[i], Hk[i]);

         // Copy the 2nd packet as the 1st packet and record inter-packet interference
         for(i=0; i<Num_slot*Training*2; i++)
         	train_sym[i] = train_sym[i+Num_slot*Training*2];

         for(i=0; i<Num_slot*Slot_length*2; i++)
         	data_bit[i] = data_bit[i+Num_slot*Slot_length*2];

         for(j=0; j<2; j++)
         	for(i=0; i<Num_slot*(Slot_length+Training); i++)
            	Yk[i][j] = Yk[i+Num_slot*(Slot_length+Training)][j];

         for(i=0; i<Num_slot*(Slot_length+Training); i++)
         	for(l=0; l<Num_path; l++)
            	for(x=0; x<2; x++)
               	fade_pattern[l][i][x] = fade_pattern[l][i+Num_slot*(Slot_length+Training)][x];

         for(j=0; j<Num_path; j++)
      		for(i=0; i<Num_slot*(Slot_length+Training)+j; i++)
	         {
   	         ch_matrix_I[j][i] = ch_matrix_I[j][i+Num_slot*(Slot_length+Training)];
      	      ch_matrix_Q[j][i] = ch_matrix_Q[j][i+Num_slot*(Slot_length+Training)];
         	}

         p++;
		} while(err_count <= 10000 && p < num_packet);

      // Statistics and records
      cout << "Error Rate = ";
      err_rate = err_count / (double)(Num_slot*Slot_length*2*p);
      printf("%e, ", err_rate);
      cout << endl;

      fprintf(ber, "%f ", Eb_No);
      fprintf(records, "%f ", Eb_No);
      fprintf(ber, "%e ", err_rate);
      fprintf(records, "%e ", err_rate);
      fprintf(ber, "\n");
      fprintf(records, "\n");
      fflush(records);
      fflush(ber);
   }

   delete data_bit;
   delete Ak;
   delete Hk;
   for(i=0; i<Num_slot*(Slot_length+Training); i++)
   	delete EQ_out[i];
   delete EQ_out;
   delete train_sym;
   for(i=0; i<Num_slot*(Slot_length+Training)*2; i++)
   	delete Yk[i];
   delete Yk;
   for(i=0; i<Num_path; i++)
   	for(l=0; l<2*Num_slot*(Slot_length+Training); l++)
   		delete fade_pattern[i][l];
   for(l=0; l<Num_path; l++)
   	delete fade_pattern[l];
	delete fade_pattern;
   delete sym_I;
   delete sym_Q;
   for(i=0; i<Num_path; i++)
   	delete ch_matrix_I[i];
   delete ch_matrix_I;
   for(i=0; i<Num_path; i++)
   	delete ch_matrix_Q[i];
   delete ch_matrix_Q;
   for(i=0; i<Num_tap; i++)
   	delete W[i];
   delete W;
   for(i=0; i<Num_tap; i++)
   	delete U[i];
   delete U;
   for(i=0; i<Num_tap; i++)
   	for(j=0; j<Num_tap; j++)
      	delete P[i][j];
   for(i=0; i<Num_tap; i++)
   	delete P[i];
   delete P;
   for(i=0; i<Num_tap; i++)
   	delete Kalman[i];
   delete Kalman;

   end = time(NULL);
   printf("Total elapsed time: %.0f(sec)\n", difftime(end,start));
   fprintf(records, "\nTotal elapsed time: %.0f(sec)\n\n", difftime(end,start));
   fclose(records);
   fclose(ber);
   printf("This program is ended. Press any key to continue.\n");
   getch();

   return 0;
}

void MMSE_EQ_RLS(double **Received, int *train_sym, double **EQ_out, int Num_tap, int delay, double lamda, int initial)
{
/*===============================================================*/
/* Adaptive MMSE EQ (Kalman Filter, RLS direct form) for 4-QAM   */
/*===============================================================*/
   int i, j, k, l, detected[2];
   double sym_I, sym_Q, e[2], Out[2], **pi, Denominator[2];

   pi = new double*[Num_tap];
   for(i=0; i<Num_tap; i++)
   	pi[i] = new double[2];

   // Input signal of adaptive RLS EQ
   if(initial == 0)
   {
   	//for(j=0; j<delay-1; j++)
      for(j=0; j<delay; j++)
	   {
   		for(l=Num_tap-1; l>0; l--)
      	{
      		U[l][0] = U[l-1][0];
	         U[l][1] = U[l-1][1];
   	   }
      	U[0][0] = Received[j][0];
	      U[0][1] = Received[j][1];
   	}
   }

   for(i=0; i<Num_slot; i++)
   {
      for(k=0; k<Training; k++)
      {
	   	// Input signal of adaptive RLS EQ
   		for(l=Num_tap-1; l>0; l--)
      	{
      		U[l][0] = U[l-1][0];
	         U[l][1] = U[l-1][1];
   	   }
      	//U[0][0] = Received[k+i*(Training+Slot_length)+delay-1][0];
	      //U[0][1] = Received[k+i*(Training+Slot_length)+delay-1][1];
         U[0][0] = Received[k+i*(Training+Slot_length)+delay][0];
	      U[0][1] = Received[k+i*(Training+Slot_length)+delay][1];

      	// Compute output
	      Out[0] = Out[1] = 0.0;
   	   for(j=0; j<Num_tap; j++)
      	{
      		Out[0] += (W[j][0]*U[j][0] + W[j][1]*U[j][1]);
	         Out[1] += (W[j][0]*U[j][1] - W[j][1]*U[j][0]);
   	   }

	      // Compute Kalman gain vector
   	   for(j=0; j<Num_tap; j++)
      		pi[j][0] = pi[j][1] = 0.0;

	      for(j=0; j<Num_tap; j++) 		// row index
   	   	for(l=0; l<Num_tap; l++)   // column index
      	   {
         		pi[j][0] += (P[j][l][0]*U[l][0] - P[j][l][1]*U[l][1]);
            	pi[j][1] += (P[j][l][1]*U[l][0] + P[j][l][0]*U[l][1]);
	         }

      	Denominator[0] = Denominator[1] = 0.0;
	      for(j=0; j<Num_tap; j++)
   	   {
      		Denominator[0] += (U[j][0]*pi[j][0] + U[j][1]*pi[j][1]);
         	Denominator[1] += (U[j][0]*pi[j][1] - U[j][1]*pi[j][0]);
	      }

   	   for(j=0; j<Num_tap; j++)
      	{
	      	Kalman[j][0] = (pi[j][0]*(Denominator[0]+lamda)+pi[j][1]*Denominator[1])
   	      					/(pow(Denominator[0]+lamda,2)+pow(Denominator[1],2));
      	   Kalman[j][1] = (pi[j][1]*(Denominator[0]+lamda)-pi[j][0]*Denominator[1])
         						/(pow(Denominator[0]+lamda,2)+pow(Denominator[1],2));
	      }

	      // Compute error signal (4-QAM)
      	// 4-QAM Mapping
         sym_I = (2*train_sym[2*(k+i*Training)]-1)/sqrt(2.0);
         sym_Q = (2*train_sym[2*(k+i*Training)+1]-1)/sqrt(2.0);

         // Error Signal
      	e[0] = sym_I - Out[0];
      	e[1] = sym_Q - Out[1];

         EQ_out[k+i*(Training+Slot_length)][0] = Out[0];
         EQ_out[k+i*(Training+Slot_length)][1] = Out[1];

         // Update Coefficients
      	for(j=0; j<Num_tap; j++)
	      {
   	   	W[j][0] += (Kalman[j][0]*e[0] + Kalman[j][1]*e[1]);
      	   W[j][1] += (Kalman[j][1]*e[0] - Kalman[j][0]*e[1]);
	      }

   	   // Update inverse of the correlation matrix
      	for(j=0; j<Num_tap; j++)
      		pi[j][0] = pi[j][1] = 0.0;

	      for(l=0; l<Num_tap; l++)      // column index
   	   	for(j=0; j<Num_tap; j++)	// row index
      	   {
         		pi[l][0] += (U[j][0]*P[j][l][0] + U[j][1]*P[j][l][1]);
            	pi[l][1] += (U[j][0]*P[j][l][1] - U[j][1]*P[j][l][0]);
	         }

   	   for(j=0; j<Num_tap; j++) 		// row index
      		for(l=0; l<Num_tap; l++)   // column index
         	{
         		P[j][l][0] = (P[j][l][0] - (Kalman[j][0]*pi[l][0] - Kalman[j][1]*pi[l][1]))/lamda;
	            P[j][l][1] = (P[j][l][1] - (Kalman[j][0]*pi[l][1] + Kalman[j][1]*pi[l][0]))/lamda;
   	      }
      }

      for(k=0; k<Slot_length; k++)
      {
         // Input signal of adaptive RLS EQ
   		for(l=Num_tap-1; l>0; l--)
      	{
      		U[l][0] = U[l-1][0];
	         U[l][1] = U[l-1][1];
   	   }
      	//U[0][0] = Received[k+i*(Training+Slot_length)+Training+delay-1][0];
	      //U[0][1] = Received[k+i*(Training+Slot_length)+Training+delay-1][1];
         U[0][0] = Received[k+i*(Training+Slot_length)+Training+delay][0];
	      U[0][1] = Received[k+i*(Training+Slot_length)+Training+delay][1];

         // Compute output
	      Out[0] = Out[1] = 0.0;
   	   for(j=0; j<Num_tap; j++)
      	{
      		Out[0] += (W[j][0]*U[j][0] + W[j][1]*U[j][1]);
	         Out[1] += (W[j][0]*U[j][1] - W[j][1]*U[j][0]);
   	   }

	      // Compute Kalman gain vector
   	   for(j=0; j<Num_tap; j++)
      		pi[j][0] = pi[j][1] = 0.0;

	      for(j=0; j<Num_tap; j++) 		// row index
   	   	for(l=0; l<Num_tap; l++)   // column index
      	   {
         		pi[j][0] += (P[j][l][0]*U[l][0] - P[j][l][1]*U[l][1]);
            	pi[j][1] += (P[j][l][1]*U[l][0] + P[j][l][0]*U[l][1]);
	         }

      	Denominator[0] = Denominator[1] = 0.0;
	      for(j=0; j<Num_tap; j++)
   	   {
      		Denominator[0] += (U[j][0]*pi[j][0] + U[j][1]*pi[j][1]);
         	Denominator[1] += (U[j][0]*pi[j][1] - U[j][1]*pi[j][0]);
	      }

   	   for(j=0; j<Num_tap; j++)
      	{
	      	Kalman[j][0] = (pi[j][0]*(Denominator[0]+lamda)+pi[j][1]*Denominator[1])
   	      					/(pow(Denominator[0]+lamda,2)+pow(Denominator[1],2));
      	   Kalman[j][1] = (pi[j][1]*(Denominator[0]+lamda)-pi[j][0]*Denominator[1])
         						/(pow(Denominator[0]+lamda,2)+pow(Denominator[1],2));
	      }

         // 4-QAM De-mapping (Soft Decision)
         sym_I = tanh(0.5*(sqrt(2.0)*Out[0]/noise_pwr))/sqrt(2.0);
         sym_Q = tanh(0.5*(sqrt(2.0)*Out[1]/noise_pwr))/sqrt(2.0);

      	// 4-QAM De-mapping (Hard Decision)
         if(Out[0] >= 0.0)
      		detected[0] = 1;
         else
	         detected[0] = 0;

         if(Out[1] >= 0.0)
      		detected[1] = 1;
         else
	         detected[1] = 0;
/*
         // 4-QAM Mapping
         sym_I = (2*detected[0]-1)/sqrt(2.0);
         sym_Q = (2*detected[1]-1)/sqrt(2.0);
  */
         // Compute error signal (4-QAM)
      	e[0] = sym_I - Out[0];
      	e[1] = sym_Q - Out[1];

         EQ_out[k+i*(Training+Slot_length)+Training][0] = Out[0];
         EQ_out[k+i*(Training+Slot_length)+Training][1] = Out[1];

         // Update Coefficients
      	for(j=0; j<Num_tap; j++)
	      {
   	   	W[j][0] += (Kalman[j][0]*e[0] + Kalman[j][1]*e[1]);
      	   W[j][1] += (Kalman[j][1]*e[0] - Kalman[j][0]*e[1]);
	      }

   	   // Update inverse of the correlation matrix
      	for(j=0; j<Num_tap; j++)
      		pi[j][0] = pi[j][1] = 0.0;

	      for(l=0; l<Num_tap; l++)      // column index
   	   	for(j=0; j<Num_tap; j++)	// row index
      	   {
         		pi[l][0] += (U[j][0]*P[j][l][0] + U[j][1]*P[j][l][1]);
            	pi[l][1] += (U[j][0]*P[j][l][1] - U[j][1]*P[j][l][0]);
	         }

   	   for(j=0; j<Num_tap; j++) 		// row index
      		for(l=0; l<Num_tap; l++)   // column index
         	{
         		P[j][l][0] = (P[j][l][0] - (Kalman[j][0]*pi[l][0] - Kalman[j][1]*pi[l][1]))/lamda;
	            P[j][l][1] = (P[j][l][1] - (Kalman[j][0]*pi[l][1] + Kalman[j][1]*pi[l][0]))/lamda;
   	      }
      }
   }

   for(i=0; i<Num_tap; i++)
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

