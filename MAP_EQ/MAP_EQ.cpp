/***************************************************/
/* Author: Chao-wang Huang                         */
/* Date: Monday, September 27, 2004                */
/* An MAP Equalizer is simulated                   */
/* Channel impulse response: {0.407, 0.815, 0.407} */
/***************************************************/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <conio.h>
#include <iostream.h>

void AWGN_noise(float, double, double *);
int Error_count(int, int);

#define Pi 	3.14159265358979
#define n 	1											// regard channel as rate-1 convolutional code
#define k	1
#define m	2											// channel memory
#define N	1024	 									// packet length
#define K	512
#define num_packet 1000								// number of packets simulated
const int num_state = 4;							// number of channel states
const int num_in_sym = 2;							// number of input symbols (BPSK)
const float CIR[3] = {0.407, 0.815, 0.407};	// Channel Impulse Response

int main(void)
{
	time_t  t, start, end;
	int i, j, l, p, input, *data_bit, *ch;
   int from_state, to_state, tmp_state, *Ak, err_count;
   float out_sym;
   double a_priori=0.5, snr, Eb_No, noise_pwr, noise[2], *Yk, err_rate;
   double p1, sum, **alpha, **beta, ***gamma, *delta, *LLR, min;
   FILE *trelis,*ber, *records;
   struct trel_diag			// Data structure of trellis diagram for each branch
   		{
         	int from;		// from_state of trellis diagram
            int to;			// to_state of trellis diagram
            int in;			// input data symbol of trellis diagram
            float out;		// output symbol of trellis diagram
         };
   struct trel_diag Trellis[4][2];		// Trellis[num_state][num_in_sym]

   start = time(NULL);
   printf("BER Performance of MAP Equalizer in ISI Channel\n");
	printf("Channel impulse response is {0.407,0.815,0.407}\n");
   printf("This program is running. Don't close, please!\n\n");

   data_bit = new int[N];
   ch = new int[m+1];
   Ak = new int[N];
   Yk = new double[N];
   delta = new double[num_in_sym];
   LLR = new double[N];
   alpha = new double*[N+1];			  		// alpha[time index][state]
   for(i=0; i<=N; i++)
      alpha[i] = new double[num_state];
   beta = new double*[N+1];       			// beta[time index][state]
   for(i=0; i<=N; i++)
      beta[i] = new double[num_state];
   gamma = new double**[N];					// gamma[time index][state][input]
   for(i=0; i<N; i++)
      gamma[i] = new double*[num_state];
   for(i=0; i<N; i++)
   	for(j=0; j<num_state; j++)
         gamma[i][j] = new double[num_in_sym];

	srand((unsigned) time(&t));

/****************************************************************/
/* Generate TrellisDiagram for {0.407,0.815,0.407} ISI channel */
/****************************************************************/
	trelis = fopen("Trellis_diagram.log", "w");
   fprintf(trelis, "Trellis diagram of (%.3f,%.3f,%.3f) ISI channel\n", CIR[0], CIR[1], CIR[2]);
   fprintf(trelis, "s(k-1) s(k) input output\n");
   for(from_state=0; from_state<num_state; from_state++) // from_state of trellis diagram
   {
   	for(input=0; input<num_in_sym; input++)		// input symbol (2*input-1) of trellis diagram
      {
      	tmp_state = from_state;
         tmp_state = (tmp_state << 1) ^ (input & 0x01);  // read input bit
         ch[0] = 2*(tmp_state & 0x01) - 1;
         ch[1] = 2*((tmp_state >> 1) & 0x01) - 1;
         ch[2] = 2*((tmp_state >> 2) & 0x01) - 1;
         out_sym = ch[0]*CIR[0] + ch[1]*CIR[1] + ch[2]*CIR[2];		// Calculate output symbol
         to_state = tmp_state & (num_state-1); 					// to_state of trellis diagram
         Trellis[from_state][input].from = from_state;
         Trellis[from_state][input].to = to_state;
         Trellis[from_state][input].in = (2*input-1);
         Trellis[from_state][input].out = out_sym;
         fprintf(trelis, "%4d %4d %5d %8.3f\n", from_state, to_state, (2*input-1), out_sym);
      }
   }
   fclose(trelis);

/************************/
/* main simulation loop */
/************************/
   for(snr=0; snr<=10; snr+=2)
   {
   	err_count = 0;
   	// noise power calculation
      Eb_No = (double)snr;
      noise_pwr = 1.0/pow(10.0, Eb_No/10.0);	// BPSK, Nyquist filter assumption
      printf("Eb_No = %f, ", Eb_No);
      from_state = random(4);
      for(p=0; p<num_packet; p++)
      {
/*   		for(i=0; i<K-2; i++)		// Generate random information bit stream
   		{
   			if(rand()/RAND_MAX>=0.5)
      			data_bit[i] = 1;
      		else
      			data_bit[i] = 0;
   		}
*/

/*****************************************************/
/* BPSK mapping and ISI Channel: {0.407,0.815,0.407} */
/*****************************************************/
   		for(i=0; i<N; i++)
   		{
         	data_bit[i] = random(2);		// Generate random information bit stream
         	tmp_state = from_state;
            tmp_state = (tmp_state << 1) ^ (data_bit[i] & 0x01);  // read input bit
            ch[0] = 2*(tmp_state & 0x01) - 1;					// input symbol
            ch[1] = 2*((tmp_state >> 1) & 0x01) - 1;        // channel memory
            ch[2] = 2*((tmp_state >> 2) & 0x01) - 1;        // channel memory
            out_sym = ch[0]*CIR[0] + ch[1]*CIR[1] + ch[2]*CIR[2];	// Calculate output symbol
            from_state = tmp_state & (num_state-1); 			// to_state of trellis diagram

            /* AWGN channel */
            AWGN_noise(0, noise_pwr, &noise[0]);
         	Yk[i] = out_sym + noise[0];
   		}

/*********************************************/
/* MAP Equalizer (normalized BCJR Algorithm) */
/*********************************************/
			// Initialization of alpha and beta
         for(i=1; i<=N; i++)				// alpha[time index][state]
         	for(j=0; j<num_state; j++)
         		alpha[i][j] = 0.0;
         for(j=0; j<num_state; j++)
         	alpha[0][j] = 1/(float)num_state;

         for(i=0; i<N; i++)           // beta[time index][state]
         	for(j=0; j<num_state; j++)
         		beta[i][j] = 0.0;
         for(j=0; j<num_state; j++)
         	beta[N][j] = 1.0;

			// calculate gamma[time index][state][input]
         for(i=0; i<N; i++)	// time index
         	for(j=0; j<num_state; j++)		// state index
            	for(l=0; l<num_in_sym; l++)	// input symbol
               {
						p1 = exp(-pow(Yk[i]-Trellis[j][l].out,2)/(2*noise_pwr))/sqrt(2*Pi*noise_pwr);
         			gamma[i][j][l] = a_priori * p1;		// gamma[time index][state][input]
               }

         // calculate alpha[time index][state]
         for(i=1; i<=N; i++)		// time index
         {
         	for(j=0; j<num_state; j++)		// from_state index
            	for(l=0; l<num_in_sym; l++)	// input bit
         			alpha[i][Trellis[j][l].to] += alpha[i-1][j] * gamma[i-1][j][l];

            sum = 0.0;		// for renormalization
            for(j=0; j<num_state; j++)		// to_state index
            	sum += alpha[i][j];

            for(j=0; j<num_state; j++)
            	alpha[i][j] = alpha[i][j] / sum;
         }

         // calculate beta[time index][state]
         for(i=N-1; i>0; i--)		// time index
         {
         	for(j=0; j<num_state; j++)		// from_state index
            	for(l=0; l<num_in_sym; l++)	// input bit
         			beta[i][j] += beta[i+1][Trellis[j][l].to] * gamma[i][j][l];

            sum = 0.0;		// for renormalization
            for(j=0; j<num_state; j++)		// from_state index
            	sum += beta[i][j];

            for(j=0; j<num_state; j++)
            	beta[i][j] = beta[i][j] / sum;
         }

         // calculate conditional LLR
         for(i=0; i<N; i++)		// time index
         {
         	min = 0.0;		// find the minimum product of alpha*gamma*beta
            for(j=0; j<num_state; j++)
            	for(l=0; l<num_in_sym; l++)
               {
                  delta[0] = alpha[i][j] * gamma[i][j][l] * beta[i+1][Trellis[j][l].to];

                  if((delta[0] < min && delta[0] != 0.0) || min == 0.0)
                  	min = delta[0];
					}

            if(min == 0.0 || min > 1.0)	// if all else fails, make min real small
            	min = 1E-100;

         	delta[0] = delta[1] = 0.0;
         	for(j=0; j<num_state; j++)		// from_state index
            	for(l=0; l<num_in_sym; l++)	// input bit
               	delta[l] += alpha[i][j] * gamma[i][j][l] * beta[i+1][Trellis[j][l].to];

            if(delta[1] == 0.0)
            	delta[1] = min;
//            else if(delta[0] == 0.0)
            if(delta[0] == 0.0)
            	delta[0] = min;

            LLR[i] = log(delta[1]/delta[0]);
         }

/*******************************************/

			for(i=0; i<N; i++)	// data decision
         {
         	if(LLR[i]>=0)
            	Ak[i] = 1;
            else
            	Ak[i] = 0;

            err_count += Error_count(data_bit[i], Ak[i]);
         }
		}

      // Statistics and records
      err_rate = err_count / (double)(N*num_packet);
      printf("Error rate = %e\n", err_rate);
      ber=fopen("ber.log","a");
   	records = fopen("record_MAP_EQ.log", "a");
      fprintf(ber, "%f %e\n", Eb_No, err_rate);
      fprintf(records, "MAP Equalizer in ISI channel {0.407,0.815,0.407} with Eb/N0 = %f dB\n", Eb_No);
   	fprintf(records, "Average bit error rate = %e\n", err_rate);
   	fprintf(records, "number of bits of simulation = %d\n\n", N*num_packet);

      fclose(ber);
	   fclose(records);
   }

   delete data_bit;
   delete ch;
   delete Ak;
   delete Yk;
   delete delta;
   delete LLR;
   delete alpha;
   delete beta;
   delete gamma;

   records = fopen("record_MAP_EQ.log", "a");
   end = time(NULL);
   printf("Total elapsed time: %.0f(sec)\n", difftime(end,start));
   fprintf(records, "Total elapsed time: %.0f(sec)\n\n", difftime(end,start));
   fclose(records);
   printf("This program is ended. Press any key to continue.\n");
   getch();

   return 0;
}

void AWGN_noise(float mu, double variance, double *noise)
{
//	const  float Pi = 3.14159265358979;
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

