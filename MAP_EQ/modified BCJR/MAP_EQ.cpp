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
void MAP_EQ_BCJR(double *, double *, double, double *, int, int, int);

#define Pi 	3.14159265358979
#define n 	1											// regard channel as rate-1 convolutional code
#define k	1
#define m	2											// channel memory
#define N	1000	 									// packet length
//#define K	512
#define num_packet 10000								// number of packets simulated
const int num_state = 4;							// number of channel states
const int num_in_sym = 2;							// number of input symbols (BPSK)
const float CIR[3] = {0.407, 0.815, 0.407};	// Channel Impulse Response

struct ch_trel_diag			// Data structure of trellis diagram for each branch
		{
      	int from;		// from_state of trellis diagram
         int to;			// to_state of trellis diagram
         int in;			// input data symbol of trellis diagram
         float out;		// output symbol of trellis diagram
      };
struct ch_trel_diag ch_Trellis[4][2];		// Trellis[num_state][num_in_sym]

int main(void)
{
	time_t  t, start, end;
	int i, p, input, *data_bit, *ch;
   int from_state, to_state, tmp_state, *Ak, err_count;
   double snr, Eb_No, noise_pwr, noise[2], *Yk, err_rate;
   double out_sym, *LLR, *intrinsic;
   FILE *trelis,*ber, *records;

   start = time(NULL);
   printf("BER Performance of MAP Equalizer in ISI Channel\n");
	printf("Channel impulse response is {0.407,0.815,0.407}\n");
   printf("This program is running. Don't close, please!\n\n");

   data_bit = new int[N];
   ch = new int[m+1];
   Ak = new int[N];
   Yk = new double[N];
   LLR = new double[N];
   intrinsic = new double[N];

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
         ch_Trellis[from_state][input].from = from_state;
         ch_Trellis[from_state][input].to = to_state;
         ch_Trellis[from_state][input].in = (2*input-1);
         ch_Trellis[from_state][input].out = out_sym;
         fprintf(trelis, "%4d %4d %5d %8.3f\n", from_state, to_state, (2*input-1), out_sym);
      }
   }
   fclose(trelis);

   for(i=0; i<N; i++)
   	intrinsic[i] = 0.0;
      
/************************/
/* main simulation loop */
/************************/
   for(snr=30; snr<=30; snr++)
   {
   	err_count = 0;
   	// noise power calculation
      Eb_No = (double)snr;
      noise_pwr = 1.0/pow(10.0, Eb_No/10.0);	// BPSK, Nyquist filter assumption
      printf("Eb_No = %f, ", Eb_No);
      from_state = random(4);
      for(p=0; p<num_packet; p++)
      {
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

         MAP_EQ_BCJR(Yk, intrinsic, noise_pwr, LLR, N, num_state, num_in_sym);

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
   delete LLR;
   delete intrinsic;

   records = fopen("record_MAP_EQ.log", "a");
   end = time(NULL);
   printf("Total elapsed time: %.0f(sec)\n", difftime(end,start));
   fprintf(records, "Total elapsed time: %.0f(sec)\n\n", difftime(end,start));
   fclose(records);
   printf("This program is ended. Press any key to continue.\n");
   getch();

   return 0;
}

void MAP_EQ_BCJR(double *Data_in, double *intrinsic, double noise_pwr,
					  double *LLR, int Packet_length, int Num_state, int Num_in)
{
/******************************************************************/
/* MAP Equalizer (modified BCJR Algorithm) for static ISI channel */
/******************************************************************/
   const float pi = 3.14159265358979;
	int i, j, l;
   double p1, *normal, **alpha, **beta, ***gamma, *delta, min, **a_priori;

   a_priori = new double*[Packet_length];			// a-priori probability, a-priori[time index][input]
   for(i=0; i<Packet_length; i++)
   	a_priori[i] = new double[Num_in];
   alpha = new double*[Packet_length+1];			// alpha[time index][state]
   for(i=0; i<=Packet_length; i++)
      alpha[i] = new double[Num_state];
   beta = new double*[Packet_length+1];       	// beta[time index][state]
   for(i=0; i<=Packet_length; i++)
      beta[i] = new double[Num_state];
   gamma = new double**[Packet_length];			// gamma[time index][state][input]
   for(i=0; i<Packet_length; i++)
      gamma[i] = new double*[Num_state];
   for(i=0; i<Packet_length; i++)
   	for(j=0; j<Num_state; j++)
         gamma[i][j] = new double[Num_in];
   normal = new double[Packet_length+1];			// renormalization of BCJR
   delta = new double[Num_in];

	// Initialization of alpha and beta
   for(i=1; i<=Packet_length; i++)					// alpha[time index][state]
   	for(j=0; j<Num_state; j++)
      	alpha[i][j] = 0.0;
   for(j=0; j<Num_state; j++)
   	alpha[0][j] = 1/(float)Num_state;

   for(i=0; i<Packet_length; i++)           		// beta[time index][state]
   	for(j=0; j<Num_state; j++)
      	beta[i][j] = 0.0;
   for(j=0; j<Num_state; j++)
   	beta[Packet_length][j] = 1.0;

   // Calculate a-priori probability from intrinsic information
   for(i=0; i<Packet_length; i++)					// time index
   	for(l=0; l<Num_in; l++)							// input symbol
      	a_priori[i][l] = exp(l*intrinsic[i]) / (1 + exp(intrinsic[i]));

   // calculate gamma[time index][state][input]
   for(i=0; i<Packet_length; i++)					// time index
   	for(j=0; j<Num_state; j++)						// state index
      	for(l=0; l<Num_in; l++)						// input symbol
         {
         	p1 = exp(-pow(Data_in[i]-ch_Trellis[j][l].out,2)/(2*noise_pwr))/sqrt(2*pi*noise_pwr);
            gamma[i][j][l] = a_priori[i][l] * p1;		// gamma[time index][state][input]
         }

   // calculate alpha[time index][state]
   for(i=1; i<=Packet_length; i++)					// time index
   {
   	for(j=0; j<Num_state; j++)						// from_state index
      	for(l=0; l<Num_in; l++)						// input bit
         	alpha[i][ch_Trellis[j][l].to] += alpha[i-1][j] * gamma[i-1][j][l];

      normal[i] = 0.0;									// for renormalization
      for(j=0; j<Num_state; j++)						// to_state index
      	normal[i] += alpha[i][j];

      for(j=0; j<Num_state; j++)
      	alpha[i][j] = alpha[i][j] / normal[i];
   }
   normal[0] = 1.0;

   // calculate beta[time index][state]
   for(i=Packet_length-1; i>0; i--)					// time index
   {
   	for(j=0; j<Num_state; j++)						// from_state index
      	for(l=0; l<Num_in; l++)						// input bit
         	beta[i][j] += beta[i+1][ch_Trellis[j][l].to] * gamma[i][j][l];

      for(j=0; j<Num_state; j++)
      	beta[i][j] = beta[i][j] / normal[i];
   }

   // calculate conditional LLR
   for(i=0; i<Packet_length; i++)					// time index
   {
   	min = 0.0;		// find the minimum product of alpha*gamma*beta
      for(j=0; j<Num_state; j++)
      	for(l=0; l<Num_in; l++)
         {
         	delta[0] = alpha[i][j] * gamma[i][j][l] * beta[i+1][ch_Trellis[j][l].to];

            if((delta[0] < min && delta[0] != 0.0) || min == 0.0)
            	min = delta[0];
         }

      if(min == 0.0 || min > 1.0)	// if all else fails, make min real small
      	min = 1E-100;

      delta[0] = delta[1] = 0.0;
      for(j=0; j<Num_state; j++)		// from_state index
      	for(l=0; l<Num_in; l++)	// input bit
         	delta[l] += alpha[i][j] * gamma[i][j][l] * beta[i+1][ch_Trellis[j][l].to];

      if(delta[1] == 0.0)
      	delta[1] = min;
      if(delta[0] == 0.0)
      	delta[0] = min;

      LLR[i] = log(delta[1]/delta[0]);
   }

   for(i=0; i<Packet_length; i++)
   	delete a_priori[i];
   delete a_priori;
   for(i=0; i<=Packet_length; i++)
      delete alpha[i];
	delete alpha;
   for(i=0; i<=Packet_length; i++)
       delete beta[i];
   delete beta;
   for(i=0; i<Packet_length; i++)
   	for(j=0; j<Num_state; j++)
         delete gamma[i][j];
   for(i=0; i<Packet_length; i++)
      delete gamma[i];
   delete gamma;
   delete normal;
   delete delta;
}
/******************************************************************************/

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

