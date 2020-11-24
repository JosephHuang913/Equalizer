/***************************************************/
/* Author: Chao-wang Huang                         */
/* Date: Monday, January 3, 2005                   */
/* An MLSE Equalizer is simulated                  */
/* Channel impulse response: {0.407, 0.815, 0.407} */
/* Algorithm: Soft Output Viterbi Algorithm (SOVA) */
/***************************************************/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <conio.h>
#include <iostream.h>
#include <limits.h>
#include <float.h>

void AWGN_noise(float, double, double *);
int Error_count(int, int);
void MLSE_SOVA(double *, int *, double *, int, int, int);

#define Pi 	3.14159265358979
#define n	1					// (n,k,m) Convolutional Code, # of output symbols
#define k	1					// (n,k,m) Convolutional Code, # of input symbols
#define m	2					// (n,k,m) Convolutional Code, # of memories
#define K	1000		  		// packet length of information bit stream
#define N	K*(n/k)	 		// packet length of coded bit stream
#define num_packet 10000		// number of packets simulated
#define N_state  4
#define num_in_sym 2
//const int num_in_sym = pow(2,k);					// number of input symbols (BPSK)
//const int ch_num_state = 4;						// number of channel states
const float CIR[3] = {0.407, 0.815, 0.407};	// Channel Impulse Response

struct ch_trel_diag			// Data structure of trellis diagram for each branch
		{
      	int from;			// from_state of trellis diagram
         int to;				// to_state of trellis diagram
         int in;				// input data bit of trellis diagram
         float out;			// output codeword symbol of trellis diagram
      };
struct ch_trel_diag ch_Trellis[4][2];			// Trellis[num_state][num_in_sym]

int main(void)
{
	time_t  t, start, end;
	int i, p, input, *data_bit, *Hk;
   int from_state, to_state, tmp_state, err_count[2], *Ak, *ch;
   float out_sym;
   double *bpsk, snr, Eb_No, noise_pwr, noise[2], *Yk, err_rate[2];
   double *Dk, *LLR;
   FILE *trelis, *ber, *records;

   start = time(NULL);
   printf("BER Performance of MLSE in ISI Channel\n");
	printf("Channel impulse response is {0.407,0.815,0.407}\n");
   printf("\nThis program is running. Don't close, please!\n\n");

   records = fopen("Records_MLSE.log", "a");
   fprintf(records, "MLSE Equalizer in ISI channel {0.407,0.815,0.407}\n");
   fprintf(records, "number of bits of simulation = %d\n\n", N*num_packet);
   fprintf(records, "Eb/No     BER\n");
   fflush(records);

   data_bit = new int[K];
   ch = new int[m+1];
   Ak = new int[K];
   Hk = new int[K];
   bpsk = new double[N];
   Yk = new double[N];		// Received signal
   Dk = new double[N];			// Soft Decision
   LLR = new double[K];		// Log-likelihood Ratio

	srand((unsigned) time(&t));

/****************************************************************/
/* Generate TrellisDiagram for {0.407,0.815,0.407} ISI channel */
/****************************************************************/
	trelis = fopen("Trellis_diagram.log", "w");
   fprintf(trelis, "Trellis diagram of (%.3f,%.3f,%.3f) ISI channel\n", CIR[0], CIR[1], CIR[2]);
   fprintf(trelis, "s(k-1) s(k) input output\n");
   for(from_state=0; from_state<N_state; from_state++) // from_state of trellis diagram
   {
   	for(input=0; input<num_in_sym; input++)		// input symbol (2*input-1) of trellis diagram
      {
      	tmp_state = from_state;
         tmp_state = (tmp_state << 1) ^ (input & 0x01);  // read input bit
         ch[0] = 2*(tmp_state & 0x01) - 1;
         ch[1] = 2*((tmp_state >> 1) & 0x01) - 1;
         ch[2] = 2*((tmp_state >> 2) & 0x01) - 1;
         out_sym = ch[0]*CIR[0] + ch[1]*CIR[1] + ch[2]*CIR[2];		// Calculate output symbol
         to_state = tmp_state & (N_state-1); 					// to_state of trellis diagram
         ch_Trellis[from_state][input].from = from_state;
         ch_Trellis[from_state][input].to = to_state;
         ch_Trellis[from_state][input].in = (2*input-1);
         ch_Trellis[from_state][input].out = out_sym;
         fprintf(trelis, "%4d %4d %5d %8.3f\n", from_state, to_state, (2*input-1), out_sym);
      }
   }
   fclose(trelis);

/************************/
/* main simulation loop */
/************************/
   ber=fopen("ber.log","w");
   for(snr=0; snr<=10; snr++)
   {
   	err_count[0] = err_count[1] = 0;
   	// noise power calculation
      Eb_No = (double)snr;
      noise_pwr = 1.0/pow(10.0, Eb_No/10.0);	// BPSK, Nyquist filter assumption
      printf("Eb_No = %f, ", Eb_No);

      from_state = random(4);						// Initial State of Channel
      for(p=0; p<num_packet; p++)
      {
/*****************************************************/
/* BPSK mapping and ISI Channel: {0.407,0.815,0.407} */
/*****************************************************/
   		for(i=0; i<K; i++)
   		{
         	data_bit[i] = random(2);		// Generate random information bit stream
         	tmp_state = from_state;
            tmp_state = (tmp_state << 1) ^ (data_bit[i] & 0x01);  // read input bit
            ch[0] = 2*(tmp_state & 0x01) - 1;					// input symbol
            ch[1] = 2*((tmp_state >> 1) & 0x01) - 1;        // channel memory
            ch[2] = 2*((tmp_state >> 2) & 0x01) - 1;        // channel memory
            out_sym = ch[0]*CIR[0] + ch[1]*CIR[1] + ch[2]*CIR[2];	// Calculate output symbol
            from_state = tmp_state & (N_state-1); 			// to_state of trellis diagram

            /* AWGN channel */
            AWGN_noise(0, noise_pwr, &noise[0]);
         	Yk[i] = out_sym + noise[0];
            // Soft Decision
//            Dk[i] = tanh(Yk[i]/noise_pwr);
            Dk[i] = Yk[i];
         }

         MLSE_SOVA(Dk, Hk, LLR, K, N_state, num_in_sym);

         // Hard Decision
         for(i=0; i<K; i++)
         {
         	if(LLR[i]>=0)
            	Ak[i] = 1;
            else
            	Ak[i] = 0;

			  	// Bit error count
            err_count[0] += Error_count(data_bit[i], Ak[i]);
            err_count[1] += Error_count(data_bit[i], Hk[i]);
         }
		}

      // Statistics and records
      err_rate[0] = err_count[0] / (double)(K*num_packet);
      err_rate[1] = err_count[1] / (double)(K*num_packet);
      printf("Error rate = %e %e\n", err_rate[0], err_rate[1]);
      fprintf(ber, "%f %e %e\n", Eb_No, err_rate[0], err_rate[1]);
      fprintf(records, "%f %e %e\n", Eb_No, err_rate[0], err_rate[1]);
      fflush(records);
   }

   fclose(ber);
   delete data_bit;
   delete ch;
   delete Ak;
   delete Hk;
   delete bpsk;
   delete Yk;
   delete Dk;
   delete LLR;

   end = time(NULL);
   printf("Total elapsed time: %.0f(sec)\n", difftime(end,start));
   fprintf(records, "\nTotal elapsed time: %.0f(sec)\n\n", difftime(end,start));
   fclose(records);
   printf("This program is ended. Press any key to continue.\n");
   getch();

   return 0;
}

void MLSE_SOVA(double *Data_in, int *Data_out, double *Soft_out, int Packet_length, int Num_state, int Num_in)
{
/*********************************************************************/
/* MLSE: Soft Output Viterbi Algorithm (SOVA) for static ISI channel */
/*********************************************************************/
	int i, j, l, q, survivor_state;
   double **mju_f, *mju, *survival_metric, metric, ***branch_metric, **mju_b;
   double mju_tmp, survivor_metric;

   struct surv        			// Data structure of survival path
   		{
         	double metric;		// Path metric
            int data_in[K];	// input bit stream, K: packet length, global
            int state[K];		// state transition sequence, K: packet length, global
         };
   struct surv survival_path[N_state], survival_temp[N_state];
   									// Survival_path[num_state], N_state: number of states of channel, globle

   survival_metric = new double[Num_state];
   mju = new double[Num_in];							// minimum path metric
   mju_f = new double*[Packet_length+1];			// forward path-metric[time index][state]
   for(i=0; i<=Packet_length; i++)
      mju_f[i] = new double[Num_state];
   mju_b = new double*[Packet_length+1];       	// backward path-metric[time index][state]
   for(i=0; i<=Packet_length; i++)
      mju_b[i] = new double[Num_state];
   branch_metric = new double**[Packet_length];	// branch[time index][state][input]
   for(i=0; i<Packet_length; i++)
      branch_metric[i] = new double*[Num_state];
   for(i=0; i<Packet_length; i++)
   	for(j=0; j<Num_state; j++)
         branch_metric[i][j] = new double[Num_in];

   // Initialize survival path
   for(i=0; i<Num_state; i++)
   {
   	survival_path[i].metric = 0.0;
      mju_f[0][i] = 0.0;
   }

/*********************/
/* Forward Recursion */
/*********************/
	for(i=0; i<Packet_length; i++)
   {
   	for(j=0; j<Num_state; j++)					// Initialize the survival path metric
      	survival_metric[j] = DBL_MAX;

      for(j=0; j<Num_state; j++)					// from_state index
      	for(l=0; l<Num_in; l++)					// input bit
         {
         	// branch metric, Euclidean Distance
            branch_metric[i][j][l] = pow(Data_in[i]-ch_Trellis[j][l].out,2);
            metric = survival_path[j].metric + branch_metric[i][j][l];

            // find the survival path metric
            if(metric < survival_metric[ch_Trellis[j][l].to])
            {
            	survival_metric[ch_Trellis[j][l].to] = metric;

               // Record and refresh the survival path
               for(q=0; q<i; q++)
               {
               	survival_temp[ch_Trellis[j][l].to].data_in[q] = survival_path[j].data_in[q];
                  survival_temp[ch_Trellis[j][l].to].state[q] = survival_path[j].state[q];
               }
               survival_temp[ch_Trellis[j][l].to].data_in[i] = l;
               survival_temp[ch_Trellis[j][l].to].state[i] = ch_Trellis[j][l].to;
            }
         }

      // Record and refresh the survival path
      for(j=0; j<Num_state; j++)					// to_state index
      {
      	survival_path[j].metric = survival_metric[j];
         mju_f[i+1][j] = survival_metric[j];
         for(q=0; q<=i; q++)
         {
         	survival_path[j].data_in[q] = survival_temp[j].data_in[q];
            survival_path[j].state[q] = survival_temp[j].state[q];
         }
      }
   }

   // Find the path with the smallest path metric
   survivor_metric = survival_path[0].metric;
   survivor_state = 0;
   for(j=1; j<Num_state; j++)
   	if(survivor_metric > survival_path[j].metric)
      {
      	survivor_metric = survival_path[j].metric;	// survivor path metric
         survivor_state = j;									// survivor state
      }

/****************************************/
/* Backward Recursion and Soft Decision */
/****************************************/
	// Initialize survival path
   for(j=0; j<Num_state; j++)
   	//mju_b[Packet_length][j] = DBL_MAX;
      mju_b[Packet_length][j] = 0.0;
   //mju_b[Packet_length][survivor_state] = 0.0;

   for(i=Packet_length-1; i>=0; i--)
   {
   	for(j=0; j<Num_state; j++)					// Initialize the survival path metric
      	survival_metric[j] = DBL_MAX;

      for(j=0; j<Num_state; j++)					// from_state index
      	for(l=0; l<Num_in; l++)					// input bit
         {
         	metric = mju_b[i+1][ch_Trellis[j][l].to] + branch_metric[i][j][l];

            // find the survival path metric
            if(metric < survival_metric[j])
            	survival_metric[j] = metric;
         }

      // Record the survival path metric
      for(j=0; j<Num_state; j++)					// from_state index
      	mju_b[i][j] = survival_metric[j];

      // LLR Calculation
      mju[survival_path[survivor_state].data_in[i]] = survivor_metric;	// mju_f[Packet_length][survivor_state];

      mju[(survival_path[survivor_state].data_in[i]+1)%2] = DBL_MAX;
      for(j=0; j<Num_state; j++)					// from_state index
      {
      	mju_tmp = mju_f[i][j] + branch_metric[i][j][(survival_path[survivor_state].data_in[i]+1)%2]
      				 + mju_b[i+1][ch_Trellis[j][(survival_path[survivor_state].data_in[i]+1)%2].to];

         if(mju_tmp < mju[(survival_path[survivor_state].data_in[i]+1)%2])
         	mju[(survival_path[survivor_state].data_in[i]+1)%2] = mju_tmp;
      }

      Soft_out[i] = mju[0] - mju[1];
      Data_out[i] = survival_path[survivor_state].data_in[i];
   }

	delete survival_metric;
   delete mju;
   for(i=0; i<=Packet_length; i++)
      delete mju_f[i];
   delete mju_f;
   for(i=0; i<=Packet_length; i++)
       delete mju_b[i];
   delete mju_b;
   for(i=0; i<Packet_length; i++)
   	for(j=0; j<Num_state; j++)
         delete branch_metric[i][j];
   for(i=0; i<Packet_length; i++)
      delete branch_metric[i];
   delete branch_metric;
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

