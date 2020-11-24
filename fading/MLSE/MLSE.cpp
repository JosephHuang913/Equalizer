/******************************************************/
/* Author: Chao-wang Huang                            */
/* Date: Friday, March 04, 2005                       */
/* An MLSE Equalizer is simulated                     */
/* 3-path Rayleigh fading channel with equal strength */
/* Multipath weighting factor: {0.577, 0.577, 0.577}  */
/* Algorithm: Soft Output Viterbi Algorithm (SOVA)    */
/******************************************************/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <conio.h>
#include <iostream.h>
#include <limits.h>
#include <float.h>

void JakesFading(double, float, double, int, double *);
void AWGN_noise(float, double, double *);
int Error_count(int, int);
void MLSE_SOVA(double **, int *, double *, int, int, int);
void CH_Trellis_diagram(int, int, float *, double *, double *, double *);

//#define Pi 	3.14159265358979
//#define n 	1											// regard channel as rate-1 convolutional code
//#define k	1
#define m	2											// channel memory
//#define N	1000	 									// packet length
const int N = 1000;
//#define K	512
#define num_packet 1000								// number of packets simulated
const int num_state = 4;							// number of channel states
const int num_in_sym = 2;							// number of input symbols (BPSK)
const int Num_path = 3;
float CIR[3] = {0.577, 0.577, 0.577};			// Channel Weighting Factor
//float CIR[3] = {0.407, 0.815, 0.407};
const float vc = 120.0; 								/* speed of vehicle in km/hr */
const double C = 3.0E8;  							/* light speed */
const double fc = 2.0e9; 							/* carrier frequency */
const double sym_rate = 1E6; 						/* Symbol rate in symbol/sec */
const double Doppler = (vc*1000.0/3600.0)/(C/fc);  // Maximum Doppler frequency (Hz)

struct ch_trel_diag			// Data structure of trellis diagram for each branch
		{
      	int from;			// from_state of trellis diagram
         int to;				// to_state of trellis diagram
         int in;				// input data symbol of trellis diagram
         float out[2];		// output symbol of trellis diagram
      };
struct ch_trel_diag ch_Trellis[1000][4][2];	// Trellis[time index][num_state][num_in_sym]

int main(void)
{
	time_t  t, start, end;
	int i, p, *data_bit, *Hk;
   int from_state, tmp_state, err_count[2], *Ak, *ch;
   double out_sym_I, out_sym_Q;
   double *bpsk, snr, Eb_No, noise_pwr, noise[2], **Yk, err_rate[2];
   double *LLR, t1, t2, t3, fade1[2], fade2[2], fade3[2];
   FILE *ber, *records;

   start = time(NULL);
   printf("BER Performance of MLSE in Multipath Fading Channel\n");
	cout << "3-path Rayleigh fading channel with equal strength" << endl;
	cout << "Multipath weighting factor: {0.577, 0.577, 0.577}" << endl;
   printf("Speed of the vehicle = %f (km/h)\n", vc);
   printf("Carrier Frequency = %e (Hz)\n", fc);
   printf("Maximum Doppler Frequency = %f (Hz)\n", Doppler);
   printf("Transmission bit Rate = %e (bps)\n", sym_rate*1);
   printf("f_d * t = %f\n", Doppler / sym_rate);
   printf("\nThis program is running. Don't close, please!\n\n");

   records = fopen("Records_MLSE.log", "a");
   fprintf(records, "BER Performance of MLSE in Multipath Fading Channel\n");
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

   data_bit = new int[N];
   ch = new int[m+1];
   Ak = new int[N];
   Hk = new int[N];
   bpsk = new double[N];
   Yk = new double*[2];
   for(i=0; i<2; i++)
   	Yk[i] = new double[N];
   LLR = new double[N];		// Log-likelihood Ratio

	srand((unsigned) time(&t));

/************************/
/* main simulation loop */
/************************/
   ber=fopen("ber_MLSE.log","w");
   for(snr=0; snr<=10; snr++)
   {
   	err_count[0] = err_count[1] = 0;
   	// noise power calculation
      Eb_No = (double)snr;
      noise_pwr = 1.0/pow(10.0, Eb_No/10.0);	// BPSK, Nyquist filter assumption
      printf("Eb_No = %f, ", Eb_No);
      t1 = 100.0;  	// Initial time of Rayleigh fading pattern
      t2 = 200.0;  	// Initial time of the 2nd path
      t3 = 300.0;		// Innitial time of the 3rd path
      from_state = random(4);						// Initial State of Channel
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

            JakesFading(fc, vc*1000/3600.0, t1, 2, &fade1[0]);
      		t1 += 1.0/sym_rate;
            JakesFading(fc, vc*1000/3600.0, t2, 2, &fade2[0]);
      		t2 += 1.0/sym_rate;
            JakesFading(fc, vc*1000/3600.0, t3, 2, &fade3[0]);
      		t3 += 1.0/sym_rate;

            // Calculate channel output symbol
         	out_sym_I = ch[0]*CIR[0]*fade1[0] + ch[1]*CIR[1]*fade2[0] + ch[2]*CIR[2]*fade3[0];
            out_sym_Q = ch[0]*CIR[0]*fade1[1] + ch[1]*CIR[1]*fade2[1] + ch[2]*CIR[2]*fade3[1];
            from_state = tmp_state & (num_state-1); 			// to_state of trellis diagram

            CH_Trellis_diagram(i, Num_path, CIR, fade1, fade2, fade3);

            /* AWGN channel */
            AWGN_noise(0, noise_pwr, &noise[0]);
         	Yk[0][i] = out_sym_I + noise[0];
            Yk[1][i] = out_sym_Q + noise[1];
         }

         MLSE_SOVA(Yk, Hk, LLR, N, num_state, num_in_sym);

         // Hard Decision
         for(i=0; i<N; i++)
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
      err_rate[0] = err_count[0] / (double)(N*num_packet);
      err_rate[1] = err_count[1] / (double)(N*num_packet);
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
   for(i=0; i<2; i++)
   	delete Yk[i];
   delete Yk;
   delete LLR;

   end = time(NULL);
   printf("Total elapsed time: %.0f(sec)\n", difftime(end,start));
   fprintf(records, "\nTotal elapsed time: %.0f(sec)\n\n", difftime(end,start));
   fclose(records);
   printf("This program is ended. Press any key to continue.\n");
   getch();

   return 0;
}

void MLSE_SOVA(double **Data_in, int *Data_out, double *Soft_out, int Packet_length, int Num_state, int Num_in)
{
/***************************************************************************/
/* MLSE: Soft Output Viterbi Algorithm (SOVA) for multipath fading channel */
/***************************************************************************/
	int i, j, l, q, survivor_state, pre_state;
   double **mju_f, *mju, *survival_metric, metric, ***branch_metric, **mju_b;
   double mju_tmp, survivor_metric;

   struct surv        			// Data structure of survival path
   		{
         	double metric;		// Path metric
            int data_in[N];	// input bit stream, K: packet length, global
            int state[N];		// state transition sequence, K: packet length, global
         };
   struct surv survival_path[num_state], survival_temp[num_state];
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
            branch_metric[i][j][l] = 0.0;
            branch_metric[i][j][l] += pow(Data_in[0][i]-ch_Trellis[i][j][l].out[0],2);
            branch_metric[i][j][l] += pow(Data_in[1][i]-ch_Trellis[i][j][l].out[1],2);
            metric = survival_path[j].metric + branch_metric[i][j][l];

            // find the survival path metric
            if(metric < survival_metric[ch_Trellis[i][j][l].to])
            {
            	survival_metric[ch_Trellis[i][j][l].to] = metric;

               // Record and refresh the survival path
               /*for(q=0; q<i; q++)
               {
               	survival_temp[ch_Trellis[i][j][l].to].data_in[q] = survival_path[j].data_in[q];
                  survival_temp[ch_Trellis[i][j][l].to].state[q] = survival_path[j].state[q];
               } */
               survival_temp[ch_Trellis[i][j][l].to].data_in[i] = l;
               survival_temp[ch_Trellis[i][j][l].to].state[i] = ch_Trellis[i][j][l].from;
            }
         }

      // Record and refresh the survival path
      for(j=0; j<Num_state; j++)					// to_state index
      {
      	survival_path[j].metric = survival_metric[j];
         mju_f[i+1][j] = survival_metric[j];
         /*for(q=0; q<=i; q++)
         {
         	survival_path[j].data_in[q] = survival_temp[j].data_in[q];
            survival_path[j].state[q] = survival_temp[j].state[q];
         } */
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

   for(j=0; j<Num_state; j++)		// to_state index
      survival_path[j].data_in[Packet_length-1] = survival_temp[j].data_in[(Packet_length-1)];

   for(j=0; j<Num_state; j++)		// to_state index
   {
   	pre_state = survival_temp[j].state[Packet_length-1];  // from state

   	for( q=Packet_length-2; q>=0; q--)
   	{
      	survival_path[j].data_in[q] = survival_temp[pre_state].data_in[q];
	      pre_state = survival_temp[pre_state].state[q];  // from state
   	}
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
         	metric = mju_b[i+1][ch_Trellis[i][j][l].to] + branch_metric[i][j][l];

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
      				 + mju_b[i+1][ch_Trellis[i][j][(survival_path[survivor_state].data_in[i]+1)%2].to];

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

void CH_Trellis_diagram(int time, int Num_path, float *CIR, double *fade1, double *fade2, double *fade3)
{
/*********************************************************/
/* Generate TrellisDiagram for time varying ISI channel */
/*********************************************************/
	int input, from_state, to_state, tmp_state, ch[3];
   int num_state = pow(2,Num_path-1), num_in_sym = 2;
   double out_sym_I, out_sym_Q;

   for(from_state=0; from_state<num_state; from_state++) // from_state of trellis diagram
   	for(input=0; input<num_in_sym; input++)				// input of trellis diagram
      {
      	tmp_state = from_state;
         tmp_state = (tmp_state << 1) ^ (input & 0x01);  // read input bit
         ch[0] = 2*(tmp_state & 0x01) - 1;
         ch[1] = 2*((tmp_state >> 1) & 0x01) - 1;
         ch[2] = 2*((tmp_state >> 2) & 0x01) - 1;
         out_sym_I = ch[0]*CIR[0]*fade1[0] + ch[1]*CIR[1]*fade2[0] + ch[2]*CIR[2]*fade3[0];
         out_sym_Q = ch[0]*CIR[0]*fade1[1] + ch[1]*CIR[1]*fade2[1] + ch[2]*CIR[2]*fade3[1];
         to_state = tmp_state & (num_state-1); 					// to_state of trellis diagram
         ch_Trellis[time][from_state][input].from = from_state;
         ch_Trellis[time][from_state][input].to = to_state;
         ch_Trellis[time][from_state][input].in = (2*input-1);
         ch_Trellis[time][from_state][input].out[0] = out_sym_I;
         ch_Trellis[time][from_state][input].out[1] = out_sym_Q;
      }
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

