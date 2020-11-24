/************************************************************/
/* Author: Chao-wang Huang                                  */
/* Date: Saturday, February 26, 2005                        */
/* An MAP Equalizer is simulated in Rayleigh fading channel */
/* 3-path Rayleigh fading channel with equal strength       */
/* Multipath weighting factor: {0.577, 0.577, 0.577}        */
/************************************************************/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <conio.h>
#include <iostream.h>

void JakesFading(double, float, double, int, double *);
void AWGN_noise(float, double, double *);
int Error_count(int, int);
void MAP_EQ_BCJR(double **, long double *, double, long double *, int, int, int);
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
const float vc = 120.0; 							/* speed of vehicle in km/hr */
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
	int i, p, *data_bit, *ch;
   int from_state, tmp_state, *Ak, err_count;
   double snr, Eb_No, noise_pwr, noise[2], **Yk, err_rate, out_sym_I, out_sym_Q;
   double t1, t2, t3, fade1[2], fade2[2], fade3[2];
   long double *LLR, *intrinsic;
   FILE *ber, *records;

   start = time(NULL);
   printf("BER Performance of MAP Equalizer in multipath fading channel\n");
   cout << "3-path Rayleigh fading channel with equal strength" << endl;
	cout << "Multipath weighting factor: {0.577, 0.577, 0.577}" << endl;
   printf("Speed of the vehicle = %f (km/h)\n", vc);
   printf("Carrier Frequency = %e (Hz)\n", fc);
   printf("Maximum Doppler Frequency = %f (Hz)\n", Doppler);
   printf("Transmission bit Rate = %e (bps)\n", sym_rate*1);
   printf("f_d * t = %f\n", Doppler / sym_rate);
   printf("number of bits of simulation = %d\n\n", N*num_packet);
   printf("This program is running. Don't close, please!\n\n");

   records = fopen("Records_MAPEQ.log", "a");
   fprintf(records, "BER Performance of MAP Equalizer in multipath fading channel\n");
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
   Yk = new double*[2];
   for(i=0; i<2; i++)
   	Yk[i] = new double[N];
   LLR = new long double[N];
   intrinsic = new long double[N];

	srand((unsigned) time(&t));

   for(i=0; i<N; i++)
   	intrinsic[i] = 0.0;

/************************/
/* main simulation loop */
/************************/
	ber=fopen("ber_MAPEQ.log","w");
   for(snr=0; snr<=10; snr++)
   {
   	err_count = 0;
   	// noise power calculation
      Eb_No = (double)snr;
      noise_pwr = 1.0/pow(10.0, Eb_No/10.0);	// BPSK, Nyquist filter assumption
      printf("Eb_No = %f, ", Eb_No);
      t1 = 100.0;  	// Initial time of Rayleigh fading pattern
      t2 = 200.0;  	// Initial time of the 2nd path
      t3 = 300.0;		// Innitial time of the 3rd path
      from_state = random(4);
      for(p=0; p<num_packet; p++)
      {
/*********************************************************************************/
/* BPSK mapping and equal strength multipath fading channel: {0.577,0.577,0.577} */
/*********************************************************************************/
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
/******************************************************************************/

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
      fprintf(ber, "%f %e\n", Eb_No, err_rate);
      fprintf(records, "%f %e\n", Eb_No, err_rate);
      fflush(records);
      fflush(ber);
   }

   delete data_bit;
   delete ch;
   delete Ak;
   for(i=0; i<2; i++)
   	delete Yk[i];
   delete Yk;
   delete LLR;
   delete intrinsic;

   end = time(NULL);
   printf("\nTotal elapsed time: %.0f(sec)\n", difftime(end,start));
   fprintf(records, "\nTotal elapsed time: %.0f(sec)\n\n", difftime(end,start));
   fclose(records);
   fclose(ber);
   printf("This program is ended. Press any key to continue.\n");
   getch();

   return 0;
}

void MAP_EQ_BCJR(double **Data_in, long double *intrinsic, double noise_pwr,
					  long double *LLR, int Packet_length, int Num_state, int Num_in)
{
/************************************************************************/
/* MAP Equalizer (modified BCJR Algorithm) for multipath fading channel */
/************************************************************************/
   const float pi = 3.14159265358979;
	int i, j, l;
   long double p_I, p_Q, *normal, **alpha, **beta, ***gamma, *delta, min, **a_priori;

   a_priori = new long double*[Packet_length];			// a-priori probability, a-priori[time index][input]
   for(i=0; i<Packet_length; i++)
   	a_priori[i] = new long double[Num_in];
   alpha = new long double*[Packet_length+1];			// alpha[time index][state]
   for(i=0; i<=Packet_length; i++)
      alpha[i] = new long double[Num_state];
   beta = new long double*[Packet_length+1];       	// beta[time index][state]
   for(i=0; i<=Packet_length; i++)
      beta[i] = new long double[Num_state];
   gamma = new long double**[Packet_length];			// gamma[time index][state][input]
   for(i=0; i<Packet_length; i++)
      gamma[i] = new long double*[Num_state];
   for(i=0; i<Packet_length; i++)
   	for(j=0; j<Num_state; j++)
         gamma[i][j] = new long double[Num_in];
   normal = new long double[Packet_length+1];			// renormalization of BCJR
   delta = new long double[Num_in];

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
         	p_I = exp(-pow(Data_in[0][i]-ch_Trellis[i][j][l].out[0],2)/noise_pwr)/sqrt(pi*noise_pwr);
            p_Q = exp(-pow(Data_in[1][i]-ch_Trellis[i][j][l].out[1],2)/noise_pwr)/sqrt(pi*noise_pwr);
            gamma[i][j][l] = a_priori[i][l] * p_I * p_Q;		// gamma[time index][state][input]
         }

   // calculate alpha[time index][state]
   for(i=1; i<=Packet_length; i++)					// time index
   {
   	for(j=0; j<Num_state; j++)						// from_state index
      	for(l=0; l<Num_in; l++)						// input bit
         	alpha[i][ch_Trellis[i-1][j][l].to] += alpha[i-1][j] * gamma[i-1][j][l];

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
         	beta[i][j] += beta[i+1][ch_Trellis[i][j][l].to] * gamma[i][j][l];

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
         	delta[0] = alpha[i][j] * gamma[i][j][l] * beta[i+1][ch_Trellis[i][j][l].to];

            if((delta[0] < min && delta[0] != 0.0) || min == 0.0)
            	min = delta[0];
         }

      if(min == 0.0 || min > 1.0)	// if all else fails, make min real small
      	min = 1E-100;

      delta[0] = delta[1] = 0.0;
      for(j=0; j<Num_state; j++)		// from_state index
      	for(l=0; l<Num_in; l++)	// input bit
         	delta[l] += alpha[i][j] * gamma[i][j][l] * beta[i+1][ch_Trellis[i][j][l].to];

      if(delta[1] == 0.0)
      	delta[1] = min;
      if(delta[0] == 0.0)
      	delta[0] = min;

      LLR[i] = logl(delta[1]/delta[0]);
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

