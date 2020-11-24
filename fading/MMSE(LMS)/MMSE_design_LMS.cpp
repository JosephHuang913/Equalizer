/***********************************************************/
/* Author: Chao-wang Huang                                 */
/* Date: Saturday, April 02, 2005                          */
/* An adaptive MMSE Equalizer is simulated                 */
/* 3-path Rayleigh fading channel with equal strength      */
/* Multipath weighting factor: {0.577, 0.577, 0.577}       */
/* Algorithm: Least Mean Square (LMS) Algorithm            */
/***********************************************************/

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

const float Pi = 3.14159265358979;
#define m	2											// channel memory
const int N = 1000;
#define num_packet 100								// number of packets simulated
const int Num_path = 3;
const int Num_state = pow(2,m);
const float mju = 0.0075;								// Step size parameter
const int Num_tap = 11;								// Tap number of adaptive MMSE EQ
//float CIR[3] = {0.577, 0.577, 0.577};			// Channel Weighting Factor
//float CIR[3] = {0.407, 0.815, 0.407};
float CIR[3] = {(1+cos(-2*Pi/3.1))/2.0, (1+cos(0.0))/2.0, (1+cos(2*Pi/3.1))/2.0};
const float vc = 0.0; 								/* speed of vehicle in km/hr */
const double C = 3.0E8;  							/* light speed */
const double fc = 2.0e9; 							/* carrier frequency */
const double sym_rate = 1E6; 						/* Symbol rate in symbol/sec */
const double Doppler = (vc*1000.0/3600.0)/(C/fc);  // Maximum Doppler frequency (Hz)

int main(void)
{
	time_t  t, start, end;
	int i, j, p, *data_bit, *ch, from_state, tmp_state;
   double out_sym_I, out_sym_Q, **MMSE, **Input, Output[2], *MSE;
   double snr, Eb_No, noise_pwr, noise[2], **Yk, e[2];
   double t1, t2, t3, fade1[2], fade2[2], fade3[2];
   FILE *mse, *records;

   start = time(NULL);
   printf("MSE Performance of adaptive MMSE EQ in Multipath Fading Channel\n");
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
   fprintf(records, "MSE Performance of adaptive MMSE EQ in Multipath Fading Channel\n");
   fprintf(records, "3-path Rayleigh fading channel with equal strength\n");
   fprintf(records, "Multipath weighting factor: {0.577, 0.577, 0.577}\n");
   fprintf(records, "Speed of the vehicle = %f (km/h)\n", vc);
   fprintf(records, "Carrier Frequency = %e (Hz)\n", fc);
   fprintf(records, "Maximum Doppler Frequency = %f (Hz)\n", Doppler);
   fprintf(records, "Transmission bit Rate = %e (bps)\n", sym_rate*1);
   fprintf(records, "f_d * t = %f\n", Doppler / sym_rate);
   fprintf(records, "number of bits of simulation = %d\n\n", N*num_packet);
   fflush(records);

   data_bit = new int[N];
   ch = new int[m+1];
   MSE = new double[N];		// Mean Square Error
   Yk = new double*[2];
   for(i=0; i<2; i++)
   	Yk[i] = new double[N];
   MMSE = new double*[2];	// Coefficients of adaptive MMSE EQ
   for(i=0; i<2; i++)
   	MMSE[i] = new double[Num_tap];
   Input = new double*[2];	// Input signal of adaptive MMSE EQ
   for(i=0; i<2; i++)
   	Input[i] = new double[Num_tap];

	srand((unsigned) time(&t));

/************************/
/* main simulation loop */
/************************/
   mse=fopen("MSE.log","w");
   for(snr=30; snr<=30; snr+=10)
   {
   	// noise power calculation
      Eb_No = (double)snr;
      noise_pwr = 1.0/pow(10.0, Eb_No/10.0);	// BPSK, Nyquist filter assumption
      printf("Eb_No = %f, ", Eb_No);
      t1 = 10.0;  	// Initial time of Rayleigh fading pattern
      t2 = 20.0;  	// Initial time of the 2nd path
      t3 = 30.0;		// Innitial time of the 3rd path
      //from_state = random(4);						// Initial State of Channel
      from_state = 0;

      for(i=0; i<N; i++)
      	MSE[i] = 0.0;

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
            from_state = tmp_state & (Num_state-1); 			// to_state of trellis diagram

            /* AWGN channel */
            AWGN_noise(0, noise_pwr, &noise[0]);
            //noise[0] = noise[1] = 0.0;
         	Yk[0][i] = out_sym_I + noise[0];
            Yk[1][i] = out_sym_Q + noise[1];
         }
/******************************************************************************/

/*****************************************************/
/* Adaptive MMSE EQ (LMS Algorithm)                  */
/*****************************************************/
			// Initialize the adaptive MMSE filter
         //for(i=0; i<N; i++)
           //	MSE[i] = 0.0;

         for(j=0; j<Num_tap; j++)
         {
         	Input[0][j] = Input[1][j] = 0.0;
            MMSE[0][j] = MMSE[1][j] = 0.0;
         }
         //MMSE[0][5] = 1.0;

     /*    cout << noise_pwr << endl;
         cout << endl;
         for(j=0; j<Num_tap; j++) 		// row index
         {
         	for(l=0; l<Num_tap; l++)   // column index
            	cout << P[0][j][l] << "+j" << P[1][j][l] << " ";
            cout << endl;
         }
         getch();
       */
			for(i=0; i<N; i++)
         {
         	// Input signal of adaptive MMSE EQ
         	for(j=Num_tap-1; j>0; j--)
            {
         		Input[0][j] = Input[0][j-1];
            	Input[1][j] = Input[1][j-1];
            }
            Input[0][0] = Yk[0][i];
            Input[1][0] = Yk[1][i];
            //Input[0][0] = 2*data_bit[i]-1;
            //Input[1][0] = 0.0;

      /*      cout << Yk[0][i] << "+j" << Yk[1][i] << endl;
            getch();

            for(j=0; j<Num_tap; j++)
            {
         		cout << Input[0][j] << "+j";
            	cout << Input[1][j] << endl;
            }
            getch();
        */
				// Compute output
            Output[0] = Output[1] = 0.0;
   	      for(j=0; j<Num_tap; j++)
      	   {
         		Output[0] += (MMSE[0][j]*Input[0][j] + MMSE[1][j]*Input[1][j]);
            	Output[1] += (MMSE[0][j]*Input[1][j] - MMSE[1][j]*Input[0][j]);
         	}

       /*     cout << Output[0] << "+j" << Output[1] << endl;
            getch();
        */
            // Compute error (BPSK)
            if(i<10)
              	//e[0] = e[1] = 0.0;
               ;
            else
            {

            e[0] = (2*data_bit[i-5]-1) - Output[0];
            e[1] = 0.0 - Output[1];

            //e[0] = Yk[0][i] - Output[0];
            //e[1] = Yk[1][i] - Output[1];

        /*    cout << (2*data_bit[i]-1) << endl;
            cout << e[0] << "+j" << e[1] << endl;
            getch();
          */
            // Mean Square Error
            MSE[i] += (pow(e[0],2) + pow(e[1],2));
            //MSE[i] += pow(e[0],2);
            //fprintf(mse, "%d %f\n", i, MSE[i]);
        /*
            cout << MSE[i] << endl;
            getch();
          */

            // Update Coefficients
            for(j=0; j<Num_tap; j++)
            {
            	MMSE[0][j] += mju*(Input[0][j]*e[0] + Input[1][j]*e[1]);
            	MMSE[1][j] += mju*(Input[1][j]*e[0] - Input[0][j]*e[1]);
            }

      /*      for(j=0; j<Num_tap; j++)
            	cout << MMSE[0][j] << "+j" << MMSE[1][j] << endl;
            cout << endl;
            getch();
        */




        		}
         }
/******************************************************************************/
    	}

      for(i=0; i<N; i++)
	      fprintf(mse, "%d %f\n", i, MSE[i]/(float)num_packet);
   }

   fclose(mse);
   delete data_bit;
   delete ch;
   delete MSE;
   for(i=0; i<2; i++)
   	delete Yk[i];
   delete Yk;
   for(i=0; i<2; i++)
   	delete MMSE[i];
   delete MMSE;
   for(i=0; i<2; i++)
   	delete Input[i];
   delete Input;
 
   end = time(NULL);
   printf("Total elapsed time: %.0f(sec)\n", difftime(end,start));
   fprintf(records, "\nTotal elapsed time: %.0f(sec)\n\n", difftime(end,start));
   fclose(records);
   printf("This program is ended. Press any key to continue.\n");
   getch();

   return 0;
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

