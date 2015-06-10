// DRB 02-Dec-2009
// Using the TI FFT function

#include <stdio.h>
#include <math.h>


#define N 4096
#define RADIX 2
#define PI 3.14159265358979

#define NS_TICK 125000 //ns in a clock tick
// complex typedef
typedef struct {
	float re,im;
} COMPLEX;

// align data (nothing works if you omit these pragma!!!!!!!!!!)
#pragma DATA_ALIGN(w,sizeof(COMPLEX))   //align w
#pragma DATA_ALIGN(x,sizeof(COMPLEX))	//align x
#pragma DATA_ALIGN(mult, sizeof(COMPLEX)) //align mult

#pragma DATA_SECTION(cosbuf,".mydata")

// function prototypes
void cfftr2_dit(COMPLEX*, COMPLEX*, short);
void bitrev(COMPLEX*, short*, int);
void digitrev_index(short*, int, int);

// global variables
COMPLEX w[N/RADIX];			    // array of complex twiddle factors
COMPLEX x[N]; 				    // array of complex FFT input/output samples

COMPLEX mult[N];  				//Array of samples that are multiplied


float DELTA = 2.0*PI/N;
short iw[N/2],ix[N];

// indices for bit reversal
int i,j = 0 ,n, imax = 0,maximum = 0, imax_shift = 0;

//I am setting this as a pointer because for some reason it isn't modifying
int maximum_shift = 0;

int add = 1;

//Keep in mind that cosbuf is i nactuality NS_TICK long.
far float cosbuf[3 * NS_TICK];
float *sinbuf;
float temp;

//debug
float a=0, c= 0, b=0, d = 0;

void main(void)
{
	sinbuf = &(cosbuf[NS_TICK]);	//initialize the sin buffer location

	int max_shift = 0;

	// compute first N/2 twiddle factors
 	for(i=0;i<N/RADIX;i++){
   		w[i].re = cos(DELTA*i);		    
   		w[i].im = sin(DELTA*i);		// negative imag component
  	}	

 	//Initialize cos and sin buffers
 	for(i = 0; i < NS_TICK; i++){
 		cosbuf[i] =cos(2 * PI * i / (NS_TICK));
 		sinbuf[i] =sin(2 * PI * i / (NS_TICK));
 	}

	// initialize complex FFT input array with modulated sinc pulse
	for (i=0;i<N;i++){
		//x[i].re = cos(PI/2*i);
		if(i != 0){
			x[i].re = cos(PI / 2 * i) *  sin(i) / i;
		}
		else
			x[i].re = 1;

		x[i].im = 0;
		mult[i].re = 0;
		mult[i].im = 0;

		//Debug if shifting proper.
		if(x[i].re > maximum){
			maximum = x[i].re;
			imax = i;
		}
	}

	//Inmaximumthe FFT
 	digitrev_index(iw,N/RADIX,RADIX);//produces index for bitrev() W
 	bitrev(w,iw,N/RADIX);		    //bit reverse W

 	//Run FFT
 	cfftr2_dit(x,w,N);			    //TI floating-pt complex FFT
 	digitrev_index(ix, N, RADIX);	    //produces index for bitrev() X 
 	bitrev(x,ix,N);			    //freq scrambled->bit-reverse X



 	for(i = 0; i < N; i++){

 		mult[i].re = x[i].re * cosbuf[j] - x[i].im * sinbuf[j];
 		mult[i].im = -(x[i].im * cosbuf[j] + x[i].re * sinbuf[j]);

 		j += add;
 		while(j >= NS_TICK)
 			j -= NS_TICK;
 	}

 	/**
 	*		The routine can be used to implement Inverse-FFT by any ONE of
 	*		the following methods:
 	*
 	*		1.Inputs (x) are replaced by their Complex-conjugate values
 	*		  Output values are divided by N
 	*		2.FFT Coefficients (w) are replaced by their Complex-conjugates
 	*		  Output values are divided by N
 	*		3.Swap Real and Imaginary values of input
 	*		4.Swap Real and Imaginary values of output
 	*/


 	cfftr2_dit(mult,w,N);			    //TI floating-pt complex FFT
 	digitrev_index(ix, N, RADIX);	    //produces index for bitrev() X
 	bitrev(mult,ix,N);			    //freq scrambled->bit-reverse X

 	for(i = 0; i < N; i++){
 		mult[i].re /= N;
 		mult[i].im /= N;
 	}

/*
 	for(i = 0; i < N; i++){
		if(mult[i].re > maximum_shift){
			maximum_shift = mult[i].re;
			imax_shift = i;
 		}
		if(mult[i].re > max_shift)
			max_shift = mult[i].re;
 	}
*/
 	i = 0;
 	while(i < N)
 	{
		if(mult[i].re > maximum_shift){
			maximum_shift = mult[i].re;
			imax_shift = i;
 		}
		if(mult[i].re > max_shift)
			max_shift = mult[i].re;

		i++;
 	}

	while(1){
		j += add;
 	}

}
