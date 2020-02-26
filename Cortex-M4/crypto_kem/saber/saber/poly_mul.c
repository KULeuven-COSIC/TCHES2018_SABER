#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include "poly_mul.h"


// define ASM for assembly, KARA for karatsuba
//#define ASM
#define KARA
// extra memory for unrolled karatsuba
uint16_t kara_tmp[16];

void karatsuba_simple(const uint16_t* a_1,const uint16_t* b_1, uint16_t* result_final);

void pol_mul(uint16_t* a, uint16_t* b, uint16_t* res, uint16_t p, uint32_t n)
{ 
	uint32_t i;

	uint16_t c[2*SABER_N];
	for (i = 0; i < 2*n; i++) c[i] = 0;

//-------------------normal multiplication-----------------

#if !defined(ASM) && !defined(KARA)

	uint32_t j;
	for (i = 0; i < n; i++) {
		for (j = 0; j < n; j++) {
			c[i+j]=( c[i+j] + (a[i] * b[j]) )&(p-1);
		}
	}

#elif defined(ASM) && defined(KARA)

	unrolled_kara_mem(a, b, c, SABER_N/2);

#else

	school_book_mul2(a, b, c, SABER_N);

#endif

	//---------------reduction-------
	for(i=n;i<2*n-1;i++){
		c[i-n]=(c[i-n]-c[i])&(p-1);
	}
	c[n-1] &= (p-1);
//----------------------copy result back----------------

	for(i=0; i<n; i++){
            res[i] = c[i];
	}

}

void karatsuba_simple(const uint16_t* a_1,const uint16_t* b_1, uint16_t* result_final){//uses 10 registers

	uint16_t N=64;
	uint16_t d01[N/2-1];
	uint16_t d0123[N/2-1];
	uint16_t d23[N/2-1];
	uint16_t result_d01[N-1];	

	int32_t i,j;

	memset(result_d01,0,(N-1)*sizeof(uint16_t));
	memset(d01,0,(N/2-1)*sizeof(uint16_t));
	memset(d0123,0,(N/2-1)*sizeof(uint16_t));
	memset(d23,0,(N/2-1)*sizeof(uint16_t));
	memset(result_final,0,(2*N-1)*sizeof(uint16_t));

	uint16_t acc1,acc2,acc3,acc4,acc5,acc6,acc7,acc8,acc9,acc10;


	for (i = 0; i < N/4; i++) {
		acc1=a_1[i];//a0
		acc2=a_1[i+N/4];//a1
		acc3=a_1[i+2*N/4];//a2
		acc4=a_1[i+3*N/4];//a3	
		for (j = 0; j < N/4; j++) {

			acc5=b_1[j];//b0
			acc6=b_1[j+N/4];//b1

			result_final[i+j+0*N/4]=result_final[i+j+0*N/4]+acc1*acc5;
			result_final[i+j+2*N/4]=result_final[i+j+2*N/4]+acc2*acc6;

			acc7=acc5+acc6;//b01
			acc8=acc1+acc2;//a01
			d01[i+j]=d01[i+j] + acc7*acc8;
	//--------------------------------------------------------

			acc7=b_1[j+2*N/4];//b2
			acc8=b_1[j+3*N/4];//b3			
			result_final[i+j+4*N/4]=result_final[i+j+4*N/4]+acc7*acc3;

			result_final[i+j+6*N/4]=result_final[i+j+6*N/4]+acc8*acc4;

			acc9=acc3+acc4;
			acc10=acc7+acc8;
			d23[i+j]=d23[i+j] + acc9*acc10;
	//--------------------------------------------------------

			acc5=acc5+acc7;//b02
			acc7=acc1+acc3;//a02
			result_d01[i+j+0*N/4]=result_d01[i+j+0*N/4]+acc5*acc7;

			acc6=acc6+acc8;//b13
			acc8=acc2+acc4;			
			result_d01[i+j+ 2*N/4]=result_d01[i+j+ 2*N/4]+acc6*acc8;

			acc5=acc5+acc6;
			acc7=acc7+acc8;
			d0123[i+j]=d0123[i+j] + acc5*acc7;
		}
	}

//------------------2nd last stage-------------------------

	for(i=0;i<N/2-1;i++){
		d0123[i]=d0123[i]-result_d01[i+0*N/4]-result_d01[i+2*N/4];
		d01[i]=d01[i]-result_final[i+0*N/4]-result_final[i+2*N/4];
		d23[i]=d23[i]-result_final[i+4*N/4]-result_final[i+6*N/4];
	}

	for(i=0;i<N/2-1;i++){
		result_d01[i+1*N/4]=result_d01[i+1*N/4]+d0123[i];
		result_final[i+1*N/4]=result_final[i+1*N/4]+d01[i];
		result_final[i+5*N/4]=result_final[i+5*N/4]+d23[i];
	}

//------------Last stage---------------------------
	for(i=0;i<N-1;i++){
		result_d01[i]=result_d01[i]-result_final[i]-result_final[i+N];
	}
	
	for(i=0;i<N-1;i++){
		result_final[i+1*N/2]=result_final[i+1*N/2]+result_d01[i];//-result_d0[i]-result_d1[i];		
	}

}


void karatsuba_simple1(const uint16_t* a_1,const uint16_t* b_1, uint16_t* result_final){ //uses 9 registers

	uint16_t N=64;
	// uint16_t d01[N/2-1];
	// uint16_t d0123[N/2-1];
	// uint16_t d23[N/2-1];
	// uint16_t result_d01[N-1];

	//d01=d[0] d0123=d[N/2-1] d23=d[N-2] result_d01=d[N+N/2-3]
	uint16_t d[(N/2-1)+(N/2-1)+(N/2-1)+(N-1)];//156 coeffs

	uint16_t * d01 = d;
	uint16_t * d0123 = &d[N/2-1];
	uint16_t * d23 = &d[N-2];
	uint16_t * result_d01 = &d[N+N/2-3];

	int32_t i,j;

	// memset(result_d01,0,(N-1)*sizeof(uint16_t));
	// memset(d01,0,(N/2-1)*sizeof(uint16_t));
	// memset(d0123,0,(N/2-1)*sizeof(uint16_t));
	// memset(d23,0,(N/2-1)*sizeof(uint16_t));
	memset(d,0,((N/2-1)+(N/2-1)+(N/2-1)+(N-1))*sizeof(uint16_t));

	memset(result_final,0,(2*N-1)*sizeof(uint16_t));

	uint16_t acc1,acc2,acc3,acc4,acc5,acc6,acc7,acc8,acc9;
	uint16_t acc11,acc21,acc31,acc41,acc51,acc61,acc71,acc81,acc91;

	// karatsuba_simple_asm( a_1, b_1, result_final, d );//not working

	for (i = 0; i < N/4; i = i+2) {
		acc1=a_1[i];//a0
		acc2=a_1[i+N/4];//a1
		acc3=a_1[i+2*N/4];//a2
		acc4=a_1[i+3*N/4];//a3
		acc11=a_1[i+1];//a0
		acc21=a_1[i+N/4+1];//a1
		acc31=a_1[i+2*N/4+1];//a2
		acc41=a_1[i+3*N/4+1];//a3
		for (j = 0; j < N/4; j = j+2) {
			acc5=b_1[j];//b0
			acc6=b_1[j+N/4];//b1
			acc51=b_1[j+1];//b0-1
			acc61=b_1[j+N/4+1];//b1-1

			result_final[i+j+0*N/4] += acc1*acc5;
			result_final[i+j+0*N/4+1] += acc1*acc51;
			result_final[i+j+0*N/4+1] += acc11*acc5;
			result_final[i+j+0*N/4+2] += acc11*acc51;

			result_final[i+j+2*N/4] += acc2*acc6;
			result_final[i+j+2*N/4+1] += acc2*acc61;
			result_final[i+j+2*N/4+1] += acc21*acc6;
			result_final[i+j+2*N/4+2] += acc21*acc61;

			acc7=b_1[j+2*N/4];//b2
			acc8=b_1[j+3*N/4];//b3
			acc71=b_1[j+2*N/4+1];//b2-1
			acc81=b_1[j+3*N/4+1];//b3-1

			result_final[i+j+4*N/4] += acc7*acc3;
			result_final[i+j+4*N/4+1] += acc71*acc3;
			result_final[i+j+4*N/4+1] += acc7*acc31;
			result_final[i+j+4*N/4+2] += acc71*acc31;

			result_final[i+j+6*N/4] += acc8*acc4;
			result_final[i+j+6*N/4+1] += acc81*acc4;
			result_final[i+j+6*N/4+1] += acc8*acc41;
			result_final[i+j+6*N/4+2] += acc81*acc41;

			acc7=acc5+acc7;//b02
			acc9=acc1+acc3;//a02
			acc71=acc51+acc71;//b02-1
			acc91=acc11+acc31;//a02-1

			//d01=d[i+j] d0123=d[N/2-1] d23=d[N-2] result_d01=d[N+N/2-3]
			d[N+N/2-3+i+j+0*N/4] += acc9*acc7;
			d[N+N/2-3+i+j+0*N/4+1] += acc9*acc71;
			d[N+N/2-3+i+j+0*N/4+1] += acc91*acc7;
			d[N+N/2-3+i+j+0*N/4+2] += acc91*acc71;

			acc8=acc8+acc6;//b13
			acc9=acc2+acc4;//a13
			acc81=acc81+acc61;//b13-1
			acc91=acc21+acc41;//a13-1

			//d01=d[i+j] d0123=d[N/2-1] d23=d[N-2] result_d01=d[N+N/2-3]			
			d[N+N/2-3+i+j+ 2*N/4] += acc8*acc9;
			d[N+N/2-3+i+j+ 2*N/4+1] += acc81*acc9;
			d[N+N/2-3+i+j+ 2*N/4+1] += acc8*acc91;
			d[N+N/2-3+i+j+ 2*N/4+2] += acc81*acc91;

			acc8=acc8-acc6;
			acc7=acc7-acc5;
			acc81=acc81-acc61;
			acc71=acc71-acc51;
			acc6=acc5+acc6;//b01
			acc9=acc1+acc2;//a01
			acc61=acc51+acc61;//b01-1
			acc91=acc11+acc21;//a01-1
			//acc5 free

			//d01=d[i+j] d0123=d[N/2-1] d23=d[N-2] result_d01=d[N+N/2-3]
			d[i+j] += acc9*acc6;
			d[i+j+1] += acc9*acc61;
			d[i+j+1] += acc91*acc6;
			d[i+j+2] += acc91*acc61;
	//--------------------------------------------------------

			acc9=acc3+acc4;
			acc8=acc7+acc8;
			acc91=acc31+acc41;
			acc81=acc71+acc81;

			//d01=d[i+j] d0123=d[N/2-1] d23=d[N-2] result_d01=d[N+N/2-3]
			d[N-2+i+j] +=acc8*acc9;
			d[N-2+i+j+1] += acc81*acc9;
			d[N-2+i+j+1] += acc8*acc91;
			d[N-2+i+j+2] += acc81*acc91;
	//--------------------------------------------------------

			acc9=acc9+acc2+acc1;
			acc8=acc8+acc6;
			acc91=acc91+acc21+acc11;
			acc81=acc81+acc61;

			//d01=d[i+j] d0123=d[N/2-1] d23=d[N-2] result_d01=d[N+N/2-3]
			d[N/2-1+i+j] += acc8*acc9;
			d[N/2-1+i+j+1] += acc81*acc9;
			d[N/2-1+i+j+1] += acc8*acc91;
			d[N/2-1+i+j+2] += acc81*acc91;
		}
	}

//------------------2nd last stage-------------------------

	for(i=0;i<N/2-1;i++){
		d0123[i]=d0123[i]-result_d01[i+0*N/4]-result_d01[i+2*N/4];
		d01[i]=d01[i]-result_final[i+0*N/4]-result_final[i+2*N/4];
		d23[i]=d23[i]-result_final[i+4*N/4]-result_final[i+6*N/4];
	}

	for(i=0;i<N/2-1;i++){
		result_d01[i+1*N/4]=result_d01[i+1*N/4]+d0123[i];
		result_final[i+1*N/4]=result_final[i+1*N/4]+d01[i];
		result_final[i+5*N/4]=result_final[i+5*N/4]+d23[i];
	}

//------------Last stage---------------------------
	for(i=0;i<N-1;i++){
		result_d01[i]=result_d01[i]-result_final[i]-result_final[i+N];
	}
	
	for(i=0;i<N-1;i++){
		result_final[i+1*N/2]=result_final[i+1*N/2]+result_d01[i];//-result_d0[i]-result_d1[i];		
	}

}

void toom_cook_4way_mem2(const uint16_t* a1,const uint16_t* b1, uint16_t* result)
{//This one reuses memory fo a and b hence less memory requiement. Suitable for m0

	int16_t i;
	int16_t small_len=SABER_N/4;
	//uint16_t result_final[2*SABER_N];

//----------------array declaration to hold results--------------------

	//uint16_t w_m[896+128];
	//uint16_t w_m[1408];
	uint16_t w_m[896+128];

//----------------array declaration to hold results ends--------------------
	//--------------------these data are created for place holding---------

	//uint16_t a[small_len],b[small_len];

	int16_t acc1, acc2, acc3, acc4, acc5, acc6, acc7, acc8;

	//--------------------these data are created for place holding ends---------

	int16_t inv3=43691;
	int16_t inv9=36409;
	int16_t inv15=61167;

	int16_t int45=45;
	int16_t int30=30;
	int16_t int0=0;

	//for(i=0; i<2*SABER_N; i++){
	//	result_final[i] = 0;
	//}

	memset(w_m,0,896*sizeof(int16_t));
//----------------------------------------------------------
	unrolled_kara_mem(a1, b1, w_m+768, small_len/2);
//----------------------------------------------------------
	unrolled_kara_mem(a1+small_len*3, b1+small_len*3, w_m, small_len/2);
//----------------------------------------------------------

//----------------------------------------------------------

	for(i=0;i<small_len;i++){
		w_m[i+896]=8*a1[i]-4*a1[i+64]+2*a1[i+128]-a1[i+192];
		w_m[i+960]=8*b1[i]-4*b1[i+64]+2*b1[i+128]-b1[i+192];
	}

	unrolled_kara_mem(w_m+896, w_m+960, w_m+640, small_len/2);

//----------------------------------------------------------
	for(i=0;i<small_len;i++){
		w_m[i+896]=8*a1[i]+4*a1[i+64]+2*a1[i+128]+a1[i+192];
		w_m[i+960]=8*b1[i]+4*b1[i+64]+2*b1[i+128]+b1[i+192];
	}
	unrolled_kara_mem(w_m+896, w_m+960, w_m+512, small_len/2);
//----------------------------------------------------------
	for(i=0;i<small_len;i++){
		w_m[i+896]=a1[i]-a1[i+64]+a1[i+128]-a1[i+192];
		w_m[i+960]=b1[i]-b1[i+64]+b1[i+128]-b1[i+192];
	}
	unrolled_kara_mem(w_m+896, w_m+960, w_m+384, small_len/2);
//----------------------------------------------------------
	for(i=0;i<small_len;i++){
		w_m[i+896]=a1[i]+a1[i+64]+a1[i+128]+a1[i+192];
		w_m[i+960]=b1[i]+b1[i+64]+b1[i+128]+b1[i+192];
	}
	unrolled_kara_mem(w_m+896, w_m+960, w_m+256, small_len/2);
//----------------------------------------------------------
	for(i=0;i<small_len;i++){
		w_m[i+896]=a1[i]+2*a1[i+64]+4*a1[i+128]+8*a1[i+192];
		w_m[i+960]=b1[i]+2*b1[i+64]+4*b1[i+128]+8*b1[i+192];
	}
	unrolled_kara_mem(w_m+896, w_m+960, w_m+128, small_len/2);
	/*************************MULTIPLICATION ENDS***********************************/
			
	//memset(w_m+896,0,512*sizeof(int16_t));

	for(i=0;i<2*small_len-1;i++){

		acc1=w_m[i];
		acc2=w_m[i+128];
		acc3=w_m[i+256];
		acc4=w_m[i+384];
		acc5=w_m[i+512];
		acc6=w_m[i+640];
		acc7=w_m[i+768];

		acc2= acc2+acc5;//w2 <- w2+w5
		acc6= acc6-acc5;// w6 <- w6-w5
		acc4= acc4-acc3;// w4 <- w4-w3
		
		acc5= acc5-acc1;// w5 <- w5-w1
		//temp = acc7<<6; //temp <- 64*w7
		acc8 = acc7<<6; //temp <- 64*w7

		acc5= acc5-acc8;// w5 <- w5-64*w7
		acc4 = acc4>>1; //w4 <- w4/2
		acc3 = acc3+acc4;//w3 <- w3+w4

		acc8 = acc5<<1; //temp <- 2*w5
		acc5= acc6+acc8;//w5 <- 2*w5+w6

		acc8 = acc3<<6; //temp <- 64*w3
		acc8 = acc3+acc8; //temp <- 65*w3
		acc2= acc2-acc8;// w2 <- w2-65*w3

		acc3= acc3-acc7;// w3 <- w3-w7

		acc3= acc3-acc1;// w3 <- w3-w1

		acc8 = acc3*int45; //temp <- 45*w3
		acc2 = acc2+acc8; //w2 <- w2+45*w3


		acc8 = acc3<<3; //temp <- 8*w3
		acc5= acc5-acc8;//w5 <- w5-8*w3

		acc5 = acc5*inv3; //w5 <- w5*1/3
		acc5 = acc5>>3; //w5 <- w5*1/8 ---> w5=w5/24

		acc6 = acc2+acc6; //w6 <- w6+w2
		acc8 = acc4<<4; //temp <- 16*w4
		acc2 = acc2+acc8; //w2 <- w2+16*w4

		acc2 = acc2*inv9; //w2 <- w2*1/9
		acc2 = acc2>>1; //w2 <- w2*1/2 ---> w2=w2/18

		acc3= acc3-acc5;//w3 <- w3-w5
		
		acc4 = acc4+acc2; //w4 <- w4+w2

		acc4 = int0-acc4; //w4 <- -(w4+w2)

		acc8 = acc2*int30; //temp <- w2*30
		acc6= acc8-acc6;//w6 <- 30*w2-w6

		acc6 = acc6*inv15; //w6 <- w6*1/15
		acc6 = acc6>>2; //w6 <- w6*1/4 ---> w6=w6/60

		acc2= acc2-acc6;//w2 <- w2-w6

		/*result_final[0*small_len+i]= (result_final[0*small_len+i]+ acc7);
		result_final[1*small_len+i]= (result_final[1*small_len+i]+ acc6);
		result_final[2*small_len+i]= (result_final[2*small_len+i]+ acc5);
		result_final[3*small_len+i]= (result_final[3*small_len+i]+ acc4);
		result_final[4*small_len+i]= (result_final[4*small_len+i]+ acc3);
		result_final[5*small_len+i]= (result_final[5*small_len+i]+ acc2);
		result_final[6*small_len+i]= (result_final[6*small_len+i]+ acc1);
		*/
		/*
		w_m[0*small_len+i+896]= (w_m[0*small_len+i+896]+ acc7-acc3);
		w_m[1*small_len+i+896]= (w_m[1*small_len+i+896]+ acc6-acc2);
		w_m[2*small_len+i+896]= (w_m[2*small_len+i+896]+ acc5-acc1);
		if(i>=64)
			w_m[3*small_len+i+896-SABER_N]= (w_m[3*small_len+i+896-SABER_N]- acc4);
		else
			w_m[3*small_len+i+896]= (w_m[3*small_len+i+896]+ acc4);
		*/

		result[0*small_len+i]= (result[0*small_len+i]+ acc7-acc3)&(SABER_Q-1);
		result[1*small_len+i]= (result[1*small_len+i]+ acc6-acc2)&(SABER_Q-1);
		result[2*small_len+i]= (result[2*small_len+i]+ acc5-acc1)&(SABER_Q-1);
		if(i>=64)
			result[3*small_len+i-SABER_N]= (result[3*small_len+i-SABER_N]- acc4)&(SABER_Q-1);
		else
			result[3*small_len+i]= (result[3*small_len+i]+ acc4)&(SABER_Q-1);
		
	}
}


//void toom_cook_4way_mem(const uint16_t* a1,const uint16_t* b1, uint16_t* result)
//{//the output is in a consecutive memory. Does not reuse memory for speed. Suitable for m4
//	uint16_t w_m[896+640];

	// int16_t inv3=43691;
	// int16_t inv9=36409;
	// int16_t inv15=61167;

	// int16_t int45=45;
	// int16_t int30=30;
	// int16_t int0=0;

	//  // int16_t i;
	//  // int16_t acc1, acc2, acc3, acc4, acc5, acc6, acc7, acc8, acc9, acc10;

//	toom_cook_4way_mem_asm(a1, b1, result, w_m);

//loop1
	// for(i=0;i<(SABER_N/4);i++){
		
	// 	acc3=a1[(SABER_N/4)*0+i];
	// 	acc4=a1[(SABER_N/4)*1+i];
	// 	acc5=a1[(SABER_N/4)*2+i];
	// 	acc6=a1[(SABER_N/4)*3+i];

	// 	acc7=acc3+acc5;//a0+a2
	// 	acc8=acc4+acc6;//a1+a3

	// 	acc9=acc7+acc8;//a0+a1+a2+a3
	// 	acc10=acc7-acc8;

	// 	w_m[i+64+896]=acc9;
	// 	w_m[i+128+896]=acc10;

	// 	acc7=acc3*4+acc5;//4*a0+a2
	// 	acc7=acc7*2;//8*a0+2*a2
	// 	acc8=acc4*4+acc6;//4*a1+a3

	// 	acc9=acc7+acc8;
	// 	acc10=acc7-acc8;				

	// 	w_m[i+192+896]=acc9;
	// 	w_m[i+256+896]=acc10;

	// 	acc7=2*acc6+acc5;
	// 	acc7=acc7*2+acc4;	
	// 	acc7=acc7*2+acc3;

	// 	w_m[i+896]=acc7;
	// }
//loop2
	// for(i=0;i<(SABER_N/4);i++){
		
	// 	acc3=b1[(SABER_N/4)*0+i];
	// 	acc4=b1[(SABER_N/4)*1+i];
	// 	acc5=b1[(SABER_N/4)*2+i];
	// 	acc6=b1[(SABER_N/4)*3+i];

	// 	acc7=acc3+acc5;//a0+a2
	// 	acc8=acc4+acc6;//a1+a3

	// 	acc9=acc7+acc8;//a0+a1+a2+a3
	// 	acc10=acc7-acc8;
		
	// 	w_m[i+64+320+896]=acc9;
	// 	w_m[i+128+320+896]=acc10;

	// 	acc7=acc3*4+acc5;//4*a0+a2
	// 	acc7=acc7*2;//8*a0+2*a2
	// 	acc8=acc4*4+acc6;//4*a1+a3

	// 	acc9=acc7+acc8;
	// 	acc10=acc7-acc8;				

	// 	w_m[i+192+320+896]=acc9;
	// 	w_m[i+256+320+896]=acc10;

	// 	acc7=2*acc6+acc5;
	// 	acc7=acc7*2+acc4;	
	// 	acc7=acc7*2+acc3;
		
	// 	w_m[i+320+896]=acc7;
	// }

	//memset(w_m,0,(896)*sizeof(int16_t));
	//unrolled_kara_mem(a1, b1, w_m+768, (SABER_N/4)/2);

	// unrolled_kara_mem(w_m+256+896, w_m+256+320+896, w_m+640, (SABER_N/4)/2);

	// unrolled_kara_mem(w_m+192+896, w_m+192+320+896, w_m+512, (SABER_N/4)/2);

	// unrolled_kara_mem(w_m+128+896,w_m+128+320+896, w_m+384, (SABER_N/4)/2);

	// unrolled_kara_mem(w_m+64+896, w_m+64+320+896, w_m+256, (SABER_N/4)/2);

	// unrolled_kara_mem(w_m+896, w_m+320+896, w_m+128, (SABER_N/4)/2);

	// unrolled_kara_mem(a1+(SABER_N/4)*3, b1+(SABER_N/4)*3, w_m, (SABER_N/4)/2);
//  // char output[32];
//loop3
	// for(i=0;i<2*(SABER_N/4)-1;i++){

	// 	acc1=w_m[i];
	// 	acc2=w_m[i+128];
	// 	acc3=w_m[i+256];
	// 	acc4=w_m[i+384];
	// 	acc5=w_m[i+512];
	// 	acc6=w_m[i+640];
	// 	acc7=w_m[i+768];

	// 	// sprintf((char *)output, "%d, ", (int)w_m[i]);
	// 	// send_USART_str(output);
	// 	// sprintf((char *)output, "%d, ", (int)w_m[i+128]);
	// 	// send_USART_str(output);
	// 	// sprintf((char *)output, "%d, ", (int)w_m[i+768]);
	// 	// send_USART_str(output);

	// 	acc2= acc2+acc5;//w2 <- w2+w5
	// 	acc6= acc6-acc5;// w6 <- w6-w5
	// 	acc4= acc4-acc3;// w4 <- w4-w3
		
	// 	acc5= acc5-acc1;// w5 <- w5-w1
	// 	//temp = acc7<<6; //temp <- 64*w7
	// 	acc8 = acc7<<6; //temp <- 64*w7

	// 	acc5= acc5-acc8;// w5 <- w5-64*w7
	// 	acc4 = acc4>>1; //w4 <- w4/2
	// 	acc3 = acc3+acc4;//w3 <- w3+w4

	// 	acc8 = acc5<<1; //temp <- 2*w5
	// 	acc5= acc6+acc8;//w5 <- 2*w5+w6

	// 	acc8 = acc3<<6; //temp <- 64*w3
	// 	acc8 = acc3+acc8; //temp <- 65*w3
	// 	acc2= acc2-acc8;// w2 <- w2-65*w3

	// 	acc3= acc3-acc7;// w3 <- w3-w7

	// 	acc3= acc3-acc1;// w3 <- w3-w1

	// 	acc8 = acc3*45; //temp <- 45*w3
	// 	acc2 = acc2+acc8; //w2 <- w2+45*w3

	// 	acc8 = acc3<<3; //temp <- 8*w3
	// 	acc5= acc5-acc8;//w5 <- w5-8*w3

	// 	acc5 = acc5*43691; //w5 <- w5*1/3
	// 	acc5 = acc5>>3; //w5 <- w5*1/8 ---> w5=w5/24

	// 	acc6 = acc2+acc6; //w6 <- w6+w2
	// 	acc8 = acc4<<4; //temp <- 16*w4
	// 	acc2 = acc2+acc8; //w2 <- w2+16*w4

	// 	acc2 = acc2*36409; //w2 <- w2*1/9
	// 	acc2 = acc2>>1; //w2 <- w2*1/2 ---> w2=w2/18

	// 	acc3= acc3-acc5;//w3 <- w3-w5
		
	// 	acc4 = acc4+acc2; //w4 <- w4+w2

	// 	acc4 = 0-acc4; //w4 <- -(w4+w2)

	// 	acc8 = acc2*30; //temp <- w2*30
	// 	acc6= acc8-acc6;//w6 <- 30*w2-w6

	// 	acc6 = acc6*61167; //w6 <- w6*1/15
	// 	acc6 = acc6>>2; //w6 <- w6*1/4 ---> w6=w6/60

	// 	acc2= acc2-acc6;//w2 <- w2-w6

	// 	result[i]= (result[i]+ acc7-acc3)&(SABER_Q-1);
	// 	result[(SABER_N/4)+i]= (result[(SABER_N/4)+i]+ acc6-acc2)&(SABER_Q-1);
	// 	result[2*(SABER_N/4)+i]= (result[2*(SABER_N/4)+i]+ acc5-acc1)&(SABER_Q-1);
	// 	// sprintf((char *)output, "%d, ", (int)result[i]);
	// 	// send_USART_str(output);
	// 	// sprintf((char *)output, "%d, ", (int)result[i+(SABER_N/4)]);
	// 	// send_USART_str(output);
	// 	// sprintf((char *)output, "%d, ", (int)result[i+2*(SABER_N/4)]);
	// 	// send_USART_str(output);
	// 	if(i>=64) {
	// 		result[3*(SABER_N/4)+i-SABER_N]= (result[3*(SABER_N/4)+i-SABER_N]- acc4)&(SABER_Q-1);
	// 		// sprintf((char *)output, "%d, ", (int)result[3*(SABER_N/4)+i-SABER_N]);
	// 		// send_USART_str(output);
	// 	} else {
	// 		result[3*(SABER_N/4)+i]= (result[3*(SABER_N/4)+i] + acc4)&(SABER_Q-1);
	// 		// sprintf((char *)output, "%d, ", (int)result[3*(SABER_N/4)+i]);
	// 		// send_USART_str(output);
	// 	}
	// }
//}


void toom_cook_4way_mem1(const uint16_t* a1,const uint16_t* b1, uint16_t* result){//the output is in a consecutive memory

	int16_t i;
	int16_t small_len=SABER_N/4;
	uint16_t result_final[2*SABER_N];

//----------------array declaration to hold results--------------------

	uint16_t w_m[896];

	int16_t acc1, acc2, acc3, acc4, acc5, acc6, acc7, acc8, acc9, acc10;

//----------------array declaration to hold results ends--------------------
	//--------------------these data are created for place holding---------


	uint16_t a2_ph[small_len],b2_ph[small_len];
	uint16_t a3_ph[small_len],b3_ph[small_len];
	uint16_t a4_ph[small_len],b4_ph[small_len];
	uint16_t a5_ph[small_len],b5_ph[small_len];
	uint16_t a6_ph[small_len],b6_ph[small_len];



	//--------------------these data are created for place holding ends---------

	int16_t inv3=43691;
	int16_t inv9=36409;
	int16_t inv15=61167;

	int16_t int45=45;
	int16_t int30=30;
	int16_t int0=0;

	for(i=0; i<2*SABER_N; i++){
		result_final[i] = 0;
	}

	for(i=0;i<small_len;i++){
		
		acc3=a1[small_len*0+i];
		acc4=a1[small_len*1+i];
		acc5=a1[small_len*2+i];
		acc6=a1[small_len*3+i];

		acc7=acc3+acc5;//a0+a2
		acc8=acc4+acc6;//a1+a3

		acc9=acc7+acc8;//a0+a1+a2+a3
		acc10=acc7-acc8;
		
		a3_ph[i]=acc9;
		a4_ph[i]=acc10;

		acc7=acc3*4+acc5;//4*a0+a2
		acc7=acc7*2;//8*a0+2*a2
		acc8=acc4*4+acc6;//4*a1+a3
		
		acc9=acc7+acc8;
		acc10=acc7-acc8;				

		a5_ph[i]=acc9;
		a6_ph[i]=acc10;

		acc7=2*acc6+acc5;
		acc7=acc7*2+acc4;	
		acc7=acc7*2+acc3;
		
		a2_ph[i]=acc7;
		//a1_ph[i]=acc6;
		//a7_ph[i]=acc3;

	}


	for(i=0;i<small_len;i++){
		
		acc3=b1[small_len*0+i];
		acc4=b1[small_len*1+i];
		acc5=b1[small_len*2+i];
		acc6=b1[small_len*3+i];

		acc7=acc3+acc5;//a0+a2
		acc8=acc4+acc6;//a1+a3

		acc9=acc7+acc8;//a0+a1+a2+a3
		acc10=acc7-acc8;
		
		b3_ph[i]=acc9;
		b4_ph[i]=acc10;

		acc7=acc3*4+acc5;//4*a0+a2
		acc7=acc7*2;//8*a0+2*a2
		acc8=acc4*4+acc6;//4*a1+a3
		
		acc9=acc7+acc8;
		acc10=acc7-acc8;				

		b5_ph[i]=acc9;
		b6_ph[i]=acc10;

		acc7=2*acc6+acc5;
		acc7=acc7*2+acc4;	
		acc7=acc7*2+acc3;
		
		b2_ph[i]=acc7;
		//b1_ph[i]=acc6;
		//b7_ph[i]=acc3;

	}


	memset(w_m,0,896*sizeof(int16_t));
//----------------------------------------------------------
	//memset(w7,0,2*small_len*sizeof(int16_t));
	//unrolled_kara_mem(a7_ph, b7_ph, w7, small_len/2);
	//unrolled_kara_mem(a1, b1, w7, small_len/2);
	unrolled_kara_mem(a1, b1, w_m+768, small_len/2);

	//karatsuba_simple1(a7_ph, b7_ph, w7);

//----------------------------------------------------------
	//memset(w6,0,2*small_len*sizeof(int16_t));
	//unrolled_kara_mem(a6_ph, b6_ph, w6, small_len/2);
	unrolled_kara_mem(a6_ph, b6_ph, w_m+640, small_len/2);

	//karatsuba_simple1(a6_ph, b6_ph, w6);

//----------------------------------------------------------
	//memset(w5,0,2*small_len*sizeof(int16_t));
	//unrolled_kara_mem(a5_ph, b5_ph, w5, small_len/2);
	unrolled_kara_mem(a5_ph, b5_ph, w_m+512, small_len/2);

	//karatsuba_simple1(a5_ph, b5_ph, w5);

//----------------------------------------------------------
	//memset(w4,0,2*small_len*sizeof(int16_t));
	//unrolled_kara_mem(a4_ph,b4_ph, w4, small_len/2);
	unrolled_kara_mem(a4_ph,b4_ph, w_m+384, small_len/2);

	//karatsuba_simple1(a4_ph,b4_ph, w4);

//----------------------------------------------------------
	//memset(w3,0,2*small_len*sizeof(int16_t));
	//unrolled_kara_mem(a3_ph, b3_ph, w3, small_len/2);
	unrolled_kara_mem(a3_ph, b3_ph, w_m+256, small_len/2);

	//karatsuba_simple1(a3_ph, b3_ph, w3);

//----------------------------------------------------------
	//memset(w2,0,2*small_len*sizeof(int16_t));
	//unrolled_kara_mem(a2_ph, b2_ph, w2, small_len/2);
	unrolled_kara_mem(a2_ph, b2_ph, w_m+128, small_len/2);

	//karatsuba_simple1(a2_ph, b2_ph, w2);

//----------------------------------------------------------
	//memset(w1,0,2*small_len*sizeof(int16_t));
	//unrolled_kara_mem(a1_ph, b1_ph, w1, small_len/2);
	//unrolled_kara_mem(a1+small_len*3, b1+small_len*3, w1, small_len/2);
	unrolled_kara_mem(a1+small_len*3, b1+small_len*3, w_m, small_len/2);

	//karatsuba_simple1(a1_ph, b1_ph, w1);

//----------------------------------------------------------

	/*************************MULTIPLICATION ENDS***********************************/
			
	for(i=0;i<2*small_len-1;i++){

		acc1=w_m[i];
		acc2=w_m[i+128];
		acc3=w_m[i+256];
		acc4=w_m[i+384];
		acc5=w_m[i+512];
		acc6=w_m[i+640];
		acc7=w_m[i+768];

		acc2= acc2+acc5;//w2 <- w2+w5
		acc6= acc6-acc5;// w6 <- w6-w5
		acc4= acc4-acc3;// w4 <- w4-w3
		
		acc5= acc5-acc1;// w5 <- w5-w1
		//temp = acc7<<6; //temp <- 64*w7
		acc8 = acc7<<6; //temp <- 64*w7

		acc5= acc5-acc8;// w5 <- w5-64*w7
		acc4 = acc4>>1; //w4 <- w4/2
		acc3 = acc3+acc4;//w3 <- w3+w4

		acc8 = acc5<<1; //temp <- 2*w5
		acc5= acc6+acc8;//w5 <- 2*w5+w6

		acc8 = acc3<<6; //temp <- 64*w3
		acc8 = acc3+acc8; //temp <- 65*w3
		acc2= acc2-acc8;// w2 <- w2-65*w3

		acc3= acc3-acc7;// w3 <- w3-w7

		acc3= acc3-acc1;// w3 <- w3-w1

		acc8 = acc3*int45; //temp <- 45*w3
		acc2 = acc2+acc8; //w2 <- w2+45*w3


		acc8 = acc3<<3; //temp <- 8*w3
		acc5= acc5-acc8;//w5 <- w5-8*w3

		acc5 = acc5*inv3; //w5 <- w5*1/3
		acc5 = acc5>>3; //w5 <- w5*1/8 ---> w5=w5/24

		acc6 = acc2+acc6; //w6 <- w6+w2
		acc8 = acc4<<4; //temp <- 16*w4
		acc2 = acc2+acc8; //w2 <- w2+16*w4

		acc2 = acc2*inv9; //w2 <- w2*1/9
		acc2 = acc2>>1; //w2 <- w2*1/2 ---> w2=w2/18

		acc3= acc3-acc5;//w3 <- w3-w5
		
		acc4 = acc4+acc2; //w4 <- w4+w2

		acc4 = int0-acc4; //w4 <- -(w4+w2)

		acc8 = acc2*int30; //temp <- w2*30
		acc6= acc8-acc6;//w6 <- 30*w2-w6

		acc6 = acc6*inv15; //w6 <- w6*1/15
		acc6 = acc6>>2; //w6 <- w6*1/4 ---> w6=w6/60

		acc2= acc2-acc6;//w2 <- w2-w6

		result_final[0*small_len+i]= (result_final[0*small_len+i]+ acc7);
		result_final[1*small_len+i]= (result_final[1*small_len+i]+ acc6);
		result_final[2*small_len+i]= (result_final[2*small_len+i]+ acc5);
		result_final[3*small_len+i]= (result_final[3*small_len+i]+ acc4);
		result_final[4*small_len+i]= (result_final[4*small_len+i]+ acc3);
		result_final[5*small_len+i]= (result_final[5*small_len+i]+ acc2);
		result_final[6*small_len+i]= (result_final[6*small_len+i]+ acc1);

	}

	  for(i=SABER_N;i<2*SABER_N-1;i++){
	    result_final[i-SABER_N]=(result_final[i-SABER_N]-result_final[i]);
	  }
	
	for(i=0;i<SABER_N;i++){
		result[i]=result_final[i]&(SABER_Q-1);
	}	

}


void toom_cook_4way(const uint16_t* a1,const uint16_t* b1, uint16_t* result_final){

	int16_t i;
	int16_t small_len=SABER_N/4;
	
	//char tempstr[50];


	/*uint64_t p_mod_or;	

	p_mod_or=p_mod;
	p_mod=p_mod*8;
	*/


//----------------array declaration to hold smaller arrays--------------------

	uint16_t a[small_len],b[small_len];
	uint16_t th_a[small_len],t_h_a[small_len];
	uint16_t th_b[small_len],t_h_b[small_len];
	uint16_t temp1[2*small_len];

//----------------array declaration to hold smaller arrays ends--------------------

	int16_t inv3=43691;
	int16_t inv9=36409;
	int16_t inv15=61167;

	int16_t int45=45;
	int16_t int30=30;
	int16_t int0=0;



//----------------array declaration to hold results--------------------
	uint16_t w1[2*small_len],w2[2*small_len],w3[2*small_len],w4[2*small_len],w5[2*small_len],w6[2*small_len],w7[2*small_len];

//----------------array declaration to hold results ends--------------------
	//--------------------these data are created for place holding---------

	uint16_t a1_ph[small_len],b1_ph[small_len];
	uint16_t a2_ph[small_len],b2_ph[small_len];
	uint16_t a3_ph[small_len],b3_ph[small_len];
	uint16_t a4_ph[small_len],b4_ph[small_len];
	uint16_t a5_ph[small_len],b5_ph[small_len];
	uint16_t a6_ph[small_len],b6_ph[small_len];

	//uint16_t temp_b[small_len];
	
	for(i=0;i<small_len;i++){

		a1_ph[i]=a1[i+0];
		b1_ph[i]=b1[i+0];

		th_a[i]= a1_ph[i]<<2; //th_x contains 4*x[0]
		th_b[i]= b1_ph[i]<<2;
	
		th_a[i]= th_a[i]+a1[small_len*2+i];//th_x contains 4*x[0]+x[2]
		th_b[i]= th_b[i]+b1[small_len*2+i];
	
		th_a[i]= th_a[i]<<1;//th_x_avx contains 8*x[0]+2*x[2]
		th_b[i]= th_b[i]<<1;
	
		t_h_a[i]= a1[small_len*1+i];//t_h_x_avx contains x[1]
		t_h_b[i]= b1[small_len*1+i];
	
		t_h_a[i]= t_h_a[i]<<2;//t_h_x_avx contains 4*x[1]
		t_h_b[i]= t_h_b[i]<<2;

	
		t_h_a[i]= t_h_a[i]+a1[small_len*3+i];//th_x_avx contains 4*x[1]+x[3]
		t_h_b[i]= t_h_b[i]+b1[small_len*3+i];

		a2_ph[i]= th_a[i]+t_h_a[i];//create th
		b2_ph[i]= th_b[i]+t_h_b[i];

		a3_ph[i]= th_a[i]-t_h_a[i];//create t_h
		b3_ph[i]= th_b[i]-t_h_b[i];

		
		th_a[i]= a1[small_len*2+i]+a1[small_len*0+i];//create partial sum for t_1 and t1
		th_b[i]= b1[small_len*2+i]+b1[small_len*0+i];//th_x_avx contains x[2]+x[0]
		
		t_h_a[i]= a1[small_len*3+i]+a1[small_len*1+i];//th_x_avx contains x[3]+x[1]
		t_h_b[i]= b1[small_len*3+i]+b1[small_len*1+i];
	
		a4_ph[i]= th_a[i]+t_h_a[i];// x[0]+x[1]+x[2]+x[3]
		b4_ph[i]= th_b[i]+t_h_b[i];

		a5_ph[i]= th_a[i]-t_h_a[i];//-x[3]+x[2]-x[1]+x[0]
		b5_ph[i]= th_b[i]-t_h_b[i];
	
		a6_ph[i]= a1[small_len*3+i];//x_avx contains x[3]
		b6_ph[i]= b1[small_len*3+i];//-------------------t_inf ends----------------------

	
		a[i]= a6_ph[i]+a1[small_len*3+i];// 2*x[3]
		b[i]= b6_ph[i]+b1[small_len*3+i];

		a[i]= a[i]+a1[small_len*2+i];// 2*x[3]+x[2]
		b[i]= b[i]+b1[small_len*2+i];
		
		a[i]= a[i]<<1;// 4*x[3]+2*x[2]
		b[i]= b[i]<<1;
		
		a[i]= a[i]+a1[small_len*1+i];// 4*x[3]+2*x[2]+x[1]
		b[i]= b[i]+b1[small_len*1+i];
		
		a[i]= a[i]<<1;// 8*x[3]+4*x[2]+2*x[1]
		b[i]= b[i]<<1;
		
		a[i]= a[i]+a1[small_len*0+i];// 8*x[3]+8*x[2]+2*x[1]+x[0]
		b[i]= b[i]+b1[small_len*0+i];
	}


	//memset(temp_b,0,small_len*sizeof(int16_t));
	memset(w7,0,2*small_len*sizeof(int16_t));

	//kara_mem(a1_ph, temp_b, b1_ph, w7, small_len/2);
	unrolled_kara_mem(a1_ph, b1_ph, w7, small_len/2);

	//memset(temp_b,0,small_len*sizeof(int16_t));
	memset(w5,0,2*small_len*sizeof(int16_t));

	//kara_mem(a2_ph, temp_b, b2_ph, w5, small_len/2);
	unrolled_kara_mem(a2_ph, b2_ph, w5, small_len/2);

	//memset(temp_b,0,small_len*sizeof(int16_t));
	memset(w6,0,2*small_len*sizeof(int16_t));

	//kara_mem(a3_ph, temp_b, b3_ph, w6, small_len/2);
	unrolled_kara_mem(a3_ph, b3_ph, w6, small_len/2);

	//memset(temp_b,0,small_len*sizeof(int16_t));
	memset(w3,0,2*small_len*sizeof(int16_t));

	//kara_mem(a4_ph, temp_b, b4_ph, w3, small_len/2);
	unrolled_kara_mem(a4_ph,b4_ph, w3, small_len/2);

	//memset(temp_b,0,small_len*sizeof(int16_t));
	memset(w4,0,2*small_len*sizeof(int16_t));

	//kara_mem(a5_ph, temp_b, b5_ph, w4, small_len/2);
	unrolled_kara_mem(a5_ph, b5_ph, w4, small_len/2);

	//memset(temp_b,0,small_len*sizeof(int16_t));
	memset(w1,0,2*small_len*sizeof(int16_t));

	//kara_mem(a6_ph, temp_b, b6_ph, w1, small_len/2);
	unrolled_kara_mem(a6_ph, b6_ph, w1, small_len/2);

	//memset(temp_b,0,small_len*sizeof(int16_t));
	memset(w2,0,2*small_len*sizeof(int16_t));

	//kara_mem(a, temp_b, b, w2, small_len/2);
	unrolled_kara_mem(a, b, w2, small_len/2);


	/*	--------------------------------------------
		---------------Solution starts--------------
		--------------------------------------------

	*/
	

	for(i=0;i<2*small_len;i++){
			
		w2[i]= w2[i]+w5[i];//w2 <- w2+w5
		w6[i]= w6[i]-w5[i];// w6 <- w6-w5
		w4[i]= w4[i]-w3[i];// w4 <- w4-w3
		
		w5[i]= w5[i]-w1[i];// w5 <- w5-w1
		temp1[i] = w7[i]<<6; //temp <- 64*w7
		w5[i]= w5[i]-temp1[i];// w5 <- w5-64*w7

		w4[i] = w4[i]>>1; //w4 <- w4/2
		w3[i] = w3[i]+w4[i];//w3 <- w3+w4

		temp1[i] = w5[i]<<1; //temp <- 2*w5
		w5[i]= w6[i]+temp1[i];//w5 <- 2*w5+w6

		temp1[i] = w3[i]<<6; //temp <- 64*w3
		temp1[i] = w3[i]+temp1[i]; //temp <- 65*w3
		w2[i]= w2[i]-temp1[i];// w2 <- w2-65*w3

		w3[i]= w3[i]-w7[i];// w3 <- w3-w7
		w3[i]= w3[i]-w1[i];// w3 <- w3-w1

		temp1[i] = w3[i]*int45; //temp <- 45*w3
		w2[i] = w2[i]+temp1[i]; //w2 <- w2+45*w3

		temp1[i] = w3[i]<<3; //temp <- 8*w3
		w5[i]= w5[i]-temp1[i];//w5 <- w5-8*w3

		w5[i] = w5[i]*inv3; //w5 <- w5*1/3
		w5[i] = w5[i]>>3; //w5 <- w5*1/8 ---> w5=w5/24

		w6[i] = w2[i]+w6[i]; //w6 <- w6+w2
		temp1[i] = w4[i]<<4; //temp <- 16*w4
		w2[i] = w2[i]+temp1[i]; //w2 <- w2+16*w4

		w2[i] = w2[i]*inv9; //w2 <- w2*1/9
		w2[i] = w2[i]>>1; //w2 <- w2*1/2 ---> w2=w2/18

		w3[i]= w3[i]-w5[i];//w3 <- w3-w5
		
		w4[i] = w4[i]+w2[i]; //w4 <- w4+w2

		w4[i] = int0-w4[i]; //w4 <- -(w4+w2)

		temp1[i] = w2[i]*int30; //temp <- w2*30
		w6[i]= temp1[i]-w6[i];//w6 <- 30*w2-w6

		w6[i] = w6[i]*inv15; //w6 <- w6*1/15
		w6[i] = w6[i]>>2; //w6 <- w6*1/4 ---> w6=w6/60

		w2[i]= w2[i]-w6[i];//w2 <- w2-w6

	}


	for(i=0; i<2*SABER_N; i++){
		result_final[i] = 0;
	}	

	for(i=0;i<2*small_len-1;i++){	
		result_final[0*small_len+i]= (result_final[0*small_len+i]+ w7[i])&(SABER_Q-1);
		result_final[1*small_len+i]= (result_final[1*small_len+i]+ w6[i])&(SABER_Q-1);
		result_final[2*small_len+i]= (result_final[2*small_len+i]+ w5[i])&(SABER_Q-1);
		result_final[3*small_len+i]= (result_final[3*small_len+i]+ w4[i])&(SABER_Q-1);
		result_final[4*small_len+i]= (result_final[4*small_len+i]+ w3[i])&(SABER_Q-1);
		result_final[5*small_len+i]= (result_final[5*small_len+i]+ w2[i])&(SABER_Q-1);
		result_final[6*small_len+i]= (result_final[6*small_len+i]+ w1[i])&(SABER_Q-1);
	}		
	


}



// void karatsuba_unroll(const uint16_t *a, const uint16_t *b, uint16_t *result){

// 	int16_t i;

// 	uint16_t result_m[468];

// 	karatsuba_asm(a, b, result, result_m);

// 	// for(i=189;i<468;i++){
// 	// 	result_m[i]=0;
// 	// }

// //loop1	
// 	// for(i=0;i<16;i++){
// 	// 	result_m[i]=a[i]+a[i+16];
// 	// 	result_m[i+16]=a[i]+a[i+32];
// 	// 	result_m[i+32]=a[i+16]+a[i+48];
// 	// 	result_m[i+48]=a[i+32]+a[i+48];
// 	// 	result_m[i+64]=result_m[i+16]+result_m[i+32];
// 	// // }


// 	// // for(i=0;i<16;i++){
// 	// 	result_m[i+80]=b[i]+b[i+16];
// 	// 	result_m[i+96]=b[i]+b[i+32];
// 	// 	result_m[i+112]=b[i+16]+b[i+48];
// 	// 	result_m[i+128]=b[i+32]+b[i+48];
// 	// 	result_m[i+144]=result_m[i+96]+result_m[i+112];
// 	// }
	
// 	// school_book_mul2_16(a, b, result_m+189);
// 	// school_book_mul2_16(a+16, b+16, result_m+31+189);
// 	// school_book_mul2_16(a+32, b+32, result_m+62+189);
// 	// school_book_mul2_16(a+48, b+48, result_m+93+189);
// 	// school_book_mul2_16(result_m, result_m+80, result_m+124+189);
// 	// school_book_mul2_16(result_m+16, result_m+96, result_m+155+189);
// 	// school_book_mul2_16(result_m+32, result_m+112, result_m+186+189);
// 	// school_book_mul2_16(result_m+48, result_m+128, result_m+217+189);
// 	// school_book_mul2_16(result_m+64, result_m+144, result_m+248+189);
	

// 	// for(i=0;i<189;i++){
// 	// 	result_m[i]=0;
// 	// }

// //-------------------------intermediate merging---------------------
// //------------------------------------------------
// 	// for(i=0;i<31;i++){	
// 	// 	result_m[i]=result_m[i+189];
// 	// 	result_m[i+63]=result_m[i+155+189];
// 	// 	result_m[i+126]=result_m[i+62+189];
// 	// }	

// 	// for(i=0;i<31;i++){	
// 	// 	result_m[i+16]=(result_m[i+16]+result_m[i+124+189]-result_m[i+189]-result_m[i+31+189]);
// 	// 	result_m[i+79]=(result_m[i+79]+result_m[i+248+189]-result_m[i+155+189]-result_m[i+186+189]);
// 	// 	result_m[i+142]=(result_m[i+142]+result_m[i+217+189]-result_m[i+62+189]-result_m[i+93+189]);				
// 	// }


// 	// for(i=0;i<31;i++){	
// 	// 	result_m[i+32]=(result_m[i+32]+result_m[i+31+189]);
// 	// 	result_m[i+95]=(result_m[i+95]+result_m[i+186+189]);
// 	// 	result_m[i+158]=(result_m[i+158]+result_m[i+93+189]);
// 	// }

// //-------------------------final merging----------------------------
// 	// for(i=0;i<63;i++){	
// 	// 	result[i]=result_m[i];
// 	// }

	
// 	// for(i=0;i<63;i++){	
// 	// 	result[i+32]=(result[i+32]+result_m[i+63]-result_m[i]-result_m[i+126]);		
// 	// // }


// 	// // for(i=0;i<63;i++){	
// 	// 	result[i+64]=(result[i+64]+result_m[i+126]);

// 	// }

// }
/*
void toom_cook_4way_mem(const uint16_t* a1,const uint16_t* b1, uint16_t* result)
{//the output is in a consecutive memory. Does not reuse memory for speed. Suitable for m4
	uint16_t w_m[896+640];

	int16_t inv3=43691;
	int16_t inv9=36409;
	int16_t inv15=61167;

	int16_t int45=45;
	int16_t int30=30;
	int16_t int0=0;

	int16_t i;
	int16_t acc1, acc2, acc3, acc4, acc5, acc6, acc7, acc8, acc9, acc10;

	for(i=0;i<(SABER_N/4);i++){
		
		acc3=a1[(SABER_N/4)*0+i];
		acc4=a1[(SABER_N/4)*1+i];
		acc5=a1[(SABER_N/4)*2+i];
		acc6=a1[(SABER_N/4)*3+i];

		acc7=acc3+acc5;//a0+a2
		acc8=acc4+acc6;//a1+a3

		acc9=acc7+acc8;//a0+a1+a2+a3
		acc10=acc7-acc8;

		w_m[i+64+896]=acc9;
		w_m[i+128+896]=acc10;

		acc7=acc3*4+acc5;//4*a0+a2
		acc7=acc7*2;//8*a0+2*a2
		acc8=acc4*4+acc6;//4*a1+a3

		acc9=acc7+acc8;
		acc10=acc7-acc8;				

		w_m[i+192+896]=acc9;
		w_m[i+256+896]=acc10;

		acc7=2*acc6+acc5;
		acc7=acc7*2+acc4;	
		acc7=acc7*2+acc3;

		w_m[i+896]=acc7;
	}

	for(i=0;i<(SABER_N/4);i++){
		
		acc3=b1[(SABER_N/4)*0+i];
		acc4=b1[(SABER_N/4)*1+i];
		acc5=b1[(SABER_N/4)*2+i];
		acc6=b1[(SABER_N/4)*3+i];

		acc7=acc3+acc5;//a0+a2
		acc8=acc4+acc6;//a1+a3

		acc9=acc7+acc8;//a0+a1+a2+a3
		acc10=acc7-acc8;
		
		w_m[i+64+320+896]=acc9;
		w_m[i+128+320+896]=acc10;

		acc7=acc3*4+acc5;//4*a0+a2
		acc7=acc7*2;//8*a0+2*a2
		acc8=acc4*4+acc6;//4*a1+a3

		acc9=acc7+acc8;
		acc10=acc7-acc8;				

		w_m[i+192+320+896]=acc9;
		w_m[i+256+320+896]=acc10;

		acc7=2*acc6+acc5;
		acc7=acc7*2+acc4;	
		acc7=acc7*2+acc3;
		
		w_m[i+320+896]=acc7;
	}

	memset(w_m,0,(896)*sizeof(int16_t));
	// unrolled_kara_mem(a1, b1, w_m+768, (SABER_N/4)/2);
	// unrolled_kara_mem(w_m+256+896, w_m+256+320+896, w_m+640, (SABER_N/4)/2);
	// unrolled_kara_mem(w_m+192+896, w_m+192+320+896, w_m+512, (SABER_N/4)/2);
	// unrolled_kara_mem(w_m+128+896,w_m+128+320+896, w_m+384, (SABER_N/4)/2);
	// unrolled_kara_mem(w_m+64+896, w_m+64+320+896, w_m+256, (SABER_N/4)/2);
	// unrolled_kara_mem(w_m+896, w_m+320+896, w_m+128, (SABER_N/4)/2);
	// unrolled_kara_mem(a1+(SABER_N/4)*3, b1+(SABER_N/4)*3, w_m, (SABER_N/4)/2);

	karatsuba_unroll(a1, b1, w_m+768);
	karatsuba_unroll(w_m+256+896, w_m+256+320+896, w_m+640);
	karatsuba_unroll(w_m+192+896, w_m+192+320+896, w_m+512);
	karatsuba_unroll(w_m+128+896,w_m+128+320+896, w_m+384);
	karatsuba_unroll(w_m+64+896, w_m+64+320+896, w_m+256);
	karatsuba_unroll(w_m+896, w_m+320+896, w_m+128);
	karatsuba_unroll(a1+(SABER_N/4)*3, b1+(SABER_N/4)*3, w_m);

	for(i=0;i<2*(SABER_N/4)-1;i++){

		acc1=w_m[i];
		acc2=w_m[i+128];
		acc3=w_m[i+256];
		acc4=w_m[i+384];
		acc5=w_m[i+512];
		acc6=w_m[i+640];
		acc7=w_m[i+768];

		// sprintf((char *)output, "%d, ", (int)w_m[i]);
		// send_USART_str(output);
		// sprintf((char *)output, "%d, ", (int)w_m[i+128]);
		// send_USART_str(output);
		// sprintf((char *)output, "%d, ", (int)w_m[i+768]);
		// send_USART_str(output);

		acc2= acc2+acc5;//w2 <- w2+w5
		acc6= acc6-acc5;// w6 <- w6-w5
		acc4= acc4-acc3;// w4 <- w4-w3
		
		acc5= acc5-acc1;// w5 <- w5-w1
		//temp = acc7<<6; //temp <- 64*w7
		acc8 = acc7<<6; //temp <- 64*w7

		acc5= acc5-acc8;// w5 <- w5-64*w7
		acc4 = acc4>>1; //w4 <- w4/2
		acc3 = acc3+acc4;//w3 <- w3+w4

		acc8 = acc5<<1; //temp <- 2*w5
		acc5= acc6+acc8;//w5 <- 2*w5+w6

		acc8 = acc3<<6; //temp <- 64*w3
		acc8 = acc3+acc8; //temp <- 65*w3
		acc2= acc2-acc8;// w2 <- w2-65*w3

		acc3= acc3-acc7;// w3 <- w3-w7

		acc3= acc3-acc1;// w3 <- w3-w1

		acc8 = acc3*45; //temp <- 45*w3
		acc2 = acc2+acc8; //w2 <- w2+45*w3

		acc8 = acc3<<3; //temp <- 8*w3
		acc5= acc5-acc8;//w5 <- w5-8*w3

		acc5 = acc5*43691; //w5 <- w5*1/3
		acc5 = acc5>>3; //w5 <- w5*1/8 ---> w5=w5/24

		acc6 = acc2+acc6; //w6 <- w6+w2
		acc8 = acc4<<4; //temp <- 16*w4
		acc2 = acc2+acc8; //w2 <- w2+16*w4

		acc2 = acc2*36409; //w2 <- w2*1/9
		acc2 = acc2>>1; //w2 <- w2*1/2 ---> w2=w2/18

		acc3= acc3-acc5;//w3 <- w3-w5
		
		acc4 = acc4+acc2; //w4 <- w4+w2

		acc4 = 0-acc4; //w4 <- -(w4+w2)

		acc8 = acc2*30; //temp <- w2*30
		acc6= acc8-acc6;//w6 <- 30*w2-w6

		acc6 = acc6*61167; //w6 <- w6*1/15
		acc6 = acc6>>2; //w6 <- w6*1/4 ---> w6=w6/60

		acc2= acc2-acc6;//w2 <- w2-w6

		result[i]= (result[i]+ acc7-acc3)&(SABER_Q-1);
		result[(SABER_N/4)+i]= (result[(SABER_N/4)+i]+ acc6-acc2)&(SABER_Q-1);
		result[2*(SABER_N/4)+i]= (result[2*(SABER_N/4)+i]+ acc5-acc1)&(SABER_Q-1);
		// sprintf((char *)output, "%d, ", (int)result[i]);
		// send_USART_str(output);
		// sprintf((char *)output, "%d, ", (int)result[i+(SABER_N/4)]);
		// send_USART_str(output);
		// sprintf((char *)output, "%d, ", (int)result[i+2*(SABER_N/4)]);
		// send_USART_str(output);
		if(i>=64) {
			result[3*(SABER_N/4)+i-SABER_N]= (result[3*(SABER_N/4)+i-SABER_N]- acc4)&(SABER_Q-1);
			// sprintf((char *)output, "%d, ", (int)result[3*(SABER_N/4)+i-SABER_N]);
			// send_USART_str(output);
		} else {
			result[3*(SABER_N/4)+i]= (result[3*(SABER_N/4)+i] + acc4)&(SABER_Q-1);
			// sprintf((char *)output, "%d, ", (int)result[3*(SABER_N/4)+i]);
			// send_USART_str(output);
		}
	}
}/**/