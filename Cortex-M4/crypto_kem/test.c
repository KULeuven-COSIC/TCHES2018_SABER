#include "api.h"
#include "SABER_params.h"
#include "randombytes.h"
#include "stm32wrapper.h"
#include <string.h>

//XXX - temporary patch for testing m0mem
#include "poly_mul.c"

#define NTESTS 2

#define BENCHMARK

// You can paste the name of your function from poly_mul.c here
// BEWARE THAT THE PROTOTYPE MATCHES!!!
//pol_mul_test(uint16_t* a, uint16_t* b, uint16_t* res, uint16_t p, uint32_t n,uint32_t start);
#define pol_mul_test pol_mul

#define REDUCTION
// by commenting REDUCTION you use schoolbook multiplication returning a 2*N polynomial

//#define TOOM_COOK
//TOOM_COOK only reduces by modulus. Not by polynomial.

#define KEM_TEST

static int32_t compare_poly(uint16_t* a, uint16_t* b, uint32_t n);
static void pol_mul_ref(uint16_t* a, uint16_t* b, uint16_t* res, uint16_t p, uint32_t n,uint32_t start);
static void get_poly(uint16_t *a, uint64_t p,int32_t n);
void printarr(uint16_t * a, int len);
#ifdef BENCHMARK
static void printcycles(const char *s, unsigned int c);
#endif


static void test_multiplication(void)
{
  send_USART_str("\nStart test_mul\n");
#ifdef BENCHMARK
  unsigned int t0, t1;
#endif

  int i;

  uint16_t a1[SABER_N] = {2, 2, 2, 2};
  uint16_t b1[SABER_N] = {1,-1, 1, 1};
  uint16_t result_test[2*SABER_N-1] = {0}, result_ref[2*SABER_N-1] = {0};

  for (i = 0; i < NTESTS; ++i)
  {
#ifndef REDUCTION
    memset(result_ref, 0, sizeof(result_ref));
    memset(result_test, 0, sizeof(result_test));
#endif

    get_poly(a1,SABER_Q,SABER_N);
    get_poly(b1,SABER_Q,SABER_N);

#ifdef BENCHMARK
    t0 = DWT_CYCCNT;
#endif

    pol_mul_ref(a1, b1, result_ref, SABER_Q, SABER_N, 0);

#ifdef BENCHMARK
    t1 = DWT_CYCCNT;
    printcycles("ref mul cycles:", t1-t0);
#endif

//####################################################################
//--------------- CALL YOUR MULTIPLICATION STARTS HERE ---------------
//####################################################################
#ifdef BENCHMARK
    t0 = DWT_CYCCNT;
#endif

  #ifdef REDUCTION
    //pol_mul_test(a1, b1, result_test, SABER_Q, SABER_N, 0);
    //toom_cook_4way_mem(a1,b1,result_test,SABER_Q,SABER_N);
	memset(result_test, 0, sizeof(result_test));
    toom_cook_4way_mem(a1,b1,result_test);

  #else
    //pol_mul_schb(a1, b1, result_test, SABER_N);
    //unrolled_kara_mem(a1,b1,result_test,SABER_N/2);
    //toom_cook_4way(a1,b1,result_test,SABER_Q,SABER_N);
    //toom_cook_4way_mem(a1,b1,result_test,SABER_Q,SABER_N);
    
  #endif
#ifdef BENCHMARK
    t1 = DWT_CYCCNT;
    printcycles("opt mul cycles:", t1-t0);
#endif

#ifdef TOOM_COOK
    int j;
    for(j=0;j<(2*SABER_N-1);j++){
	result_ref[j]=result_ref[j]&(SABER_Q-1);
    }
#endif

//##################################################################
//--------------- CALL YOUR MULTIPLICATION ENDS HERE ---------------
//##################################################################
//printarr(&result_ref[254], 8);
//printarr(&result_test[254], 8);
#ifdef REDUCTION
	int j;
	for(j=SABER_N;j<2*SABER_N-1;j++){
		result_test[j]=0;	
	}

#endif

    if(compare_poly(result_ref,result_test,2*SABER_N-1)==1){
      send_USART_str("\n-----ERROR-----\n");
      
    }
  }
}

static int test_keys(void)
{
  unsigned char key_a[CRYPTO_BYTES], key_b[CRYPTO_BYTES];
  unsigned char pk[CRYPTO_PUBLICKEYBYTES];
  unsigned char sendb[CRYPTO_CIPHERTEXTBYTES];
  unsigned char sk_a[CRYPTO_SECRETKEYBYTES];
  int i;

  for(i=0; i<NTESTS; i++)
  {
    //Alice generates a public key
    crypto_kem_keypair(pk, sk_a);
    send_USART_str("PASS key pair generation!");

    //Bob derives a secret key and creates a response
    crypto_kem_enc(sendb, key_b, pk);
    send_USART_str("PASS encapsulation!");

    //Alice uses Bobs response to get her secret key
    crypto_kem_dec(key_a, sendb, sk_a);
    send_USART_str("PASS decapsulation!");

    if(memcmp(key_a, key_b, CRYPTO_BYTES))
      send_USART_str("ERROR KEYS\n");
  }

  return 0;
}


static int test_invalid_sk_a(void)
{
  unsigned char sk_a[CRYPTO_SECRETKEYBYTES];
  unsigned char key_a[CRYPTO_BYTES], key_b[CRYPTO_BYTES];
  unsigned char pk[CRYPTO_PUBLICKEYBYTES];
  unsigned char sendb[CRYPTO_CIPHERTEXTBYTES];
  int i;

  for(i=0; i<NTESTS; i++)
  {
    //Alice generates a public key
    crypto_kem_keypair(pk, sk_a);

    //Bob derives a secret key and creates a response
    crypto_kem_enc(sendb, key_b, pk);

    //Replace secret key with random values
    randombytes(sk_a, CRYPTO_SECRETKEYBYTES);

    //Alice uses Bobs response to get her secre key
    crypto_kem_dec(key_a, sendb, sk_a);

    if(!memcmp(key_a, key_b, CRYPTO_BYTES))
      send_USART_str("ERROR invalid sk_a\n");
  }

  return 0;
}


static int test_invalid_ciphertext(void)
{
  unsigned char sk_a[CRYPTO_SECRETKEYBYTES];
  unsigned char key_a[CRYPTO_BYTES], key_b[CRYPTO_BYTES];
  unsigned char pk[CRYPTO_PUBLICKEYBYTES];
  unsigned char sendb[CRYPTO_CIPHERTEXTBYTES];
  int i;
  size_t pos;

  for(i=0; i<NTESTS; i++)
  {
    randombytes((unsigned char *)&pos, sizeof(size_t));

    //Alice generates a public key
    crypto_kem_keypair(pk, sk_a);

    //Bob derives a secret key and creates a response
    crypto_kem_enc(sendb, key_b, pk);

    //Change some byte in the ciphertext (i.e., encapsulated key)
    sendb[pos % CRYPTO_CIPHERTEXTBYTES] ^= 23;

    //Alice uses Bobs response to get her secret key
    crypto_kem_dec(key_a, sendb, sk_a);

    if(!memcmp(key_a, key_b, CRYPTO_BYTES))
      send_USART_str("ERROR invalid ciphertext\n");
  }

  return 0;
}

static void get_poly(uint16_t *a, uint64_t p,int32_t n)
{
  int32_t i, r;
  for(i=0; i<n; i += 2){
    r = rng_get_random_blocking();
    a[i] = r & (p-1);
    a[i+1] = (r>>16) & (p-1);
  }
}

static int32_t compare_poly(uint16_t* a, uint16_t* b, uint32_t n)
{ // Comparing two polynomials over GF(p), a[x]=b[x]? : (0) a=b, (1) a!=b
  // SECURITY NOTE: TO BE USED FOR TESTING ONLY.
    uint32_t i;

    char output[32];
    for (i = 0; i < n; i++)
    {
        if (a[i] != b[i]) {
	sprintf((char *)output, "%d, ", (int)i);
	send_USART_str(output);
	sprintf((char *)output, "%d, ", (int)a[i]);
	send_USART_str(output);
	sprintf((char *)output, "%d, ", (int)b[i]);
	send_USART_str(output);
            return 1;
        }
    }

    return 0; 
}

void
printarr(uint16_t * a, int len)
{
  int i;
  char output[32];
  send_USART_str("arr = [");
  for (i = 0; i < len; i++) {
    sprintf((char *)output, "%d, ", (int)a[i]);
    send_USART_str(output);
  }
  send_USART_str("];\n");
}

#ifdef BENCHMARK
static void printcycles(const char *s, unsigned int c)
{
  char outs[12];
  send_USART_str(s);
  sprintf(outs,"%d\n",c);
  send_USART_str(outs);
}
#endif

static void pol_mul_ref(uint16_t* a, uint16_t* b, uint16_t* res, uint16_t p, uint32_t n,uint32_t start)
{
#ifndef REDUCTION
  int i,j;
  for (i = 0; i < n; i++) {
    for (j = 0; j < n; j++) {
      res[i+j] += a[i]*b[j];
    }
  }
#else
  // Polynomial multiplication using the schoolbook method, c[x] = a[x]*b[x] 
  // SECURITY NOTE: TO BE USED FOR TESTING ONLY.  
  uint32_t i, j;

//-------------------normal multiplication-----------------
  uint16_t c[2*SABER_N];
  for (i = 0; i < 2*n; i++) c[i] = 0;

  for (i = start; i < start+n; i++) {
    for (j = start; j < start+n; j++) {
      c[i+j-2*start]=( c[i+j-2*start] + (a[i] * b[j]) )&(p-1);
    }
  }

  //---------------reduction-------
  for(i=n;i<2*n-1;i++){
    c[i-n]=(c[i-n]-c[i])&(p-1);
  }
  
  c[n]=0; //256th coefficient=0;
//----------------------copy result back----------------
  for(i=0; i<n; i++){
            res[i] = c[i];
  }
#endif
}

int main(void)
{
#ifndef BENCHMARK
  clock_setup(CLOCK_FAST);
#else
  clock_setup(CLOCK_BENCHMARK);
#endif
  gpio_setup();
  usart_setup(115200);
#ifdef BENCHMARK
  cyccnt_setup();
#endif
  rng_enable();

#ifdef KEM_TEST
  test_keys();
  test_invalid_sk_a();
  test_invalid_ciphertext();
#endif

  test_multiplication();

  send_USART_str("\ndone!");

// send # to terminate the screen script
  send_USART_str("#");

  while(1);

  return 0;
}
