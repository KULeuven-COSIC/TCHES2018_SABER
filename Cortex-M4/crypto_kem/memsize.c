#include "api.h"
#include "stm32wrapper.h"

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>



#define MAXSTACK 50000
unsigned char key_a[CRYPTO_BYTES], key_b[CRYPTO_BYTES];
unsigned char pk[CRYPTO_PUBLICKEYBYTES];
unsigned char sendb[CRYPTO_CIPHERTEXTBYTES];
unsigned char sk_a[CRYPTO_SECRETKEYBYTES];
  

char output[32];
unsigned int ctr;
unsigned char canary;
volatile unsigned char *p;
extern unsigned char _end; 

static unsigned int stack_count(unsigned char canaryv,volatile unsigned char *a)
{
  //volatile unsigned char *p = (a-MAXSTACK);
  p = (a-MAXSTACK);
  unsigned int c = 0;
  while(*p == canaryv && p < a)
  {
    p++;
    c++;
  }
  return c;
} 

#define WRITE_CANARY(X) {p=X;while(p>= (X-MAXSTACK)) *(p--) = canary;}
 
int main(void)
{
    clock_setup(CLOCK_FAST);
    gpio_setup();
    usart_setup(115200);
    rng_enable();


  volatile unsigned char a; /* Mark the beginning of the stack */
  //int i;
  //poly sk;
    canary = 42;

    WRITE_CANARY(&a);
  crypto_kem_keypair(pk, sk_a);
    ctr = MAXSTACK - stack_count(canary,&a);
  sprintf((char *)output, "RAM usage of keygen: %d",ctr);
    send_USART_str(output);

    WRITE_CANARY(&a);
  crypto_kem_keypair(pk, sk_a);
    ctr = MAXSTACK - stack_count(canary,&a);
  sprintf((char *)output, "RAM usage of keygen: %d",ctr);
    send_USART_str(output);
  
  WRITE_CANARY(&a);
  crypto_kem_enc(sendb, key_b, pk);
    ctr = MAXSTACK - stack_count(canary,&a);
  sprintf((char *)output, "RAM usage of sharedb: %d",ctr);
    send_USART_str(output);
     


  WRITE_CANARY(&a);   
  crypto_kem_dec(key_a, sendb, sk_a);
    ctr = MAXSTACK - stack_count(canary,&a);
  sprintf((char *)output, "RAM usage of shareda: %d",ctr);
    send_USART_str(output);
        
  WRITE_CANARY(&a);
  crypto_kem_keypair(pk, sk_a);
  crypto_kem_enc(sendb, key_b, pk);
    crypto_kem_dec(key_a, sendb, sk_a);
    ctr = MAXSTACK - stack_count(canary,&a);
  sprintf((char *)output, "RAM usage of KEM: %d",ctr);
    send_USART_str(output);
  
    sprintf((char *)output, "done!");
    send_USART_str(output);

// send # to terminate the screen script
    send_USART_str("#");

    return 0;
}
