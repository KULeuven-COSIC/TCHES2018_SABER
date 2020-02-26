#include "api.h"
#include "stm32wrapper.h"

#include <stdio.h>
#include <stdint.h>
#include <string.h>

static void printcycles(const char *s, unsigned int c)
{
  char outs[12];
  send_USART_str(s);
  sprintf(outs,"%d\n",c);
  send_USART_str(outs);
}


int main(void)
{
  unsigned char ss[CRYPTO_BYTES];
  unsigned char sk[CRYPTO_SECRETKEYBYTES];
  unsigned char pk[CRYPTO_PUBLICKEYBYTES];
  unsigned char ct[CRYPTO_CIPHERTEXTBYTES];
  unsigned int t0, t1;

  clock_setup(CLOCK_BENCHMARK);
  gpio_setup();
  usart_setup(115200);
  cyccnt_setup();
  rng_enable();

  send_USART_str("==========================");

  // Key-pair generation
  t0 = DWT_CYCCNT;
  crypto_kem_keypair(pk, sk);
  t1 = DWT_CYCCNT;
  printcycles("keypair cycles:", t1-t0);

  // Encapsulation
  t0 = DWT_CYCCNT;
  crypto_kem_enc(ct, ss, pk);
  t1 = DWT_CYCCNT;
  printcycles("encaps cycles: ", t1-t0);

  // Decapsulation
  t0 = DWT_CYCCNT;
  crypto_kem_dec(ss, ct, sk);
  t1 = DWT_CYCCNT;
  printcycles("decaps cycles: ", t1-t0);

  send_USART_str("#");
  return 0;
}
