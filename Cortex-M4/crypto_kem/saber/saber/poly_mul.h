#include "SABER_params.h"

void pol_mul(uint16_t* a, uint16_t* b, uint16_t* res, uint16_t p, uint32_t n);

// void toom_cook_4way(const uint16_t* a1,const uint16_t* b1, uint16_t* result_final, uint64_t p_mod, uint16_t n);
// void toom_cook_4way_mem(const uint16_t* a1,const uint16_t* b1, uint16_t* result_final, uint64_t p_mod, uint16_t n);
void toom_cook_4way(const uint16_t* a1,const uint16_t* b1, uint16_t* result_final);

//void toom_cook_4way_mem(const uint16_t* a1,const uint16_t* b1, uint16_t* result_final);
#define toom_cook_4way_mem toom_cook_4way_mem_asm

extern void pol_mul_schb(const uint16_t * a, const uint16_t * b, uint16_t * c, const uint16_t len);

extern void school_book_mul1(const uint16_t* a, const uint16_t* b, uint16_t *c, uint16_t len);

extern void school_book_mul2(const uint16_t* a, const uint16_t* b, uint16_t *c, uint16_t len);

extern void school_book_mul2_16(const uint16_t* a, const uint16_t* b, uint16_t *c);

extern void unrolled_kara_mem(const uint16_t * a, const uint16_t * c, uint16_t * d, const uint16_t k);

//void karatsuba_unroll(const uint16_t *a, const uint16_t *b, uint16_t *result);
#define karatsuba_unroll karatsuba_asm

extern void toom_cook_4way_mem_asm(const uint16_t* a1,const uint16_t* b1, uint16_t* result);

extern void karatsuba_asm(const uint16_t *a, const uint16_t *b, uint16_t *result);
