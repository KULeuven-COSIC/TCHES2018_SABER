	.thumb
	.syntax unified
.text

/*************************************************
* void pol_mul_schb(const uint16_t * a, const uint16_t * b, uint16_t * c, const uint16_t len);
**************************************************/
.global pol_mul_schb
.func pol_mul_schb, pol_mul_schb
.type pol_mul_schb, %function
pol_mul_schb:
	cbz		r3, exit_pol_mul_schb
	stmdb	sp!, {r4, r5, r6, r7, r8, r9, lr}
	lsl 	r3, r3, #1
	mov		r4, r2
	add		r2, r2, r3
outloop_pol_mul_schb:
	mov		r5, #0
	ldrh	r6, [r0], #2
innloop_pol_mul_schb:
	ldrh	r8, [r4], #2
	ldr		r7, [r1, r5]
	ldrh	r9, [r4], #-2
	add		r5, r5, #4
	smlabb	r8, r6, r7, r8
	smlabt	r9, r6, r7, r9
	strh	r8, [r4], #2
	strh	r9, [r4], #2
	cmp		r5, r3
	bne		innloop_pol_mul_schb
	sub		r4, r4, r3
	add		r4, r4, #2
	cmp		r4, r2
	bne		outloop_pol_mul_schb
	ldmia.w	sp!, {r4, r5, r6, r7, r8, r9, pc}
exit_pol_mul_schb:
	bx	lr
	nop
.endfunc

.global school_book_mul1 @this school book reads word by word
.func school_book_mul1, school_book_mul1
.type school_book_mul1, %function
school_book_mul1:
	cbz		r3, exit_pol_mul1_schb
	stmdb		sp!, {r4, r5, r6, r7, r8, r9, r10, r11, lr}
	lsl 		r3, r3, #1
	mov		r4, r2
	add		r2, r2, r3
outloop_pol_mul1_schb:
	mov		r5, #0
	ldr		r6, [r0], #4
innloop_pol_mul1_schb:
	ldr	r8, [r4], #4
	ldr	r7, [r1, r5]
	ldrh	r9, [r4], #-4
	add	r5, r5, #4
	
	lsr	r11,r8,#16 		@extract the top half-word
	@and	r10,r8,#0xFFFF		@extract the bottom half-word
	@mov	r10,r8

	smlabb	r8, r6, r7, r8
	smlabt	r11, r6, r7, r11
	smlatb	r11, r6, r7, r11
	smlatt	r9, r6, r7, r9

	lsl	r11, #16
	pkhbt	r8, r8, r11		@pack

	str	r8, [r4], #4
	strh	r9, [r4], #-4
	add	r4,r4,#4
	cmp		r5, r3
	bne		innloop_pol_mul1_schb
	sub		r4, r4, r3
	add		r4, r4, #4
	cmp		r4, r2
	bne		outloop_pol_mul1_schb
	ldmia.w	sp!, {r4, r5, r6, r7, r8, r9, r10, r11, pc}
exit_pol_mul1_schb:
	bx	lr
	nop
.endfunc

.global school_book_mul2 @this school book reads 2 words
.func school_book_mul2, school_book_mul2
.type school_book_mul2, %function
school_book_mul2:
	@cbz		r3, exit_pol_mul2_schb
	stmdb		sp!, {r4, r5, r6, r7, r8, r9, r10, r11, r12, lr}
	lsl 		r3, r3, #1
	mov		r4, r2
	add		r2, r2, r3
outloop_pol_mul2_schb:
	mov		r5, #0
	ldr		r6, [r0], #4 	@Load first two multiplicands
	ldr		r12, [r0], #4	@
innloop_pol_mul2_schb:

	ldr	r9, [r4], #0	@Load the result
	ldr	r7, [r1, r5]	@Load other multiplicands
	add	r5, r5, #4
	ldr	r8, [r1, r5]	@Load other multiplicands
	add	r5, r5, #4
	
	lsr	r11,r9,#16 		@extract the top half-word

	smlabb	r9, r7, r6, r9
	smladx	r11, r7, r6, r11

	pkhbt	r9, r9, r11, lsl #16	@pack
	str	r9, [r4], #4
@-------------------------------------------
	ldr	r9, [r4], #0
	lsr	r11,r9,#16 		@extract the top half-word

	smlabb	r9, r7, r12, r9
	smlatt	r9, r7, r6, r9
	smlabb	r9, r8, r6, r9

	smladx	r11, r7, r12, r11
	smladx	r11, r8, r6, r11

	pkhbt	r9, r9, r11, lsl #16		@pack
	str	r9, [r4], #4
@-------------------------------------------
	ldr	r9, [r4], #0
	lsr	r11,r9,#16 		@extract the top half-word
	
	smlatt r9, r7, r12, r9
	smlabb	r9, r8, r12, r9
	smlatt	r9, r8, r6, r9

	smladx	r11, r8, r12, r11

	pkhbt	r9, r9, r11, lsl #16	@pack
	str	r9, [r4], #4
@-------------------------------------------
	ldrh	r9, [r4], #0

	smlatt	r9, r8, r12, r9
	strh	r9, [r4], #-4

@-------------------------------------------

	cmp		r5, r3
	bne		innloop_pol_mul2_schb
	sub		r4, r4, r3
	add		r4, r4, #8
	cmp		r4, r2
	bne		outloop_pol_mul2_schb
	ldmia.w	sp!, {r4, r5, r6, r7, r8, r9, r10, r11, r12, pc}
exit_pol_mul2_schb:
	bx	lr
	nop
.endfunc

/*
.global school_book_mul2_16 @this school book reads 2 words. CAUTION : TAKES ONLY 16x16 polynomials.
.func school_book_mul2_16, school_book_mul2_16
.type school_book_mul2_16, %function
school_book_mul2_16:
	@cbz		r3, exit_pol_mul2_16_schb
	stmdb		sp!, {r3, r4, r5, r6, r7, r8, r9, r10, r11, r12, lr}
	@lsl 		r3, r3, #1
	mov		r4, r2
	add		r2, r2, #32
outloop_pol_mul2_16_schb:
	mov		r5, #0
	ldr		r6, [r0], #4 	@Load first two multiplicands
	ldr		r12, [r0], #4
	ldr		r3, [r0], #4
	ldr		r10, [r0], #4
	@ldrd		r6, r12, [r0], #8

innloop_pol_mul2_16_schb:

	ldr	r9, [r4], #0 	@Load the result
	@ldr	r7, [r1, r5]	@Load other multiplicands
	@add	r5, r5, #4
	@ldr	r8, [r1, r5]	@Load other multiplicands
	@add	r5, r5, #4
	ldr	r7, [r1], #4
	ldr	r8, [r1], #4

	lsr	r11,r9,#16 		@extract the top half-word

	smlabb	r9, r7, r6, r9
	smladx	r11, r7, r6, r11

	pkhbt	r9, r9, r11, lsl #16	@pack
	str	r9, [r4], #4
@-------------------------------------------
	ldr	r9, [r4], #0
	lsr	r11,r9,#16 		@extract the top half-word

	smlabb	r9, r7, r12, r9
	smlatt	r9, r7, r6, r9
	smlabb	r9, r8, r6, r9

	smladx	r11, r7, r12, r11
	smladx	r11, r8, r6, r11

	pkhbt	r9, r9, r11, lsl #16		@pack
	str	r9, [r4], #4
@-------------------------------------------
	ldr	r9, [r4], #0
	lsr	r11,r9,#16 		@extract the top half-word

	smlabb	r9, r3, r7, r9	
	smlatt r9, r7, r12, r9
	smlabb	r9, r8, r12, r9
	smlatt	r9, r8, r6, r9

	smladx	r11, r3, r7, r11
	smladx	r11, r8, r12, r11

	pkhbt	r9, r9, r11, lsl #16	@pack
	str	r9, [r4], #4
@-------------------------------------------
	ldr	r9, [r4], #0
	lsr	r11, r9, #16

	smlabb	r9, r10, r7, r9
	smlatt	r9, r3, r7, r9
	smlabb	r9, r3, r8, r9
	smlatt	r9, r8, r12, r9

	smladx	r11, r10, r7, r11
	smladx	r11, r3, r8, r11

	pkhbt	r9, r9, r11, lsl #16	@pack

	str	r9, [r4], #4
@-------------------------------------------
	ldr	r9, [r4], #0
	lsr	r11, r9, #16

	smlatt	r9, r10, r7, r9
	smlabb	r9, r10, r8, r9
	smlatt	r9, r3, r8, r9

	smladx	r11, r10, r8, r11

	pkhbt	r9, r9, r11, lsl #16	@pack

	str	r9, [r4], #4
@-------------------------------------------
	ldrh	r9, [r4], #0

	smlatt	r9, r10, r8, r9
	strh	r9, [r4], #-12
@-------------------------------------------
	add		r5, r5, #8
	cmp		r5, #32
	bne		innloop_pol_mul2_16_schb
	sub		r4, r4, #16
	sub		r1, r1, #32
	@sub		r4, r4, #32
	@add		r4, r4, #16
	cmp		r4, r2
	bne		outloop_pol_mul2_16_schb
	ldmia.w	sp!, {r3, r4, r5, r6, r7, r8, r9, r10, r11, r12, pc}
exit_pol_mul2_16_schb:
	bx	lr
	nop
.endfunc
*/


.global school_book_mul2_16 @this school book reads 2 words. CAUTION : TAKES ONLY 16x16 polynomials.
.func school_book_mul2_16, school_book_mul2_16
.type school_book_mul2_16, %function
school_book_mul2_16:
	@cbz		r3, exit_pol_mul2_16_schb
	stmdb		sp!, {r3, r4, r5, r6, r7, r8, r9, r10, r11, r12, lr}
	mov		r4, r2
	add		r2, r2, #32
outloop_pol_mul2_16_schb:
	mov		r5, #0
	ldr		r6, [r0], #4 	@Load first two multiplicands
	ldr		r12, [r0], #4
	ldr		r3, [r0], #4
	ldr		r10, [r0], #4

innloop_pol_mul2_16_schb:

	ldr	r9, [r4], #4 	@Load the result
	ldr	lr, [r4], #-4

	ldr	r7, [r1], #4
	ldr	r8, [r1], #4
	lsr	r11,r9,#16 		@extract the top half-word
	smlabb	r9, r7, r6, r9
	smladx	r11, r7, r6, r11
	pkhbt	r9, r9, r11, lsl #16	@pack
	str	r9, [r4], #4
@-------------------------------------------
	@ldr	r9, [r4], #0

	lsr	r11,lr,#16 		@extract the top half-word
	pkhbt	r9, r8, r7

	smladx	r11, r7, r12, r11
	smladx	r11, r8, r6, r11

	smlabb	lr, r7, r12, lr

	@smlatt	lr, r7, r6, lr
	@smlabb	lr, r8, r6, lr
	smlad	lr, r9, r6, lr

	pkhbt	lr, lr, r11, lsl #16		@pack
	str	lr, [r4], #4
@-------------------------------------------
	ldr	lr, [r4], #0
	lsr	r11,lr,#16 		@extract the top half-word

	smlabb	lr, r3, r7, lr	
	@smlatt lr, r7, r12, lr
	@smlabb	lr, r8, r12, lr
	smlad	lr, r9, r12, lr	

	smlatt	lr, r8, r6, lr

	smladx	r11, r3, r7, r11
	smladx	r11, r8, r12, r11

	pkhbt	lr, lr, r11, lsl #16	@pack
	str	lr, [r4], #4
@-------------------------------------------
	ldr	lr, [r4], #0
	lsr	r11, lr, #16

	smlabb	lr, r10, r7, lr
	@smlatt	lr, r3, r7, lr
	@smlabb	lr, r3, r8, lr
	smlad   lr, r3, r9, lr
	smlatt	lr, r8, r12, lr

	smladx	r11, r10, r7, r11
	smladx	r11, r3, r8, r11

	pkhbt	lr, lr, r11, lsl #16	@pack

	str	lr, [r4], #4
@-------------------------------------------
	ldr	lr, [r4], #0
	lsr	r11, lr, #16

	@smlatt	lr, r10, r7, lr
	@smlabb	lr, r10, r8, lr
	smlad	lr, r9, r10, lr
	smlatt	lr, r3, r8, lr

	smladx	r11, r10, r8, r11

	pkhbt	lr, lr, r11, lsl #16	@pack

	str	lr, [r4], #4
@-------------------------------------------
	ldrh	lr, [r4], #0

	smlatt	lr, r10, r8, lr
	strh	lr, [r4], #-12
@-------------------------------------------
	add		r5, r5, #8
	cmp		r5, #32
	bne		innloop_pol_mul2_16_schb
	sub		r4, r4, #16
	sub		r1, r1, #32
	@sub		r4, r4, #32
	@add		r4, r4, #16
	cmp		r4, r2
	bne		outloop_pol_mul2_16_schb
	ldmia.w	sp!, {r3, r4, r5, r6, r7, r8, r9, r10, r11, r12, pc}
exit_pol_mul2_16_schb:
	bx	lr
	nop
.endfunc

.global school_book_mul1_16 @this school book reads word by word
.func school_book_mul1_16, school_book_mul1_16
.type school_book_mul1_16, %function
school_book_mul1_16:

	stmdb		sp!, {r3,r4, r5, r6, r7, r8, r9, r10, r11, lr}

	mov		r3,#32
	mov		r4, r2
	add		r2, r2, r3
outloop_pol_mul1_16_schb:
	mov		r5, #0
	ldr		r6, [r0], #4
innloop_pol_mul1_16_schb:
	ldr	r8, [r4], #4
	ldr	r7, [r1, r5]
	ldrh	r9, [r4], #-4
	add	r5, r5, #4
	
	lsr	r11,r8,#16 		@extract the top half-word

	smlabb	r8, r6, r7, r8
	smladx	r11, r7, r6, r11
	smlatt	r9, r6, r7, r9

	lsl	r11, #16
	pkhbt	r8, r8, r11		@pack

	str	r8, [r4], #4
	strh	r9, [r4], #-4
	add	r4,r4,#4
	cmp		r5, r3
	bne		innloop_pol_mul1_16_schb
	sub		r4, r4, r3
	add		r4, r4, #4
	cmp		r4, r2
	bne		outloop_pol_mul1_16_schb
	ldmia.w	sp!, {r3, r4, r5, r6, r7, r8, r9, r10, r11, pc}
exit_pol_mul1_16_schb:
	bx	lr
	nop
.endfunc


//*************************************************
// unrolled_kara_mem(uint16_t* a, uint16_t* c, uint16_t* d, uint16_t k)
//*************************************************
.global unrolled_kara_mem
.func unrolled_kara_mem, unrolled_kara_mem
.type unrolled_kara_mem, %function
unrolled_kara_mem:
	stmdb	sp!, {r4, r5, r6, r7, r8, r9, r10, r11, r12, lr}
	mov		r9, r3, lsr #1
	mov		r11, r0
	mov		r12, r1
	add		r3, r3, r9
	ldr		r0, =kara_tmp
	add		r1, r2, r3, lsl #2
	mov		r4, #0
	sub		r1, r1, #2
// r11=a r12=c
	mov		r10, #0
//r0=kara_tmp r1=&d[3k-1] r2=d r3=k+k/2 r4=i r9=k/2
unrolled_kara_mem_loop1:
	ldr 	r5, [r2, r9, lsl #2]
	ldr 	r6, [r2]
	ldr 	r7, [r2, r3, lsl #1]
	ldr 	r8, [r2, r9, lsl #1]
	add		r4, r4, #2
	uadd16	r5, r5, r6
	uadd16	r7, r7, r8
	str		r10, [r2, r9, lsl #3]
	uadd16	r7, r7, r5
	str		r5, [r2, r9, lsl #2]
	str		r7, [r2, r3, lsl #1]
	ldr		r8, [r11, r3, lsl #1]
	ldr		r7, [r11, r9, lsl #1]
	ldr		r6, [r11, r9, lsl #2]
	ldr		r5, [r11], #4
	uadd16	r7, r7, r8
	uadd16	r5, r5, r6
	str		r7, [r1, r9, lsl #1]
	uadd16	r7, r5, r7
	str		r5, [r1], #4
	str		r7, [r0], #4
	cmp		r4, r9
	add		r2, r2, #4
	bne		unrolled_kara_mem_loop1
	sub		r1, r1, r9, lsl #2
	sub		r0, r0, r9, lsl #1
	add		r2, r2, r9, lsl #2
	sub		r11, r11, r9, lsl #1
	mov		r4, #0
//r0=kara_tmp r1=&d[2k+k/2-1] r2=&d[k+k/2] r3=k+k/2 r4=i
unrolled_kara_mem_loop1_2:
	ldr 	r8, [r12, r3, lsl #1]
	ldr 	r7, [r12, r9, lsl #2]
	ldr 	r6, [r12, r9, lsl #1]
	ldr		r5, [r12], #4
	uadd16 	r7, r7, r8
	uadd16 	r5, r5, r6
	uadd16	r5, r5, r7
	str 	r5, [r1, r4, lsl #1]
	add		r4, r4, #2
	cmp		r4, r9
	bne		unrolled_kara_mem_loop1_2
	mov		r3, r9
	sub		r12, r12, r9, lsl #1
//r0=kara_tmp r1=&d[2k+k/2-1] r2=&d[k+k/2] r3=k/2 r4=i
	#bl 		school_book_mul2
	bl 		school_book_mul2_16
	sub 	r0, r0, r4, lsl #1
	sub 	r2, r2, r4, lsl #1
	# mov	r3, r4
	mov		r4, #0
unrolled_kara_mem_loop2:
	ldr 	r8, [r2, r3, lsl #1]
	ldr 	r7, [r2]
	ldr 	r6, [r12, r3, lsl #2]
	ldr 	r5, [r12], #4
	add		r4, r4, #2
	uadd16	r7, r7, r8
	uadd16	r5, r5, r6
	str 	r10, [r2], #4
	str 	r7, [r1], #4
	str 	r5, [r0], #4
	cmp		r4, r3
	bne		unrolled_kara_mem_loop2
//undo last iteration for one of the segments
	usub16 	r7, r7, r8
	sub 	r2, r2, r3, lsl #2
	lsr 	r7, r7, #16
	sub 	r0, r0, r3, lsl #1
	strh 	r7, [r1, #-2]
//r0=kara_tmp r1=&d[3k-1] r2=&d[k] r3=k/2 r4=i r11=a r12=&c[k/2]
	# bl 		school_book_mul2
	bl 		school_book_mul2_16
	sub 	r0, r0, r4, lsl #1
	sub 	r1, r1, r4, lsl #1
	sub 	r2, r2, r4, lsl #1
	mov		r3, r9
	mov		r4, #0
//r0=kara_tmp r1=&d[2k+k/2-1] r2=&d[k] r3=k/2 r4=i r11=a r12=&c[k/2]
unrolled_kara_mem_loop3:
	ldr 	r8, [r2]
	ldr 	r7, [r1]
	ldr 	r6, [r2, r3, lsl #1]
	ldr 	r5, [r2, r3, lsl #2]
	usub16	r7, r7, r8
	ldr 	r8, [r12, r3, lsl #2]
	ldr 	r9, [r12], #4
	usub16	r5, r5, r6
	uadd16	r8, r8, r9
	add		r4, r4, #2
	str 	r5, [r2, r3, lsl #2]
	str 	r7, [r2, r3, lsl #1]
	str 	r10, [r1], #4
	str 	r8, [r0], #4
	cmp		r4, r3
	add		r2, r2, #4
	bne		unrolled_kara_mem_loop3
//undo last iteration for one of the segments
//prepare call to subroutine
	uadd16	r5, r5, r6
	add		r2, r2, r3, lsl #1
	lsr		r5, r5, #16
	sub 	r9, r3, #1
	add		r1, r1, r3, lsl #1
	sub 	r0, r0, r3, lsl #1
	strh	r5, [r2, r9, lsl #1]
//r0=kara_tmp r1=&d[3k+k/2-1] r2=&d[2k] r3=k/2 r4=i r11=a r12=&c[k]
	# bl 		school_book_mul2
	bl 		school_book_mul2_16
	add		r9, r4, r4, lsl #1
	sub 	r2, r2, r4, lsl #1
	sub 	r0, r0, r4, lsl #1
	sub 	r1, r1, r4, lsl #2
	sub 	r2, r2, r4, lsl #3
	# mov		r3, r4
	mov		r4, #0
	sub 	r12, r12, r3, lsl #2
//r0=kara_tmp r1=&d[2k+k/2-1] r2=d r3=k/2 r4=i r9=k+k/2 r11=a r12=c
unrolled_kara_mem_loop4:
	ldr 	r5, [r2]
	ldr 	r6, [r2, r3, lsl #1]
	ldr 	r7, [r2, r3, lsl #3]
	ldr 	r8, [r1, #2]
	uadd16	r5, r5, r6
	str 	r5, [r2, r3, lsl #1]
	usub16	r7, r7, r8
	ldr 	r5, [r2, r3, lsl #2]
	str 	r7, [r2, r3, lsl #3]
	str 	r10, [r2, r3, lsl #2]
	uadd16	r6, r5, r7
	str 	r6, [r1, r3, lsl #1]
	ldr 	r6, [r2, r9, lsl #1]
	add		r4, r4, #2
	add		r2, r2, #4
	usub16	r6, r6, r7
	cmp		r4, r3
	str 	r6, [r1, r3, lsl #2]
	add		r1, r1, #4
	bne 	unrolled_kara_mem_loop4
//undo last iteration for one of the arrays
	sub 	r1, r1, r9, lsl #1
	uadd16	r7, r7, r8
	usub16	r6, r6, r8
	uadd16	r8, r5, r7
	lsr		r5, #16
	strh 	r5, [r1]
	add		r1, r1, r3, lsl #2
	lsr		r7, #16
	strh	r7, [r1]
	add		r1, r1, r3, lsl #2
	lsr 	r8, #16
	strh 	r8, [r1, #-2]!
	add		r1, r1, r3, lsl #1
	lsr		r6, #16
	strh	r6, [r1], #2
	sub 	r1, r1, r9, lsl #2
	mov		r4, #0
	add		r1, r1, r3, lsl #1
//r0=kara_tmp r1=&d[k+k/2-1] r2=&d[k/2] r3=k/2 r4=i r9=k+k/2 r11=a r12=c
unrolled_kara_mem_loop4_2:
	ldr 	r5, [r11, r3, lsl #1]
	ldr 	r6, [r11], #4
	ldr 	r7, [r12, r3, lsl #1]
	ldr 	r8, [r12], #4
	uadd16	r5, r5, r6
	uadd16	r7, r7, r8
	str 	r5, [r1, r4, lsl #1]
	str 	r7, [r0, r4, lsl #1]
	add 	r4, r4, #2
	cmp 	r4, r3
	bne 	unrolled_kara_mem_loop4_2
	sub 	r11, r11, r3, lsl #1
	sub 	r12, r12, r3, lsl #1
	# bl 		school_book_mul2
	bl 		school_book_mul2_16
	sub 	r2, r2, r4, lsl #1
	# mov		r3, r4
	mov 	r4, #0
//r0=&kara_tmp[k/2] r1=&d[k+k/2-1] r2=&d[k/2] r3=k/2 r4=i r9=k+k/2 r11=a r12=c
unrolled_kara_mem_loop5:
	ldr 	r5, [r2, r3, lsl #1]
	ldr 	r6, [r2]
	str 	r10, [r2], #4
	add 	r4, r4, #2
	uadd16 	r5, r5, r6
	cmp 	r4, r3
	str 	r5, [r1], #4
	bne 	unrolled_kara_mem_loop5
	lsr 	r6, r6, #16
	sub 	r2, r2, r3, lsl #2
	strh 	r6, [r1, #-2]
	sub 	r9, r1, r3, lsl #1
	mov		r0, r11
	mov		r1, r12
//r0=a r1=c r2=d r3=k/2 r4=i r9=&d[k+k/2-1] r11=a r12=c
	# bl 		school_book_mul2
	bl 		school_book_mul2_16
	sub 	r2, r2, r4, lsl #1
	mov 	r3, #0
//r0=&a[k/2] r1=c r2=d r3=i r4=k/2 r9=&d[k+k/2-1] r11=a r12=c
unrolled_kara_mem_loop6:
	ldr 	r5, [r2, r4, lsl #1]
	ldr 	r6, [r2, r4, lsl #2]
	ldr 	r7, [r9]
	ldr 	r8, [r2]
	add 	r3, r3, #2
	usub16	r5, r6, r5
	usub16	r7, r7, r8
	str 	r5, [r2, r4, lsl #2]
	str 	r7, [r2, r4, lsl #1]
	cmp		r3, r4
	add 	r2, r2, #4
	str 	r10, [r9], #4
	bne 	unrolled_kara_mem_loop6
//r0=&a[k/2] r1=c r2=&d[k/2] r3=k/2 r4=k/2 r9=&d[2k-1] r11=a r12=c
//undo last iteration for one of the arrays
	sub 	r9, r9, r3, lsl #1
	add 	r2, r2, r3, lsl #1
	add 	r1, r1, r3, lsl #1
	strh 	r10, [r9]
//r0=&a[k/2] r1=&c[k/2] r2=&d[k] r3=k/2 r4=k/2 r9=&d[k+k/2-1] r11=a r12=c
	# bl 		school_book_mul2
	bl 		school_book_mul2_16
	add 	r1, r4, r4, lsl #1
	mov 	r3, #0
	sub 	r2, r2, r1, lsl #1
	add 	r9, r9, r4, lsl #2
// loop7 much optmzd!
//r0=&a[k] r1=k+k/2 r2=d r3=i r4=k/2 r9=&d[2k+k/2-1] r10=0 r12=c
unrolled_kara_mem_loop7:
	ldr 	r5, [r2, r4, lsl #1]
	ldr 	r6, [r2, r4, lsl #2]
	ldr 	r7, [r2, r4, lsl #3]
	ldr 	r11, [r2, r1, lsl #1]
	usub16	r5, r5, r6
	usub16 	r7, r7, r6
	ldr 	r8, [r9, r4, lsl #2]
	ldr 	r6, [r2]
	str 	r5, [r2, r4, lsl #1]
	usub16 	r8, r8, r5
	ldr 	r5, [r9, r4, lsl #1]
	str 	r8, [r2, r1, lsl #1]
	add 	r3, r3, #2
	usub16 	r5, r5, r6
	ldr 	r8, [r9, #2]
	str 	r5, [r2, r4, lsl #2]
	str 	r10, [r9, r4, lsl #1]
	uadd16	r5, r7, r11
	uadd16	r8, r7, r8
	str 	r5, [r2, r4, lsl #3]
	str 	r8, [r9, #2]!
	add 	r2, r2, #4
	cmp 	r4, r3
	add 	r9, r9, #2
	bne 	unrolled_kara_mem_loop7
	lsr 	r7, r7, #16
	add 	r2, r2, r4, lsl #3
	mov 	r3, #0
	strh 	r7, [r2, #-2]
	add 	r1, r9, r4, lsl #1
	mov		r11, r0
	add 	r12, r12, r4, lsl #2
	ldr 	r0, =kara_tmp
//r0=kara_tmp r1=&d[3k+k/2-1] r2=&d[2k+k/2] r3=i r4=k/2 r11=&a[k] r12=&c[k]
unrolled_kara_mem_loop7_2:
	ldr 	r5, [r11, r4, lsl #1]
	ldr 	r6, [r11], #4
	ldr 	r7, [r12, r4, lsl #1]
	ldr 	r8, [r12], #4
	add 	r3, r3, #4
	uadd16	r5, r5, r6
	uadd16	r7, r7, r8
	ldr 	r6, [r11, r4, lsl #1]
	ldr 	r8, [r11], #4
	ldr 	r9, [r12, r4, lsl #1]
	ldr 	r10, [r12], #4
	str 	r5, [r1], #4
	str 	r7, [r0], #4
	uadd16	r6, r6, r8
	uadd16	r9, r9, r10
	cmp 	r3, r4
	str 	r6, [r1], #4
	str 	r9, [r0], #4
	bne 	unrolled_kara_mem_loop7_2
	sub 	r0, r0, r4, lsl #1
	sub 	r1, r1, r4, lsl #1
	mov 	r10, #0
//r0=kara_tmp r1=&d[3k+k/2-1] r2=&d[2k+k/2] r3=k/2 r4=k/2 r10=0 r11=&a[k+k/2] r12=&c[k+k/2]
	bl 		school_book_mul2_16
	sub 	r2, r2, r4, lsl #1
	sub 	r0, r11, r4, lsl #1
	mov		r9, r1
	mov 	r3, #0
	sub 	r1, r12, r4, lsl #1
//r0=&a[k] r1=&c[k] r2=&d[2k+k/2] r3=i r4=k/2 r9=&d[3k+k/2-1] r10=0 r11=&a[k+k/2] r12=&c[k+k/2]
unrolled_kara_mem_loop8:
	ldr 	r5, [r2, r4, lsl #1]
	ldr 	r6, [r2], #4
	ldr 	r7, [r2, r4, lsl #1]
	ldr 	r8, [r2], #4
	add 	r3, r3, #4
	uadd16	r5, r5, r6
	uadd16	r7, r7, r8
	str 	r5, [r9], #4
	str 	r10, [r2, #-8]
	cmp		r3, r4
	str 	r7, [r9], #4
	str 	r10, [r2, #-4]
	bne		unrolled_kara_mem_loop8
	lsr		r8, r8, #16
	sub 	r2, r2, r4, lsl #2
	strh	r8, [r9, #-2]
//r0=&a[k] r1=&c[k] r2=&d[2k] r3=i r4=k/2 r9=&d[4k-1] r10=0 r11=&a[k+k/2] r12=&c[k+k/2]
	bl 		school_book_mul2_16
	sub 	r9, r9, r4, lsl #1
	sub 	r2, r2, r4, lsl #1
	mov		r3, #0
	mov		r1, r12
//r0=&a[k+k/2] r1=&c[k+k/2] r2=&d[2k] r3=i r4=k/2 r9=&d[3k+k/2-1] r10=0
unrolled_kara_mem_loop9:
	ldr 	r5, [r2, r4, lsl #1]
	ldr 	r6, [r2, r4, lsl #2]
	ldr 	r7, [r2], #4
	ldr 	r8, [r9], #4
	ldr 	r11, [r2, r4, lsl #1]
	ldr 	r12, [r2, r4, lsl #2]
	usub16	r5, r6, r5
	usub16	r7, r8, r7
	ldr 	r6, [r2], #-4
	ldr 	r8, [r9], #4
	str 	r5, [r2, r4, lsl #2]
	str 	r7, [r2, r4, lsl #1]
	usub16	r5, r12, r11
	usub16	r7, r8, r6
	add		r2, r2, #4
	add		r3, r3, #4
	str 	r5, [r2, r4, lsl #2]
	str 	r7, [r2, r4, lsl #1]
	add		r2, r2, #4
	cmp 	r3, r4
	str 	r10, [r9, #-8]
	str 	r10, [r9, #-4]
	bne 	unrolled_kara_mem_loop9
	mov		r12, r2
	add 	r2, r2, r4, lsl #1
	sub 	r11, r2, r4, lsl #3
//r0=&a[k+k/2] r1=&c[k+k/2] r2=&d[3k] r3=k/2 r4=k/2 r9=&d[4k-1] r10=0 r11=d[k] r12=d[2k+k/2]
	bl 		school_book_mul2_16
//r0=&a[2k] r1=&c[k+k/2] r2=&d[3k+k/2] r3=k/2 r4=k/2 r9=&d[4k-1] r10=0 r11=d[k] r12=d[2k+k/2]
unrolled_kara_mem_loop10:
	ldr 	r7, [r12]
	ldr 	r8, [r12, r3, lsl #1]
	ldr 	r5, [r11, r3, lsl #1]
	ldr 	r4, [r11]
	usub16 	r7, r7, r8
	ldr 	r6, [r11, r3, lsl #2]
	usub16	r5, r5, r7
	ldr 	r9, [r12, r3, lsl #2]
	str 	r5, [r11, r3, lsl #1]
	usub16	r4, r4, r6
	usub16	r5, r8, r9
	add		r10, r10, #2
	usub16	r6, r6, r5
	usub16 	r9, r7, r9
	str 	r5, [r12, r3, lsl #1]
	str 	r6, [r11, r3, lsl #2]
	cmp 	r3, r10
	str 	r4, [r11], #4
	str 	r9, [r12], #4
	bne 	unrolled_kara_mem_loop10
//r0=&a[2k] r1=&c[k+k/2] r2=&d[3k+k/2] r3=k/2 r4=k/2 r9=&d[4k-1] r10=0 r11=d[k+k/2] r12=d[3k]
	lsr 	r7, r7, #16
	lsr 	r8, r8, #16
	strh 	r7, [r12, #-2]
	strh	r8, [r2, #-2]
	ldmia.w	sp!, {r4, r5, r6, r7, r8, r9, r10, r11, r12, pc}
exit_unrolled_kara_mem:
	bx	lr
	nop
.endfunc

//*************************************************
// bigloop( a_1, b_1, result_final, result_d01, d01, d23, d0123 )
//*************************************************
.global karatsuba_simple_asm
.func karatsuba_simple_asm, karatsuba_simple_asm
.type karatsuba_simple_asm, %function
karatsuba_simple_asm:
	stmdb	sp!, {r4, r5, r6, r7, r8, r9, r10, r11, r12, lr}
//d01=d[0] d0123=d[N/2-1] d23=d[N-2] result_d01=d[N+N/2-3]
// N/2-1 = 31(62)   N-2 = 62(124)  N+N/2-3 = 93(186)
//r0=a_1 r1=b_1 r2=result_final r3=d r14=lr=...
//N=64(128), N/4=16(32) 2*N/4=32(64) 3*N/4=48(96) 4*N/4=64(128) 6*N/4=96(192)
	add 	r14, r0, #32
	push 	{r14}
outloop_karatsuba_simple_asm:
	ldr 	r7, [r0, #96]
	ldr 	r6, [r0, #64]
	ldr 	r5, [r0, #32]
	ldr 	r4, [r0], #4
	push 	{r0}
	add 	r0, r1, #32
	push 	{r0}
inloop_karatsuba_simple_asm:
	ldr 	r11, [r1, #96]
	ldr 	r10, [r1, #64]
	ldr 	r9, [r1, #32]
	ldr 	r8, [r1], #4
	ldr 	r12, [r2, #192]
	ldr 	r0, [r2, #196]
	lsr 	r14, r12, #16
	smlabb 	r12, r7, r11, r12
	smladx 	r14, r7, r11, r14
	smlatt 	r0, r7, r11, r0
	pkhbt	r12, r12, r14, lsl #16	@pack
	strh 	r0, [r2, #196]
	str 	r12, [r2, #192]
	ldr 	r12, [r2, #128]
	ldr 	r0, [r2, #132]
	lsr 	r14, r12, #16
	smlabb 	r12, r6, r10, r12
	smladx 	r14, r6, r10, r14
	smlatt 	r0, r6, r10, r0
	pkhbt	r12, r12, r14, lsl #16	@pack
	strh 	r0, [r2, #132]
	str 	r12, [r2, #128]
	ldr 	r12, [r2, #64]
	ldr 	r0, [r2, #68]
	lsr 	r14, r12, #16
	smlabb 	r12, r5, r9, r12
	smladx 	r14, r5, r9, r14
	smlatt 	r0, r5, r9, r0
	pkhbt	r12, r12, r14, lsl #16	@pack
	strh 	r0, [r2, #68]
	str 	r12, [r2, #64]
	ldr 	r12, [r2]
	ldr 	r0, [r2, #4]
	lsr 	r14, r12, #16
	smlabb 	r12, r4, r8, r12
	smladx 	r14, r4, r8, r14
	smlatt 	r0, r4, r8, r0
	pkhbt	r12, r12, r14, lsl #16	@pack
	strh 	r0, [r2, #4]
	str 	r12, [r2], #4
//N+N/2-3 = 93(186) 2*N/4=32(64) ++=125(250)
	ldr 	r14, [r3, #186]
	uadd16 	r10, r8, r10
	uadd16 	r12, r4, r6
	lsr 	r0, r14, #16
	smlabb 	r14, r10, r12, r14
	smladx 	r0, r10, r12, r0
	pkhbt 	r14, r14, r0, lsl #16
	str 	r14, [r3, #186]
	ldr 	r0, [r3, #190]
	ldr 	r14, [r3, #250]
	smlatt 	r0, r12, r10, r0
	strh	r0, [r3, #190]
	lsr		r0, r14, #16
	uadd16 	r11, r11, r9
	uadd16	r12, r5, r7
	smlabb 	r14, r11, r12, r14
	smladx 	r0, r11, r12, r0
	pkhbt 	r14, r14, r0, lsl #16
	str 	r14, [r3, #250]
	ldr 	r0, [r3, #254]
	usub16 	r11, r11, r9
	usub16	r10, r10, r8
	smlatt 	r0, r11, r12, r0
	strh 	r0, [r3, #254]
	uadd16 	r9, r8, r9
	uadd16	r12, r4, r5
//free acc5(r8)!
	ldr 	r14, [r3]
	ldr 	r8, [r3, #4]
	lsr 	r0, r14, #16
	smlabb 	r14, r12, r9, r14
	smladx 	r0, r12, r9, r0
	smlatt 	r8, r12, r9, r8
	pkhbt 	r14, r14, r0, lsl #16
	strh	r8, [r3, #4]
	str 	r14, [r3], #4
//beware r3 already incremented
// N-2 = 62(124) N/2-1 = 31(62)
	ldr 	r14, [r3, #120]
	ldr 	r8, [r3, #124]
	uadd16 	r12, r6, r7
	uadd16 	r11, r10, r11
	lsr 	r0, r14, #16
//free acc7(r10)!
	smlabb 	r14, r11, r12, r14
	smladx	r0, r11, r12, r0
	smlatt 	r8, r11, r12, r8
	pkhbt 	r14, r14, r0, lsl #16
	str 	r14, [r3, #120]
	strh	r8, [r3, #124]
	ldr 	r14, [r3, #58]
	ldr 	r8, [r3, #62]
	uadd16	r12, r12, r5
	uadd16 	r11, r11, r9
	uadd16 	r12, r12, r4
	lsr 	r0, r14, #16
//free acc6(r9)!
	ldr 	r9, [sp]
	smlabb 	r14, r11, r12, r14
	smladx	r0, r11, r12, r0
	smlatt 	r8, r11, r12, r8
	pkhbt 	r14, r14, r0, lsl #16
	cmp 	r9, r1
	strh 	r8, [r3, #62]
	str 	r14, [r3, #58]
	bne 	inloop_karatsuba_simple_asm
	pop 	{r9}
	pop 	{r0}
	pop 	{r14}
	sub 	r1, r1, #32
	sub 	r2, r2, #28
	sub 	r3, r3, #28
	cmp 	r0, r14
	push 	{r14}
	bne 	outloop_karatsuba_simple_asm
	pop 	{r14}
//r0=a_1[N/4] r1=b_1 r2=result_final[N/4] r3=d[N/4] r14=lr=...
	ldmia.w	sp!, {r4, r5, r6, r7, r8, r9, r10, r11, r12, pc}
exit_karatsuba_simple_asm:
	bx	lr
	nop
.endfunc


//function for a TOOM-COOK in assembly
//*************************************************
// toom_cook_4way_mem_asm(&a1, &b1, &result);
//*************************************************
.global toom_cook_4way_mem_asm
.func toom_cook_4way_mem_asm, toom_cook_4way_mem_asm
.type toom_cook_4way_mem_asm, %function
toom_cook_4way_mem_asm:
	stmdb	sp!, {r3, r4, r5, r6, r7, r8, r9, r10, r11, r12, lr}
	subw	sp, sp, #3076
	mov 	r3, sp
	mov 	r11, r2
	add 	r2, r3, #1792
//r0=&a1[i] r1=b1 r2=&w_m[896+i]
//loop1
	add 	r12, r0, #128
toom_cook_4way_mem_asm_loop1:
	ldr 	r6, [r0, #384]
	ldr 	r4, [r0, #128]
	ldr 	r5, [r0, #256]
	ldr 	r3, [r0], #4
	cmp 	r0, r12
	uadd16 	r8, r4, r6
	uadd16 	r7, r3, r5
	uadd16 	r9, r7, r8
	usub16 	r10, r7, r8
	str 	r9, [r2, #128]
	str 	r10, [r2, #256]
	uadd16 	r7, r3, r3
	uadd16	r7, r7, r7
	uadd16	r7, r7, r5
	uadd16 	r7, r7, r7
	uadd16 	r8, r4, r4
	uadd16 	r8, r8, r8
	uadd16 	r8, r8, r6
	uadd16 	r9, r7, r8
	usub16 	r10, r7, r8
	str 	r9, [r2, #384]
	str 	r10, [r2, #512]
	uadd16 	r7, r6, r6
	uadd16 	r7, r7, r5
	uadd16 	r7, r7, r7
	uadd16 	r7, r7, r4
	uadd16 	r7, r7, r7
	uadd16 	r7, r7 ,r3
	str 	r7, [r2], #4
	bne 	toom_cook_4way_mem_asm_loop1
//loop2
	add 	r12, r1, #128
	add 	r2, r2, #512
toom_cook_4way_mem_asm_loop2:
	ldr 	r6, [r1, #384]
	ldr 	r4, [r1, #128]
	ldr 	r5, [r1, #256]
	ldr 	r3, [r1], #4
	cmp 	r1, r12
	uadd16 	r8, r4, r6
	uadd16 	r7, r3, r5
	uadd16 	r9, r7, r8
	usub16 	r10, r7, r8
	str 	r9, [r2, #128]
	str 	r10, [r2, #256]
	uadd16 	r7, r3, r3
	uadd16	r7, r7, r7
	uadd16	r7, r7, r5
	uadd16 	r7, r7, r7
	uadd16 	r8, r4, r4
	uadd16 	r8, r8, r8
	uadd16 	r8, r8, r6
	uadd16 	r9, r7, r8
	usub16 	r10, r7, r8
	str 	r9, [r2, #384]
	str 	r10, [r2, #512]
	uadd16 	r7, r6, r6
	uadd16 	r7, r7, r5
	uadd16 	r7, r7, r7
	uadd16 	r7, r7, r4
	uadd16 	r7, r7, r7
	uadd16 	r7, r7 ,r3
	str 	r7, [r2], #4
	bne 	toom_cook_4way_mem_asm_loop2
//r0=&a1[i+64] r1=&b1[i+64] r2=&w_m[896+320+64+i] r3-r10 r12
	mov 	r9, r11
	mov 	r11, r0
	mov 	r4, r1
	sub 	r0, r2, #2560
	mov 	r1, #0
	mov 	r2, #1792
	bl 		memset
	mov 	r10, r0
	sub 	r0, r11, #128
	sub 	r1, r4, #128
	add		r2, r10, #1536
	mov 	r3, #32
	bl 		unrolled_kara_mem
	add 	r0, r10, #2304
	add 	r1, r10, #2944
	add 	r2, r10, #1280
	mov 	r3, #32
	bl 		unrolled_kara_mem
	add 	r0, r10, #2176
	add 	r1, r10, #2816
	add 	r2, r10, #1024
	mov 	r3, #32
	bl 		unrolled_kara_mem
	add 	r0, r10, #2048
	add 	r1, r10, #2688
	add 	r2, r10, #768
	mov 	r3, #32
	bl 		unrolled_kara_mem
	add 	r0, r10, #1920
	add 	r1, r10, #2560
	add 	r2, r10, #512
	mov 	r3, #32
	bl 		unrolled_kara_mem
	add 	r0, r10, #1792
	add 	r1, r10, #2432
	add 	r2, r10, #256
	mov 	r3, #32
	bl 		unrolled_kara_mem
	add 	r0, r11, #256
	add 	r1, r4, #256
	mov 	r2, r10
	mov 	r3, #32
	bl 		unrolled_kara_mem
//r0=i r9=result r10=w_m
//loop3
	movw	r14, #0x1fff
	mov 	r0, #0
	mov 	r11, #128
	add 	r14, r14, r14, lsl #16
toom_cook_4way_mem_asm_loop3:
	ldr 	r5, [r10, #1024]
	ldr 	r3, [r10, #512]
	ldr 	r2, [r10, #256]
	ldr 	r7, [r10, #1536]
	ldr 	r6, [r10, #1280]
	ldr 	r4, [r10, #768]
	ldr 	r1, [r10], #4
	uadd16	r2, r2, r5
	usub16	r6, r6, r5
	usub16	r4, r4, r3
	usub16	r5, r5, r1
	lsr 	r11, r7, #16
	lsl		r12, r7, #6
	pkhbt	r8, r12, r11, lsl #22
	usub16	r5, r5, r8
	lsl		r11, r4, #16
	lsr	 	r12, r4, #1
	pkhtb 	r4, r12, r11, asr #17
	uadd16	r3, r3, r4
	lsr 	r11, r5, #16
	lsl		r12, r5, #1
	pkhbt	r8, r12, r11, lsl #17
	uadd16	r5, r6, r8
	lsr 	r11, r3, #16
	lsl		r12, r3, #6
	pkhbt	r8, r12, r11, lsl #22
	uadd16	r8, r3, r8
	usub16	r2, r2, r8
	usub16	r3, r3, r7
	usub16 	r3, r3, r1
	mov 	r8, #45
	smulbb 	r11, r3, r8
	smultb 	r12, r3, r8
	pkhbt 	r8, r11, r12, lsl #16
	uadd16	r2, r2, r8
	lsr 	r11, r3, #16
	lsl		r12, r3, #3
	pkhbt	r8, r12, r11, lsl #19
	usub16	r5, r5, r8
	movw 	r8, #43691
	smulbb 	r11, r5, r8
	smultb 	r12, r5, r8
	pkhbt 	r8, r11, r12, lsl #16
	lsl		r11, r8, #16
	lsr	 	r12, r8, #3
	pkhtb 	r5, r12, r11, asr #19
	uadd16 	r6, r2, r6
	lsr 	r11, r4, #16
	lsl		r12, r4, #4
	pkhbt	r8, r12, r11, lsl #20
	uadd16	r2, r2, r8
	movw 	r8, #36409
	smulbb 	r11, r2, r8
	smultb 	r12, r2, r8
	pkhbt 	r8, r11, r12, lsl #16
	lsl		r11, r8, #16
	lsr	 	r12, r8, #1
	pkhtb 	r2, r12, r11, asr #17
	usub16	r3, r3, r5
	uadd16	r4, r4, r2
	mov		r8, #0
	usub16	r4, r8, r4
	mov 	r8, #30
	smulbb 	r11, r2, r8
	smultb 	r12, r2, r8
	pkhbt 	r8, r11, r12, lsl #16
	usub16	r6, r8, r6
	movw 	r8, #61167
	smulbb 	r11, r6, r8
	smultb 	r12, r6, r8
	pkhbt 	r8, r11, r12, lsl #16
	lsl		r11, r8, #16
	lsr	 	r12, r8, #1
	pkhtb 	r6, r12, r11, asr #17
//r0=i r9=result r10=w_m
	usub16	r7, r7, r3
	usub16	r6, r6, r2
	usub16	r5, r5, r1
	ldr 	r3, [r9]
	ldr 	r2, [r9, #128]
	ldr 	r1, [r9, #256]
	uadd16	r7, r7, r3
	uadd16	r6, r6, r2
	uadd16	r5, r5, r1
	#ubfx	r7, r7, #0, #13
	#ubfx	r6, r6, #0, #13
	#ubfx	r5, r5, #0, #13
	and		r7, r7, r14
	and		r6, r6, r14
	and 	r5, r5, r14
	cmp 	r0, #63
	add		r0, r0, #2
	sub 	r1, r9, #128
	it 		ls
	addls	r1, r1, #512
	ldr 	r8, [r1]
	str 	r5, [r9, #256]
	str 	r6, [r9, #128]
	str 	r7, [r9], #4
	usub16	r3, r8, r4
	it 		ls
	uadd16ls r3, r8, r4
	cmp 	r0, #126
	#ubfx 	r3, r3, #0, #13
	and 	r3, r3, r14
	str 	r3, [r1]
	bne		toom_cook_4way_mem_asm_loop3
//for(i=0;i<2*(SABER_N/4)-1;i++){ last iteration only halfword!!!
	ldrh 	r5, [r10, #1024]
	ldrh 	r3, [r10, #512]
	ldrh 	r2, [r10, #256]
	ldrh 	r7, [r10, #1536]
	ldrh 	r6, [r10, #1280]
	ldrh 	r4, [r10, #768]
	ldrh 	r1, [r10], #4
	uadd16	r2, r2, r5
	usub16	r6, r6, r5
	usub16	r4, r4, r3
	usub16	r5, r5, r1
	lsl 	r8, r7, #6
	usub16	r5, r5, r8
	lsr		r4, r4, #1
	uadd16	r3, r3, r4
	lsl		r8, r5, #1
	uadd16	r5, r6, r8
	lsl		r8, r3, #6
	uadd16	r8, r3, r8
	usub16	r2, r2, r8
	usub16	r3, r3, r7
	usub16 	r3, r3, r1
	mov 	r8, #45
	smulbb 	r8, r3, r8
	uadd16	r2, r2, r8
	lsl		r8, r3, #3
	usub16	r5, r5, r8
	movw 	r8, #43691
	smulbb 	r8, r5, r8
	lsr		r5, r8, #3
	uadd16 	r6, r2, r6
	lsl		r8, r4, #4
	uadd16	r2, r2, r8
	movw 	r8, #36409
	smulbb 	r8, r2, r8
	lsr 	r2, r8, #1
	usub16	r3, r3, r5
	uadd16	r4, r4, r2
	mov		r8, #0
	usub16	r4, r8, r4
	mov 	r8, #30
	smulbb 	r8, r2, r8
	usub16	r6, r8, r6
	movw 	r8, #61167
	smulbb 	r8, r6, r8
	lsr		r6, r8, #1
//r0=i r9=result r10=w_m
	usub16	r7, r7, r3
	usub16	r6, r6, r2
	usub16	r5, r5, r1
	ldrh 	r3, [r9]
	ldrh 	r2, [r9, #128]
	ldrh 	r1, [r9, #256]
	uadd16	r7, r7, r3
	uadd16	r6, r6, r2
	uadd16	r5, r5, r1
	ubfx	r7, r7, #0, #13
	ubfx	r6, r6, #0, #13
	ubfx	r5, r5, #0, #13
	cmp 	r0, #63
	add		r0, r0, #2
	sub 	r1, r9, #128
	it 		ls
	addls	r1, r1, #512
	ldr 	r8, [r1]
	strh 	r5, [r9, #256]
	strh 	r6, [r9, #128]
	strh 	r7, [r9], #4
	usub16	r3, r8, r4
	it 		ls
	uadd16ls r3, r8, r4
	cmp 	r0, #126
	ubfx 	r3, r3, #0, #13
	strh 	r3, [r1]
exit_toom_cook_4way_mem_asm:
	addw	sp, sp, #3076
	ldmia.w	sp!, {r3, r4, r5, r6, r7, r8, r9, r10, r11, r12, pc}
	bx	lr
	nop
.endfunc


//function for classical Karatsuba in assembly
//*************************************************
// karatsuba_asm(&a, &b, &result);
//*************************************************
.global karatsuba_asm
.func karatsuba_asm, karatsuba_asm
.type karatsuba_asm, %function
karatsuba_asm:
	stmdb	sp!, {r3, r4, r5, r6, r7, r8, r9, r10, r11, r12, lr}
	#mov 	r9, r3
	subw	sp, sp, #940
	mov 	r9, sp
	mov 	r10, r0
	mov 	r11, r1
	mov 	r12, r2
	# mov 	r0, r9
	add 	r0, r9, #378
	mov 	r1, #0
	# movw	r2, #936
	movw 	r2, #558
	push {r12}
	bl 		memset
	pop {r12}
//r10=&a[i] r11=&b[i] r12=&result[i] r9=&result_m[i] r0=&b[16]
	add 	r0, r11, #32
//loop1
karatsuba_asm_loop1:
	ldr 	r1, [r11, #96]
	ldr 	r2, [r11, #64]
	ldr 	r4, [r11, #32]
	ldr 	r5, [r11], #4
	uadd16 	r6, r1, r2
	uadd16 	r7, r1, r4
	uadd16 	r8, r2, r5
	uadd16 	r3, r4, r5
	uadd16 	r5, r7, r8
	str 	r8, [r9, #192]
	str 	r7, [r9, #224]
	str 	r6, [r9, #256]
	str 	r5, [r9, #288]
	str 	r3, [r9, #160]
	cmp 	r11, r0
	ldr 	r1, [r10, #96]
	ldr 	r2, [r10, #64]
	ldr 	r4, [r10, #32]
	ldr 	r5, [r10], #4
	uadd16 	r6, r1, r2
	uadd16 	r7, r1, r4
	uadd16 	r8, r2, r5
	uadd16 	r3, r4, r5
	uadd16 	r5, r7, r8
	str 	r8, [r9, #32]
	str 	r7, [r9, #64]
	str 	r6, [r9, #96]
	str 	r5, [r9, #128]
	str 	r3, [r9], #4
	bne 	karatsuba_asm_loop1
//r10=&a[i+16] r11=&b[i+16] r12=&result[i] r9=&result_m[i+16] r0=&b[16]
	sub 	r0, r10, #32
	sub 	r1, r11, #32
	add 	r2, r9, #346
	bl 		school_book_mul2_16
	add 	r1, r1, #32
	add 	r2, r2, #30
	bl 		school_book_mul2_16
	add 	r1, r1, #32
	add 	r2, r2, #30
	bl 		school_book_mul2_16
	add 	r1, r1, #32
	add 	r2, r2, #30
	bl 		school_book_mul2_16
	sub 	r0, r9, #32
	add 	r1, r9, #128
	add 	r2, r2, #30
	bl 		school_book_mul2_16
	add 	r1, r1, #32
	add 	r2, r2, #30
	bl 		school_book_mul2_16
	add 	r1, r1, #32
	add 	r2, r2, #30
	bl 		school_book_mul2_16
	add 	r1, r1, #32
	add 	r2, r2, #30
	bl 		school_book_mul2_16
	add 	r1, r1, #32
	add 	r2, r2, #30
	bl 		school_book_mul2_16
//r10=&a[i+16] r11=&b[i+16] r12=&result[i] r9=&result_m[i+16]
	sub 	r0, r9, #32
	mov 	r1, #0
	mov 	r2, #378
	mov 	r11, r0
	push {r12}
	bl 		memset
	pop {r12}
	add 	r10, r11, #60
karatsuba_asm_inter1:
	ldr 	r0, [r11, #502]
	ldr 	r1, [r11, #506]
	ldr 	r2, [r11, #510]
	ldr 	r3, [r11, #688]
	ldr 	r4, [r11, #692]
	ldr 	r5, [r11, #696]
	ldr 	r6, [r11, #378]
	ldr 	r7, [r11, #382]
	ldr 	r8, [r11, #386]
	str 	r0, [r11, #252]
	str 	r1, [r11, #256]
	str 	r2, [r11, #260]
	str 	r3, [r11, #126]
	str 	r4, [r11, #130]
	str 	r5, [r11, #134]
	str 	r8, [r11, #8]
	str 	r7, [r11, #4]
	str 	r6, [r11], #12
	cmp 	r10, r11
	bne 	karatsuba_asm_inter1
	ldrh 	r0, [r11, #502]
	ldrh 	r3, [r11, #688]
	ldrh 	r6, [r11, #378]
	strh 	r0, [r11, #252]
	strh 	r3, [r11, #126]
	strh 	r6, [r11], #-28
//r11=&result_m[i+16] r12=&result[i] r10=&result_m[i+46]
	add 	r10, r10, #32
karatsuba_asm_inter2:
	ldr 	r0, [r11, #252]
	ldr 	r1, [r11, #780]
	ldr 	r2, [r11, #470]
	ldr 	r3, [r11, #532]
	ldr 	r4, [r11, #126]
	ldr 	r5, [r11, #842]
	ldr 	r6, [r11, #656]
	ldr 	r7, [r11, #718]
	uadd16 	r8, r0, r1
	uadd16 	r9, r2, r3
	ldr 	r0, [r11]
	ldr 	r1, [r11, #594]
	ldr 	r2, [r11, #346]
	ldr 	r3, [r11, #408]
	usub16	r8, r8, r9
	uadd16 	r4, r4, r5
	uadd16 	r6, r6, r7
	usub16 	r4, r4, r6
	uadd16 	r0, r0, r1
	uadd16 	r2, r2, r3
	usub16 	r0, r0, r2
	str 	r8, [r11, #252]
	str 	r4, [r11, #126]
	str 	r0, [r11], #4
	cmp 	r10, r11
	bne 	karatsuba_asm_inter2
	ldrh 	r0, [r11, #252]
	ldrh 	r1, [r11, #780]
	ldrh 	r2, [r11, #470]
	ldrh 	r3, [r11, #532]
	ldrh 	r4, [r11, #126]
	ldrh 	r5, [r11, #842]
	ldrh 	r6, [r11, #656]
	ldrh 	r7, [r11, #718]
	uadd16 	r8, r0, r1
	uadd16 	r9, r2, r3
	ldrh 	r0, [r11]
	ldrh 	r1, [r11, #594]
	ldrh 	r2, [r11, #346]
	ldrh 	r3, [r11, #408]
	usub16	r8, r8, r9
	uadd16 	r4, r4, r5
	uadd16 	r6, r6, r7
	usub16 	r4, r4, r6
	uadd16 	r0, r0, r1
	uadd16 	r2, r2, r3
	usub16 	r0, r0, r2
	strh 	r8, [r11, #252]
	strh 	r4, [r11, #126]
	strh 	r0, [r11], #-28
//r11=&result_m[i+32] r12=&result[i] r10=&result_m[i+62]
	add 	r10, r10, #32
karatsuba_asm_inter3:
	ldr 	r0, [r11, #500]
	ldr 	r1, [r11, #252]
	ldr 	r2, [r11, #686]
	ldr 	r3, [r11, #126]
	ldr 	r4, [r11, #376]
	ldr 	r5, [r11]
	uadd16 	r6, r0, r1
	uadd16 	r7, r2, r3
	uadd16 	r8, r4, r5
	str 	r6, [r11, #252]
	str 	r7, [r11, #126]
	str 	r8, [r11], #4
	cmp 	r10, r11
	bne 	karatsuba_asm_inter3
	ldrh 	r0, [r11, #500]
	ldrh 	r1, [r11, #252]
	ldrh 	r2, [r11, #686]
	ldrh 	r3, [r11, #126]
	ldrh 	r4, [r11, #376]
	ldrh 	r5, [r11]
	uadd16 	r6, r0, r1
	uadd16 	r7, r2, r3
	uadd16 	r8, r4, r5
	strh 	r6, [r11, #252]
	strh 	r7, [r11, #126]
	strh 	r8, [r11], #-124
//r11=&result_m[i] r12=&result[i] r10=&result_m[i+62]
//TODO: i think memcpy is not much optimized and it can be improved
	mov 	r0, r12
	mov 	r1, r11
	mov 	r2, #126
	bl 		memcpy
	add 	r12, r12, #64
//r11=&result_m[i] r12=&result[i+32] r10=&result_m[i+62]
karatsuba_asm_final:
	ldr 	r0, [r12, #64]
	ldr 	r1, [r11, #252]
	ldr 	r2, [r12]
	ldr 	r3, [r11, #126]
	ldr 	r4, [r11], #4
	uadd16 	r5, r0, r1
	uadd16 	r6, r2, r3
	uadd16 	r7, r4, r1
	cmp 	r10, r11
	usub16 	r8, r6, r7
	str 	r5, [r12, #64]
	str 	r8, [r12], #4
	bne 	karatsuba_asm_final
	ldrh 	r0, [r12, #64]
	ldrh 	r1, [r11, #252]
	ldrh 	r2, [r12]
	ldrh 	r3, [r11, #126]
	ldrh 	r4, [r11], #4
	uadd16 	r5, r0, r1
	uadd16 	r6, r2, r3
	uadd16 	r7, r4, r1
	usub16 	r8, r6, r7
	strh 	r5, [r12, #64]
	strh 	r8, [r12], #4
exit_karatsuba_asm:
	addw	sp, sp, #940
	ldmia.w	sp!, {r3, r4, r5, r6, r7, r8, r9, r10, r11, r12, pc}
	bx	lr
	nop
.endfunc

/*
//function for a TOOM-COOK in assembly
//*************************************************
// toom_cook_4way_mem_asm(&a1, &b1, &result);
//*************************************************
.global toom_cook_4way_mem_asm
.func toom_cook_4way_mem_asm, toom_cook_4way_mem_asm
.type toom_cook_4way_mem_asm, %function
toom_cook_4way_mem_asm:
	stmdb	sp!, {r3, r4, r5, r6, r7, r8, r9, r10, r11, r12, lr}
	subw	sp, sp, #3076
	mov 	r3, sp
	mov 	r11, r2
	add 	r2, r3, #1792
//r0=&a1[i] r1=b1 r2=&w_m[896+i]
//loop1
	add 	r12, r0, #128
toom_cook_4way_mem_asm_loop1:
	ldr 	r6, [r0, #384]
	ldr 	r4, [r0, #128]
	ldr 	r5, [r0, #256]
	ldr 	r3, [r0], #4
	cmp 	r0, r12
	uadd16 	r8, r4, r6
	uadd16 	r7, r3, r5
	uadd16 	r9, r7, r8
	usub16 	r10, r7, r8
	str 	r9, [r2, #128]
	str 	r10, [r2, #256]
	uadd16 	r7, r3, r3
	uadd16	r7, r7, r7
	uadd16	r7, r7, r5
	uadd16 	r7, r7, r7
	uadd16 	r8, r4, r4
	uadd16 	r8, r8, r8
	uadd16 	r8, r8, r6
	uadd16 	r9, r7, r8
	usub16 	r10, r7, r8
	str 	r9, [r2, #384]
	str 	r10, [r2, #512]
	uadd16 	r7, r6, r6
	uadd16 	r7, r7, r5
	uadd16 	r7, r7, r7
	uadd16 	r7, r7, r4
	uadd16 	r7, r7, r7
	uadd16 	r7, r7 ,r3
	str 	r7, [r2], #4
	bne 	toom_cook_4way_mem_asm_loop1
//loop2
	add 	r12, r1, #128
	add 	r2, r2, #512
toom_cook_4way_mem_asm_loop2:
	ldr 	r6, [r1, #384]
	ldr 	r4, [r1, #128]
	ldr 	r5, [r1, #256]
	ldr 	r3, [r1], #4
	cmp 	r1, r12
	uadd16 	r8, r4, r6
	uadd16 	r7, r3, r5
	uadd16 	r9, r7, r8
	usub16 	r10, r7, r8
	str 	r9, [r2, #128]
	str 	r10, [r2, #256]
	uadd16 	r7, r3, r3
	uadd16	r7, r7, r7
	uadd16	r7, r7, r5
	uadd16 	r7, r7, r7
	uadd16 	r8, r4, r4
	uadd16 	r8, r8, r8
	uadd16 	r8, r8, r6
	uadd16 	r9, r7, r8
	usub16 	r10, r7, r8
	str 	r9, [r2, #384]
	str 	r10, [r2, #512]
	uadd16 	r7, r6, r6
	uadd16 	r7, r7, r5
	uadd16 	r7, r7, r7
	uadd16 	r7, r7, r4
	uadd16 	r7, r7, r7
	uadd16 	r7, r7 ,r3
	str 	r7, [r2], #4
	bne 	toom_cook_4way_mem_asm_loop2
//r0=&a1[i+64] r1=&b1[i+64] r2=&w_m[896+320+64+i] r3-r10 r12
	mov 	r9, r11
	mov 	r11, r0
	mov 	r12, r1
	sub 	r0, r2, #2560
	mov 	r1, #0
	mov 	r2, #1792
	bl 		memset
	mov 	r10, r0
	sub 	r0, r11, #128
	sub 	r1, r12, #128
	add		r2, r10, #1536
	mov 	r3, #32
	bl 		karatsuba_asm
	add 	r0, r10, #2304
	add 	r1, r10, #2944
	add 	r2, r10, #1280
	mov 	r3, #32
	bl 		karatsuba_asm
	add 	r0, r10, #2176
	add 	r1, r10, #2816
	add 	r2, r10, #1024
	mov 	r3, #32
	bl 		karatsuba_asm
	add 	r0, r10, #2048
	add 	r1, r10, #2688
	add 	r2, r10, #768
	mov 	r3, #32
	bl 		karatsuba_asm
	add 	r0, r10, #1920
	add 	r1, r10, #2560
	add 	r2, r10, #512
	mov 	r3, #32
	bl 		karatsuba_asm
	add 	r0, r10, #1792
	add 	r1, r10, #2432
	add 	r2, r10, #256
	mov 	r3, #32
	bl 		karatsuba_asm
	add 	r0, r11, #256
	add 	r1, r12, #256
	mov 	r2, r10
	mov 	r3, #32
	bl 		karatsuba_asm
//r0=i r9=result r10=w_m
//loop3
	movw	r14, #0x1fff
	mov 	r0, #0
	mov 	r11, #128
	add 	r14, r14, r14, lsl #16
toom_cook_4way_mem_asm_loop3:
	ldr 	r5, [r10, #1024]
	ldr 	r3, [r10, #512]
	ldr 	r2, [r10, #256]
	ldr 	r7, [r10, #1536]
	ldr 	r6, [r10, #1280]
	ldr 	r4, [r10, #768]
	ldr 	r1, [r10], #4
	uadd16	r2, r2, r5
	usub16	r6, r6, r5
	usub16	r4, r4, r3
	usub16	r5, r5, r1
	lsr 	r11, r7, #16
	lsl		r12, r7, #6
	pkhbt	r8, r12, r11, lsl #22
	usub16	r5, r5, r8
	lsl		r11, r4, #16
	lsr	 	r12, r4, #1
	pkhtb 	r4, r12, r11, asr #17
	uadd16	r3, r3, r4
	lsr 	r11, r5, #16
	lsl		r12, r5, #1
	pkhbt	r8, r12, r11, lsl #17
	uadd16	r5, r6, r8
	lsr 	r11, r3, #16
	lsl		r12, r3, #6
	pkhbt	r8, r12, r11, lsl #22
	uadd16	r8, r3, r8
	usub16	r2, r2, r8
	usub16	r3, r3, r7
	usub16 	r3, r3, r1
	mov 	r8, #45
	smulbb 	r11, r3, r8
	smultb 	r12, r3, r8
	pkhbt 	r8, r11, r12, lsl #16
	uadd16	r2, r2, r8
	lsr 	r11, r3, #16
	lsl		r12, r3, #3
	pkhbt	r8, r12, r11, lsl #19
	usub16	r5, r5, r8
	movw 	r8, #43691
	smulbb 	r11, r5, r8
	smultb 	r12, r5, r8
	pkhbt 	r8, r11, r12, lsl #16
	lsl		r11, r8, #16
	lsr	 	r12, r8, #3
	pkhtb 	r5, r12, r11, asr #19
	uadd16 	r6, r2, r6
	lsr 	r11, r4, #16
	lsl		r12, r4, #4
	pkhbt	r8, r12, r11, lsl #20
	uadd16	r2, r2, r8
	movw 	r8, #36409
	smulbb 	r11, r2, r8
	smultb 	r12, r2, r8
	pkhbt 	r8, r11, r12, lsl #16
	lsl		r11, r8, #16
	lsr	 	r12, r8, #1
	pkhtb 	r2, r12, r11, asr #17
	usub16	r3, r3, r5
	uadd16	r4, r4, r2
	mov		r8, #0
	usub16	r4, r8, r4
	mov 	r8, #30
	smulbb 	r11, r2, r8
	smultb 	r12, r2, r8
	pkhbt 	r8, r11, r12, lsl #16
	usub16	r6, r8, r6
	movw 	r8, #61167
	smulbb 	r11, r6, r8
	smultb 	r12, r6, r8
	pkhbt 	r8, r11, r12, lsl #16
	lsl		r11, r8, #16
	lsr	 	r12, r8, #1
	pkhtb 	r6, r12, r11, asr #17
//r0=i r9=result r10=w_m
	usub16	r7, r7, r3
	usub16	r6, r6, r2
	usub16	r5, r5, r1
	ldr 	r3, [r9]
	ldr 	r2, [r9, #128]
	ldr 	r1, [r9, #256]
	uadd16	r7, r7, r3
	uadd16	r6, r6, r2
	uadd16	r5, r5, r1
	#ubfx	r7, r7, #0, #13
	#ubfx	r6, r6, #0, #13
	#ubfx	r5, r5, #0, #13
	and		r7, r7, r14
	and		r6, r6, r14
	and 	r5, r5, r14
	cmp 	r0, #63
	add		r0, r0, #2
	sub 	r1, r9, #128
	it 		ls
	addls	r1, r1, #512
	ldr 	r8, [r1]
	str 	r5, [r9, #256]
	str 	r6, [r9, #128]
	str 	r7, [r9], #4
	usub16	r3, r8, r4
	it 		ls
	uadd16ls r3, r8, r4
	cmp 	r0, #126
	#ubfx 	r3, r3, #0, #13
	and 	r3, r3, r14
	str 	r3, [r1]
	bne		toom_cook_4way_mem_asm_loop3
//for(i=0;i<2*(SABER_N/4)-1;i++){ last iteration only halfword!!!
	ldrh 	r5, [r10, #1024]
	ldrh 	r3, [r10, #512]
	ldrh 	r2, [r10, #256]
	ldrh 	r7, [r10, #1536]
	ldrh 	r6, [r10, #1280]
	ldrh 	r4, [r10, #768]
	ldrh 	r1, [r10], #4
	uadd16	r2, r2, r5
	usub16	r6, r6, r5
	usub16	r4, r4, r3
	usub16	r5, r5, r1
	lsl 	r8, r7, #6
	usub16	r5, r5, r8
	lsr		r4, r4, #1
	uadd16	r3, r3, r4
	lsl		r8, r5, #1
	uadd16	r5, r6, r8
	lsl		r8, r3, #6
	uadd16	r8, r3, r8
	usub16	r2, r2, r8
	usub16	r3, r3, r7
	usub16 	r3, r3, r1
	mov 	r8, #45
	smulbb 	r8, r3, r8
	uadd16	r2, r2, r8
	lsl		r8, r3, #3
	usub16	r5, r5, r8
	movw 	r8, #43691
	smulbb 	r8, r5, r8
	lsr		r5, r8, #3
	uadd16 	r6, r2, r6
	lsl		r8, r4, #4
	uadd16	r2, r2, r8
	movw 	r8, #36409
	smulbb 	r8, r2, r8
	lsr 	r2, r8, #1
	usub16	r3, r3, r5
	uadd16	r4, r4, r2
	mov		r8, #0
	usub16	r4, r8, r4
	mov 	r8, #30
	smulbb 	r8, r2, r8
	usub16	r6, r8, r6
	movw 	r8, #61167
	smulbb 	r8, r6, r8
	lsr		r6, r8, #1
//r0=i r9=result r10=w_m
	usub16	r7, r7, r3
	usub16	r6, r6, r2
	usub16	r5, r5, r1
	ldrh 	r3, [r9]
	ldrh 	r2, [r9, #128]
	ldrh 	r1, [r9, #256]
	uadd16	r7, r7, r3
	uadd16	r6, r6, r2
	uadd16	r5, r5, r1
	ubfx	r7, r7, #0, #13
	ubfx	r6, r6, #0, #13
	ubfx	r5, r5, #0, #13
	cmp 	r0, #63
	add		r0, r0, #2
	sub 	r1, r9, #128
	it 		ls
	addls	r1, r1, #512
	ldr 	r8, [r1]
	strh 	r5, [r9, #256]
	strh 	r6, [r9, #128]
	strh 	r7, [r9], #4
	usub16	r3, r8, r4
	it 		ls
	uadd16ls r3, r8, r4
	cmp 	r0, #126
	ubfx 	r3, r3, #0, #13
	strh 	r3, [r1]
exit_toom_cook_4way_mem_asm:
	addw	sp, sp, #3076
	ldmia.w	sp!, {r3, r4, r5, r6, r7, r8, r9, r10, r11, r12, pc}
	bx	lr
	nop
.endfunc
*/

