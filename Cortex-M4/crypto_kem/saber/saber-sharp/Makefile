CC     = arm-none-eabi-gcc
CFLAGS = -Wall -Wextra -O3 -mthumb -mcpu=cortex-m4 -mfloat-abi=hard -mfpu=fpv4-sp-d16 
AR     = arm-none-eabi-gcc-ar

HEADERS = api.h cbd.h kem.h pack_unpack.h poly.h poly_mul.h encode.h SABER_indcpa.h SABER_params.h verify.h
SOURCES = cbd.c kem.c pack_unpack.c poly.c poly_mul.c encode.c SABER_indcpa.c verify.c asm_funcs.S
OBJECTS = cbd.o kem.o pack_unpack.o poly.o poly_mul.o encode.o SABER_indcpa.o verify.o asm_funcs.o

libpqm4.a: $(OBJECTS)
	$(AR) rcs $@ $(OBJECTS)

%.o: %.c $(HEADERS)
	$(CC) -I$(INCPATH) $(CFLAGS) -c -o $@ $<

%.o: %.S $(HEADERS)
	$(CC) -I$(INCPATH) $(CFLAGS) -c -o $@ $<

.PHONY: clean

clean:
	-rm -f libpqm4.a
	-rm -f *.o
