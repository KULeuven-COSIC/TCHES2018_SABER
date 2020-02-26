#STFLASH=$(HOME)/stlink/st-flash
STFLASH=st-flash

OPENCM3DIR  = ./libopencm3
OPENCM3NAME = opencm3_stm32f4
OPENCM3FILE = $(OPENCM3DIR)/lib/lib$(OPENCM3NAME).a
LDSCRIPT    = $(OPENCM3DIR)/lib/stm32/f4/stm32f405x6.ld

PREFIX     ?= arm-none-eabi
CC          = $(PREFIX)-gcc
LD          = $(PREFIX)-gcc
OBJCOPY     = $(PREFIX)-objcopy
OBJDUMP     = $(PREFIX)-objdump
GDB         = $(PREFIX)-gdb

ARCH_FLAGS  = -mthumb -mcpu=cortex-m4 -mfloat-abi=hard -mfpu=fpv4-sp-d16
DEFINES     = -DSTM32F4
OBJS        = obj/stm32f4_wrapper.o obj/fips202.o obj/keccakf1600.o
RANDOMBYTES = obj/randombytes.o

CFLAGS     += -O3 \
              -Wall -Wextra -Wimplicit-function-declaration \
              -Wredundant-decls -Wmissing-prototypes -Wstrict-prototypes \
              -Wundef -Wshadow \
              -I$(OPENCM3DIR)/include \
              -fno-common $(ARCH_FLAGS) -MD $(DEFINES)
LDFLAGS    += --static -Wl,--start-group -lc -lgcc -lnosys -Wl,--end-group \
              -T$(LDSCRIPT) -nostartfiles -Wl,--gc-sections \
               $(ARCH_FLAGS) -L$(OPENCM3DIR)/lib


KEMLIBS=$(wildcard crypto_kem/*/*)
KEMTESTS=$(patsubst %,bin/%,$(patsubst %,%_test.bin,$(subst /,_,$(KEMLIBS))))
KEMTESTVECTORS=$(patsubst %,bin/%,$(patsubst %,%_testvectors.bin,$(subst /,_,$(KEMLIBS))))
KEMSPEEDS=$(patsubst %,bin/%,$(patsubst %,%_speed.bin,$(subst /,_,$(KEMLIBS))))
KEMMEMSIZES=$(patsubst %,bin/%,$(patsubst %,%_memsize.bin,$(subst /,_,$(KEMLIBS))))

OWNDIR=$(shell pwd)
INCPATH=$(OWNDIR)/common

all: tests testvectors speeds memsizes

libs: $(KEMLIBS) 
tests: libs $(KEMTESTS) 
testvectors: libs $(KEMTESTVECTORS) 
speeds: libs $(KEMSPEEDS) 
memsizes: libs $(KEMMEMSIZES)

# Memory utilization for Saber
%Memsize: bin/crypto_kem_saber_%_memsize.bin
	@echo "Writing Memsize binary to devices" 
	$(STFLASH) write $^ 0x8000000
	python3 ./hostside/host_unidirectional.py

# Test program for Saber
%Test: bin/crypto_kem_saber_%_test.bin
	@echo "Writing Test binary to devices" 
	$(STFLASH) write $^ 0x8000000
	python3 ./hostside/host_unidirectional.py

# Cycle count for Saber
%Speed: bin/crypto_kem_saber_%_speed.bin
	@echo "Writing Speed binary to devices"
	$(STFLASH) write $^ 0x8000000
	python3 ./hostside/host_unidirectional.py



.PHONY: force
export INCPATH

$(KEMLIBS): force
	make -C $@

$(SIGNLIBS): force
	make -C $@


bin/%.bin: elf/%.elf
	mkdir -p bin
	$(OBJCOPY) -Obinary $^ $@



elf/crypto_kem_%_test.elf: $(OBJS) $(RANDOMBYTES) $(LDSCRIPT) obj/$(patsubst %,crypto_kem_%_test.o,%) $(OPENCM3FILE)
	mkdir -p elf
	$(LD) -o $@ \
	$(patsubst elf/%.elf,obj/%.o,$@) \
	$(patsubst %test.elf,%libpqm4.a,$(patsubst elf/crypto/kem%,crypto_kem/%,$(subst _,/,$@))) \
	$(OBJS) $(RANDOMBYTES) $(LDFLAGS) -l$(OPENCM3NAME) -lm

elf/crypto_kem_%_testvectors.elf: $(OBJS) $(LDSCRIPT) obj/$(patsubst %,crypto_kem_%_testvectors.o,%) $(OPENCM3FILE)
	mkdir -p elf
	$(LD) -o $@ \
	$(patsubst elf/%.elf,obj/%.o,$@) \
	$(patsubst %testvectors.elf,%libpqm4.a,$(patsubst elf/crypto/kem%,crypto_kem/%,$(subst _,/,$@))) \
	$(OBJS) $(LDFLAGS) -l$(OPENCM3NAME) -lm

elf/crypto_kem_%_speed.elf: $(OBJS) $(LDSCRIPT) obj/$(patsubst %,crypto_kem_%_speed.o,%) $(OPENCM3FILE)
	mkdir -p elf
	$(LD) -o $@ \
	$(patsubst elf/%.elf,obj/%.o,$@) \
	$(patsubst %speed.elf,%libpqm4.a,$(patsubst elf/crypto/kem%,crypto_kem/%,$(subst _,/,$@))) \
	$(OBJS) $(RANDOMBYTES) $(LDFLAGS) -l$(OPENCM3NAME) -lm

elf/crypto_kem_%_memsize.elf: $(OBJS) $(LDSCRIPT) obj/$(patsubst %,crypto_kem_%_memsize.o,%) $(OPENCM3FILE)
	mkdir -p elf
	$(LD) -o $@ \
	$(patsubst elf/%.elf,obj/%.o,$@) \
	$(patsubst %memsize.elf,%libpqm4.a,$(patsubst elf/crypto/kem%,crypto_kem/%,$(subst _,/,$@))) \
	$(OBJS) $(RANDOMBYTES) $(LDFLAGS) -l$(OPENCM3NAME) -lm




obj/crypto_kem_%_test.o: crypto_kem/test.c $(patsubst %,%/api.h,$(patsubst %,crypto_kem/%,$(subst _,/,$%)))
	mkdir -p obj 
	$(CC) $(CFLAGS) -o $@ -c $< \
	-I$(patsubst %test.o,%,$(patsubst obj/%,%,$(subst crypto/kem,crypto_kem,$(subst _,/,$@)))) \
	-I./common/

obj/crypto_kem_%_testvectors.o: crypto_kem/testvectors.c $(patsubst %,%/api.h,$(patsubst %,crypto_kem/%,$(subst _,/,$%)))
	mkdir -p obj 
	$(CC) $(CFLAGS) -o $@ -c $< \
	-I$(patsubst %testvectors.o,%,$(patsubst obj/%,%,$(subst crypto/kem,crypto_kem,$(subst _,/,$@)))) \
	-I./common/

obj/crypto_kem_%_speed.o: crypto_kem/speed.c $(patsubst %,%/api.h,$(patsubst %,crypto_kem/%,$(subst _,/,$%)))
	mkdir -p obj 
	$(CC) $(CFLAGS) -o $@ -c $< \
	-I$(patsubst %speed.o,%,$(patsubst obj/%,%,$(subst crypto/kem,crypto_kem,$(subst _,/,$@)))) \
	-I./common/

obj/crypto_kem_%_memsize.o: crypto_kem/memsize.c $(patsubst %,%/api.h,$(patsubst %,crypto_kem/%,$(subst _,/,$%)))
	mkdir -p obj 
	$(CC) $(CFLAGS) -o $@ -c $< \
	-I$(patsubst %memsize.o,%,$(patsubst obj/%,%,$(subst crypto/kem,crypto_kem,$(subst _,/,$@)))) \
	-I./common/




obj/randombytes.o: common/randombytes.c
	mkdir -p obj 
	$(CC) $(CFLAGS) -o $@ -c $^

obj/stm32f4_wrapper.o:  common/stm32f4_wrapper.c $(OPENCM3FILE)
	mkdir -p obj 
	$(CC) $(CFLAGS) -o $@ -c $^

obj/fips202.o:  common/fips202.c
	mkdir -p obj 
	$(CC) $(CFLAGS) -o $@ -c $^

obj/keccakf1600.o:  common/keccakf1600.S
	mkdir -p obj 
	$(CC) $(CFLAGS) -o $@ -c $^

$(OPENCM3FILE):
	@if [ ! "`ls -A $(OPENCM3_DIR)`" ] ; then \
		printf "######## ERROR ########\n"; \
		printf "\tlibopencm3 is not initialized.\n"; \
		printf "\tPlease run (in the root directory):\n"; \
		printf "\t$$ git submodule init\n"; \
		printf "\t$$ git submodule update\n"; \
		printf "\tbefore running make.\n"; \
		printf "######## ERROR ########\n"; \
		exit 1; \
		fi
	make -C $(OPENCM3DIR)



.PHONY: clean libclean

clean:
	find . -name \*.o -type f -exec rm -f {} \;
	find . -name \*.d -type f -exec rm -f {} \;
	find crypto_kem -name \*.a -type f -exec rm -f {} \;
	rm -rf elf/
	rm -rf bin/
	rm -rf obj/

libclean:
	make -C $(OPENCM3DIR) clean
