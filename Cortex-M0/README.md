
Implementation of the CCA-secure module lattice-based key encapsulation Saber suitable for ARM Cortex-M0 processors. The goal of this implementation is to reduce the required memory as much as possible. In addition to the cryptosystem as presented to the standardization contest by NIST, we introduced a slightly modified version called Saber# that enables a further decrease of the memory requirements.

## Setup

#### Requirements

- The cross compiler for ARM Cortex-M series is required in order to build the executable. The toolchain can be found in https://launchpad.net/gcc-arm-embedded and the path to the libraries for the linker should be `/usr/arm-none-eabi/lib/armv6-m`. If you choose a custom location for your installation, the path to the libraries should be updated accordingly in the Makefiles.

- The JLink debugging software tool is required in order to flash the executable into the board. The package can be found in https://www.segger.com/downloads/jlink/ and it is necessary to update the path in the line 42 of the Makefiles to match your installation.

#### Benchmarks

After plugging the XMC2Go Infineon board it is enough to run the script `evaluate.py` in the corresponding directory to get all timing and memory results of that scheme. The following table shows the data for the available Saber and Saber# implementations.

Scheme    | key_gen [clock cycles]/[bytes] | encaps [clock cycles]/[bytes] | decaps [clock cycles]/[bytes]
--------- | ------------------------------ | ----------------------------- | -----------------------------
Saber     | 4.786.727 / 5.031              | 6.328.445 / 5.119             | 7.509.801 / 6.215
Saber#    | 4.782.306 / 5.007              | 6.325.718 / 4.079             | 7.507.074 / 5.175
