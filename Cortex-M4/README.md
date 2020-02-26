
Implementation of the CCA-secure module lattice-based key encapsulation Saber suitable for ARM Cortex-M4 processors. The goal of this implementation is to get high performance, competing with similar state-of-the-art schemes, while keeping a reasonably low memory utilization. In addition to the cryptosystem as presented to the standardization contest by NIST, we introduced a slightly modified version called Saber# that enables a further decrease of the memory requirements with no cost (even a slightly improve) in performance.

## Setup

After cloning or downloading the repository, it is necessary to add the libraries to work with the microcontroller. The most straightforward way to do it is to run the following commands since the library has been added as a submodule. More information about this library can be found in https://github.com/libopencm3/libopencm3.

```sh
$ git submodule init
$ git submodule update
```

#### Requirements

- The cross compiler for ARM Cortex-M series is required in order to build the executable. The toolchain can be found in https://launchpad.net/gcc-arm-embedded.

- The stlink software tool is required in order to flash the executables into the board. The source code for this tool as well as installation instructions can be found in https://github.com/texane/stlink/blob/master/doc/compiling.md.

#### Benchmarks

In order to communicate with the board through the serial port it is necessary a UART-USB connector. The `RX` and `TX` pins of the connector must be connected to the `PA2` and `PA3` pins of the board, respectively. Currently, there are two commands of interest for benchmarking purposes:

- `make saberSpeed` (alternatively `make saber-sharpSpeed`). This command will print out the number of clock cycles required to perform each of the operations of the key exchange mechanism (KEM), namely key generation, encapsulation and decapsulation.

- `make saberMemsize` (alternatively `make saber-sharpMemsize`). This command will print out the memory utilization of each of the operations of the KEM.

The table below summarizes the results of the available implementations.

Scheme    | key_gen [clock cycles]/[bytes] | encaps [clock cycles]/[bytes] | decaps [clock cycles]/[bytes]
--------- | ------------------------------ | ----------------------------- | -----------------------------
Saber     | 1.165.845 / 6.931              | 1.530.745 / 7.019             | 1.635.720 / 8.115
Saber#    | 1.162.680 / 6.923              | 1.528.061 / 5.987             | 1.633.035 / 7.163

