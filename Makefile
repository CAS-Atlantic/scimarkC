# SciMark2: A Java numerical benchmark measuring performance
# of computational kernels for FFTs, Monte Carlo simulation,
# sparse matrix computations, Jacobi SOR, and dense LU matrix
# factorizations.
#
# Authors: Roldan Pozo (pozo@nist.gov) and Bruce Miller (bruce.miller@nist.gov)
#
# Downloaded from: https://math.nist.gov/scimark2/download.html
#
# Modified by: Aaron Graham (aaron.graham@unb.ca, aarongraham9@gmail.com) and
#               Jean-Philippe Legault (jlegault@unb.ca, jeanphilippe.legault@gmail.com)
#                for the Centre for Advanced Studies - Atlantic (CAS-Atlantic) at the
#                 Univerity of New Brunswick in Fredericton, New Brunswick, Canada

CC ?= gcc
SRC := $(wildcard src/*.c)
CFLAGS := -O3 -fno-unroll-loops -Isrc/include
BM := fft sor lu monte_carlo sparse 
.PHONY: scimark $(BM)

all: $(BM)

$(BM): $(SRC)
	mkdir -p bin
	$(CC) $(CFLAGS) -o bin/$@ $(SRC) src/$@/kernel.c -lm

