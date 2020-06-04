CC ?= gcc
SRC := $(wildcard src/*.c)
CFLAGS := -O3 -fno-unroll-loops -Isrc/include
BM := FFT SOR LU MonteCarlo SparseMatMult 
.PHONY: scimark $(BM)

all: $(BM)

$(BM): $(SRC)
	mkdir -p bin
	$(CC) $(CFLAGS) -o bin/$@ $(SRC) src/$@/kernel.c -lm

