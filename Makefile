
CC ?= gcc

OBJ := \
	build/FFT.o \
	build/kernel.o \
	build/Stopwatch.o \
	build/Random.o \
	build/SOR.o \
	build/SparseCompRow.o \
	build/array.o \
	build/MonteCarlo.o \
	build/LU.o \
	build/scimark4.o

CFLAGS := -O3 -march=native -flto -fwhole-program -ffast-math -lm -Isrc/include

.SUFFIXES: .o .c
.PHONY: build

build:
	mkdir -p build
	$(MAKE) scimark

build/%.o: src/%.c
	$(CC) $(CFLAGS) -c $< -o $@

scimark: $(OBJ)
	$(CC) $(CFLAGS) -o scimark $(OBJ)

clean:
	rm -Rf build

wipe: clean
	rm scimark

