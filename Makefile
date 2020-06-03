CC ?= gcc
SRC := $(wildcard src/*.c) scimark4.c
CFLAGS := -O3 -fno-unroll-loops -Isrc/include

.PHONY: scimark

scimark: 
	$(CC) $(CFLAGS) -o $@ $(SRC) -lm

