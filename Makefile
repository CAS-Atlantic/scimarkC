CC ?= gcc
SRC := $(wildcard src/*.c) scimark4.c
CFLAGS := -O3 -march=native -flto -fwhole-program -ffast-math -Isrc/include

.PHONY: scimark

scimark: 
	$(CC) $(CFLAGS) -o $@ $(SRC) -lm

