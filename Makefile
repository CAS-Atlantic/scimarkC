CC ?= gcc
SRC := $(wildcard src/*.c) scimark4.c
CFLAGS := -O3 -flto -fwhole-program -lm -Isrc/include

.PHONY: scimark

scimark: 
	$(CC) $(CFLAGS) -o $@ $(SRC)

