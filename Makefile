build:
	gcc -O3 -march=native -flto -fwhole-program -ffast-math -lm -Isrc/include -o scimark $(wildcard src/*.c)

