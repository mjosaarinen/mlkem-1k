#   Makefile for the minimal x86-64 ML-KEM (mlkem.S) and its test harness.
#   make test   -- build and run ML-KEM-512, -768 and -1024 tests
#   make size   -- code size of each parameter set and per-function sizes
#   make kat    -- regenerate the vector block (paste it into mlkem_test.c)

CC      = gcc
CFLAGS  = -O2 -Wall -Wextra -Wpedantic -Wno-overlength-strings
BINS    = mlkem512_test mlkem768_test mlkem1024_test

all:    $(BINS)

mlkem512_test:  mlkem_test.c mlkem.S
	$(CC) $(CFLAGS) -DMLKEM_K=2 -o $@ mlkem_test.c mlkem.S

mlkem768_test:  mlkem_test.c mlkem.S
	$(CC) $(CFLAGS) -DMLKEM_K=3 -o $@ mlkem_test.c mlkem.S

mlkem1024_test: mlkem_test.c mlkem.S
	$(CC) $(CFLAGS) -DMLKEM_K=4 -o $@ mlkem_test.c mlkem.S

test:   $(BINS)
	./mlkem512_test && ./mlkem768_test && ./mlkem1024_test

size:
	@for k in 2 3 4; do \
		$(CC) -c -DMLKEM_K=$$k mlkem.S -o mlkem_k$$k.o && \
		size mlkem_k$$k.o | tail -1; done
	@nm -S --size-sort mlkem_k2.o
	@rm -f mlkem_k*.o

kat:
	python3 gen_kat.py ../py-acvp-pqc/json-copy

clean:
	rm -f $(BINS) *.o

.PHONY: all test size kat clean
