COPTS = -std=c99 -g -march=native -O3 -DNDEBUG -Wall -Wextra -Werror

all: lib ber_trials throughput_trials

lib:
	mkdir -p build
	gcc $(COPTS) -Iinclude/ -c src/ldpc_decoder.c -o build/ldpc_decoder.o
	gcc $(COPTS) -Iinclude/ -c src/ldpc_encoder.c -o build/ldpc_encoder.o
	gcc $(COPTS) -Iinclude/ -c src/ldpc_codes.c -o build/ldpc_codes.o
	ar -rcs build/liblabradorldpc.a build/ldpc_decoder.o build/ldpc_encoder.o build/ldpc_codes.o

ber_trials: lib
	gcc $(COPTS) -std=gnu99 -fopenmp -Iinclude/ bin/ber_trials.c -Lbuild -llabradorldpc -lm -o build/ber_trials

throughput_trials: lib
	gcc $(COPTS) -std=gnu99 -fopenmp -Iinclude/ bin/throughput_trials.c -Lbuild -llabradorldpc -lm -o build/throughput_trials

test: lib
	gcc $(COPTS) -std=gnu99 -fopenmp -Iinclude/ bin/test.c -Lbuild -llabradorldpc -lm -o build/test
	./build/test
