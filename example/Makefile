CC = g++
CPPFLAGS = --std=c++11 -Wall -I../src/

all: gen_rand_HKC gen_canon

gen_rand_HKC: gen_random_HKC
	$(CC) $(CPPFLAGS) gen_random_HKC.cpp -o gen_random_HKC

gen_canon: canon
	$(CC) $(CPPFLAGS) canon.cpp -o canon
