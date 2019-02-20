CC = g++
CPPFLAGS = --std=c++11 -Wall -I../src/ -O3

all: gen_rand_HKC gen_canon gen_HKC_rule print_DBG

print_DBG_rule: print_DBG
	$(CC) $(CPPFLAGS) print_DBG.cpp -o print_DBG

gen_HKC_rule: gen_HKC
	$(CC) $(CPPFLAGS) gen_HKC.cpp -o gen_HKC

gen_rand_HKC: gen_random_HKC
	$(CC) $(CPPFLAGS) gen_random_HKC.cpp -o gen_random_HKC

gen_canon: canon
	$(CC) $(CPPFLAGS) canon.cpp -o canon
