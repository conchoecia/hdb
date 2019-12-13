CC = g++
CPPFLAGS = --std=c++11 -Wall -I../src/ -O3
CPPFLAGSDEBUG = --std=c++11 -Wall -I../src/


#all: gen_rand_HKC gen_canon gen_HKC_rule gen_HKC_nodelete_rule print_DBG
all: gen_HKC gen_HKC_nodelete

#print_DBG_rule: print_DBG
#	$(CC) $(CPPFLAGS) print_DBG.cpp -o print_DBG

gen_HKC: gen_HKC.cpp
	$(CC) $(CPPFLAGS) gen_HKC.cpp -o gen_HKC

gen_HKC_nodelete: gen_HKC_nodelete.cpp
	$(CC) $(CPPFLAGSDEBUG) gen_HKC_nodelete.cpp -o gen_HKC_nodelete

#gen_rand_HKC: gen_random_HKC
#	$(CC) $(CPPFLAGS) gen_random_HKC.cpp -o gen_random_HKC
#
#gen_canon: canon
#	$(CC) $(CPPFLAGS) canon.cpp -o canon
