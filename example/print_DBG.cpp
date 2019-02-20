#include <iostream>
#include <fstream>
#include <unistd.h> //needed for kseq.h and getopt
#include "khash.h"
#include "gen_HKC.hpp"
#include "hdb.hpp"

int main(int argc, char **argv) {
  Vars vars = process_cl_args(argc, argv);
  DBG G = DBG(vars.k, 1);
  //std::cout << " - Reading kmers from stdin." << std::endl;
  parse_kmer_dump(G, std::cin, vars.k, 0 );
  std::cout << "source\ttarget\n";
  G.print_graph();
  return 0;
}

void help(void){
  std::cout << "help statement needs to be rewritten\n";
}

Vars process_cl_args(int argc, char **argv){
  Vars vars;
  int c;
  if ( argc == 1 ) {
    help();
    exit(1);
  }
  //printf("argv: %s\n", argv);
  //printf("argc: %s\n", argc);
  while( (c=getopt( argc, argv, "k:m:M:o:" )) != -1 ) {
    switch(c) {
    case 'o' :
      vars.out_prefix = optarg;
      break;
    case 'k' :
      vars.k = static_cast<uint32_t>(atoi(optarg));
      break;
    case 'h' :
      help();
    }
  }
  // check mandatory
  // first check that k is not zero, then that it is odd
  if (vars.k == 0) {
    printf("-k is mandatory.\n");
    printf("    Goodbye.\n");
    exit(1);
  }
  if (vars.k % 2 == 0){
    printf("-k must be an odd number.\n");
    printf( "    %d is not an odd number.\n", vars.k);
    printf( "    Goodbye.\n");
    exit(1);
  }
  return(vars);
}
