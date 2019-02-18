#include <iostream>
#include <string>
#include "hdb.hpp"

int main(int argc, char *argv[]){
  std::string kmer(argv[1]);
  uint32_t k = kmer.size();
  uint64_t key = kmer_to_uint64(kmer, k);
  uint64_t can = canon(key, k);
  std::cout << "canonical kmer of " << kmer << ": " << uint64_to_kmer(can, k) << '\n';
  return 0;
}
