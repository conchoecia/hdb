#include <iostream>
#include <random>
#include <set>
#include <string>
#include "hdb.hpp"

/* this program generates a string that contains unique kmers along its entirety
   arg 1 is size k
   arg 2 is how many kmers long to make the string
   */
int main (int argc, char* argv[]){
  uint32_t k = atoi(argv[1]);
  if (k > 32){
    throw std::runtime_error("Goodbye. k must be 32 or smaller.\n");
  }
  if (k%2 == 0){
    throw std::runtime_error("Goodbye. k must be an odd number.\n");
  }

  uint32_t num_kmers = atoi(argv[2]);
  uint32_t str_length = k+num_kmers-1;
  std::cout << "k is " << k << std::endl;
  std::cout << "num kmers is " << k << std::endl;
  std::cout << "length of string is " << str_length << std::endl;

  uint64_t kmer;
  uint64_t temp;
  uint64_t mask = (1ULL << (k*2))-1;
  std::default_random_engine generator (time(0));
  std::uniform_int_distribution<uint64_t> uniform (0,3);
  uint32_t l = 0;
  std::set<uint64_t> visited;
  std::string string_out;
  uint64_t iloop_counter = 0;
  while(string_out.size() < str_length){
    temp = uniform(generator);
    kmer =  ( (kmer << 2) | temp ) & mask;
    if (l < k){
     string_out = string_out + int_base_to_str[temp];
     l++;
     if (l == k){
        uint64_t canonical = canon(kmer, k);
        visited.insert(canonical);
     }
    }
    else {
      if (iloop_counter < 4000){
        uint64_t canonical = canon(kmer, k);
        if (visited.count(canonical) == 0){ //haven't encountered this kmer before
          string_out = string_out + int_base_to_str[temp];
        }
        visited.insert(canonical);
        iloop_counter++;
      }
      else {
        throw std::runtime_error("Goodbye. couldn't make a unique string. Please try again or increase k/decrease num_kmers.\n");
      }
    }
  }

  //now we just verify that all of the kmers are unique
  visited.clear();
  for (uint32_t i = 0; i < string_out.size() - k + 1; i++){
    temp = kmer_to_uint64(string_out.substr(i, k), k);
    if (visited.count(temp) == 1){
      throw std::runtime_error("Goodbye. Thought we had made a unique string but verification failed.\n");
    }
    visited.insert(temp);
  }

  std::cout << "final string is: " << string_out << std::endl;

  //now generate a kmer in front and at the end
  std::string kmer_front = string_out.substr(0, k);
  std::vector<std::string> vec = {"A", "C", "G", "T"};
  for (const auto & n: vec){
    kmer_front.replace(0,1,n);
    temp = kmer_to_uint64(kmer_front, k);
    if(visited.count(temp) == 0){
      std::cout << "A branching 5p kmer is: " << kmer_front << std::endl;
      break;
    }
  }
  //now generate a kmer at the end
  std::string kmer_end = string_out.substr(string_out.size()-k, k);
  for (const auto & n: vec){
    kmer_end.replace(kmer_end.size()-1,1,n);
    temp = kmer_to_uint64(kmer_end, k);
    if(visited.count(temp) == 0){
      std::cout << "A branching 3p kmer is: " << kmer_end << std::endl;
      break;
    }
  }
}
