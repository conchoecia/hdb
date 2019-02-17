#ifndef DE_BRUIJN_H
#define DE_BRUIJN_H

#include <fstream>
#include <map>
#include <stdexcept> //std::runtime_error
#include <utility>
#include <vector>
#include "khash.h"
#include <iostream>

template<class T>
uint64_t canon(const uint64_t & kmer, T k);

/* The de brujin graph stored in this struct
   Flags
    0 - fiveprime_branching
    1 - threeprime_branching
    2 - visited
    3 - delete this
    4 -
    5 -
    6 -
    7 -
    8 -
    9 -
   10 -
   11 -
   12 -
   13 -
   14 -
   15 -
 */
typedef struct DBnode {
  uint16_t count = 0;
  uint16_t flag = 0;
  template<class T> inline int is_flag_on(T pos){ return !!(flag & (1<<pos));};
  template<class T> inline void bit_on(T pos){flag |= (1U << pos);};
  template<class T> inline void bit_off(T pos){flag &= ~(1U << pos);};
  template<class T> inline void bit_toggle(T pos){flag ^= (1U << pos);};
  template<class T1, class T2> inline void bit_set(T1 pos, T2 val){flag ^= (-(uint32_t)val ^ flag) & (1U << pos);};
} DBnode;

static inline uint64_t hash_64(uint64_t key)
{ // more sophisticated hash function to reduce collisions
  key = (~key + (key << 21)); // key = (key << 21) - key - 1;
  key = key ^ key >> 24;
  key = ((key + (key << 3)) + (key << 8)); // key * 265
  key = key ^ key >> 14;
  key = ((key + (key << 2)) + (key << 4)); // key * 21
  key = key ^ key >> 28;
  key = (key + (key << 31));
  return key;
}

void print_uint64_t(uint64_t s){
  for (uint32_t j = 63; (j >= 0) && (j < 64); j--){
    std::cout << !!(s & (1UL<<j));
  }
  std::cout << '\n';
}

KHASH_INIT(64, khint64_t, DBnode, 1, hash_64, kh_int64_hash_equal)

unsigned char seq_nt4_table[128] = { // Table to change "ACGTN" to 01234
  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
  4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
  4, 4, 4, 4,  3, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
  4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
  4, 4, 4, 4,  3, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4
};

// converts 0123 to ACGT
std::vector<std::string> int_base_to_str = { "A", "C", "G", "T"};


/** This function takes a sequence as a char array,
        and adds all the kmers into the hash table, and counts **/
template<class T>
uint64_t kmer_to_uint64(std::string str, T inp_k){
  uint32_t k = static_cast<uint32_t>(inp_k);
  if (str.size() != k){
    return 0;
  }
  int i, l;
  uint64_t x[2], mask = (1ULL<<k*2) - 1, shift = (k - 1) * 2;
  for (i = l = 0, x[0] = x[1] = 0; i < k; ++i) {
    int c = (uint8_t)str[i] < 128? seq_nt4_table[(uint8_t)str[i]] : 4;
    //printf("The char is %c. The value is %u. l is %d\n", seq[i], c, l);
    if (c < 4) { // not an "N" base
      x[0] = (x[0] << 2 | c) & mask;                  // forward strand
      x[1] = x[1] >> 2 | (uint64_t)(3 - c) << shift;  // reverse strand
      if (++l >= k) { // we find a k-mer
        return x[0] < x[1]? x[0] : x[1];
      }
    } else l = 0, x[0] = x[1] = 0; // if there is an "N", restart
  }
  return 0;
}

template<typename T>
std::string uint64_to_kmer(uint64_t kmer, T inp_k){
  uint64_t mask = (1ULL<<1)+1;
  uint64_t result = 0;
  std::string out("");
  while (inp_k-- > 0){
    result = (kmer & mask);
    switch(result){
      case(0):
        out.insert(0, 1, 'A');
        break;
      case(1):
        out.insert(0, 1, 'C');
        break;
      case(2):
        out.insert(0, 1, 'G');
        break;
      case(3):
        out.insert(0, 1, 'T');
        break;
    }
    kmer >>= 2;
  }
  return out;
}

/*
  0 - strand A
  1 - strand C
  2 - strand G
  3 - strand T
  4 - rc A
  5 - rc C
  6 - rc G
  7 - rc T
  */
template<typename T>
std::vector<uint64_t> get_extensions(uint64_t & kmer, T k){
  uint64_t shift = (k - 1) * 2;
  uint64_t mask = (1ULL<<k*2) - 1;
  uint64_t temp;
  std::vector<uint64_t> e;
  for (uint64_t i = 0; i < 4; i++){
    temp = (kmer >> 2) | (i << shift);
    e.push_back(canon(temp, k));
  }
  for (uint64_t i = 0; i < 4; i++){
    temp = ( (kmer << 2) | i) & mask;
    e.push_back( canon(temp, k));
  }
  return e;
}

/*overload for strings*/
template<typename T>
std::vector<uint64_t> get_extensions(std::string kmer, T inp_k){
  uint32_t k = static_cast<uint32_t>(inp_k);
  if (kmer.size() != k){
    std::vector<uint64_t> ret = {};
    return ret ;
  }
  int i, l;
  uint64_t x, mask = (1ULL<<k*2) - 1;
  for (i = l = 0, x = 0; i < k; ++i) {
    int c = (uint8_t)kmer[i] < 128? seq_nt4_table[(uint8_t)kmer[i]] : 4;
    //printf("The char is %c. The value is %u. l is %d\n", seq[i], c, l);
    if (c < 4) { // not an "N" base
      x = (x << 2 | c) & mask;                  // forward strand
    } else{ // if there is an N, throw an error
      std::cout << kmer << '\n';
      throw std::runtime_error("Can't get the extensions of a kmer with an N.\n");
    }
  }
  return get_extensions(x, k);
}


/* just get the rc of the seq */
template<class T>
uint64_t get_rc(uint64_t  kmer, T k){
  uint64_t rc = 0;
  uint64_t mask = 3;
  for (uint64_t i = 0; i < k*2; i+=2){
    rc = (rc << 2) | (3 - (kmer & mask));
    kmer >>= 2;
  }
  return rc;
}

template<class T>
uint64_t canon(const uint64_t & kmer, T k){
  uint64_t rc = get_rc(kmer, k);
  if (kmer < rc){
    return kmer;
  } else { return rc; }
}

class DBG {
  public:                      // begin public section
    //constructors and setup
    template<class T> DBG(T init_k);      // constructor
    ~DBG();
    // standard getters
    template<class T> DBnode* access_node(uint64_t index, T put );
    template<class T> DBnode* access_node(std::string str, T put );
    int remove_node(uint64_t index);
    // getter
    uint32_t get_k();
    uint64_t size();
    // navigate
    int mark_branching();
    int mark_all_as_unvisited();
    int delete_flagged();
    int count_nulls();
    int gen_HKCs(std::string filename);
    std::string _gen_HKC_helper(uint64_t origin, uint32_t char_cat, uint32_t extend);
    int mark_non_het_for_deletion();
    int _mark_nhfd_helper(uint64_t source, uint32_t dir);
    // just does stuff
    int print_graph();
  private:                // begin private section
    uint32_t class_k;                // member variable
    khash_t(64) *class_h; //initialise the hash table of kmers to index
};

/* constructor of HKCgraph,
    just set k and initialize the empty HKC map*/
template<class T>
DBG::DBG(T init_k) {
  class_k = static_cast<uint32_t>(init_k);
  class_h = kh_init(64);
}

/* destructor. just tear down the hash table */
DBG::~DBG() {
  kh_destroy(64, class_h);              // deallocate the hash table
}

/* just returns k*/
uint32_t DBG::get_k(){
  return class_k;
}

/* Returns the total size of the hash
   */
uint64_t DBG::size(){
  return static_cast<uint64_t>(kh_size(class_h));
}

/* If put is 1, we are adding to the graph
     - if it doesn't exist, adds it to the graph and return pointer to object
     - if it does exist, throw error since it already exists
   If put is 0, we are just accessing the node.
     - If it doesn't exist then it returns a nulptr
     - If it exists it returns the pointer to the DBnode
  */
template<class T>
DBnode* DBG::access_node(uint64_t index, T put ){
  uint64_t lookup = canon(index, class_k);
  khint_t k;
  k = kh_get(64, class_h, lookup); // query the hash table
  int is_missing, absent;
  is_missing = (int)(k == kh_end(class_h));
  if (put == 1){ //insert
    if (is_missing){
      k = kh_put(64, class_h, lookup, &absent);
      kh_val(class_h, k).count = 1;
      kh_val(class_h, k).flag = 0;
      //when we insert we also need to add edges
    } else {
      kh_val(class_h, k).count += 1;
      //don't mess with the flag here. we already instantiated it
    }
  }
  else{
    if (is_missing) return nullptr;
  }
  return &( kh_val(class_h, k) );
}

template<class T>
DBnode* DBG::access_node(std::string str, T put ){
  //std::cout << "I'm in str access node\n";
  uint64_t index = kmer_to_uint64(str, class_k);
  //std::cout << "I'm returning access node\n";
  return access_node(index, put);
}


/* This removes a kmer from the hash table
   Returns 0 if didn't remove (already gone)
   Returns 1 if removed
*/
int DBG::remove_node(uint64_t index){
  uint64_t lookup = canon(index, class_k);
  khint_t k;
  k = kh_get(64, class_h, lookup); // query the hash table
  int is_missing;
  is_missing = (int)(k == kh_end(class_h));
  if (is_missing){
    return 0;
  } else {
    kh_del(64, class_h, k);// remove a key-value pair
    return 1;
  }
  return 0;
}

/* count nulls
   */
int DBG::count_nulls(){
  khint_t k;
  uint64_t key;
  DBnode * pNode;
  int counter = 0;
  for (k = kh_begin(class_h); k != kh_end(class_h); ++k){  // traverse
    if (kh_exist(class_h, k)){            // test if a bucket contains data
      key = kh_key(class_h, k);
      pNode = access_node(key, 0);
      if (pNode == NULL){
        counter++;
      }
    }
  }
  return counter;
}

/* Delete all of the kmers that are flagged for deletion   */
int DBG::delete_flagged(){
  khint_t k;
  for (k = kh_begin(class_h); k != kh_end(class_h); ++k){  // traverse
    if (kh_exist(class_h, k)){            // test if a bucket contains data
      if (kh_val(class_h, k).is_flag_on(3) == 1){
        kh_del(64, class_h, k);// remove a key-value pair
      }
    }
  }
  return 0;
}

/* Mark all of the branching nodes
   TODO - This could be optimized if I use more flags
   */
int DBG::mark_branching(){
  khint_t k;
  uint64_t key;
  DBnode * pNode;
  uint64_t counter = 0;
  uint16_t fiveprime_counter = 0;
  uint16_t threeprime_counter = 0;
  for (k = kh_begin(class_h); k != kh_end(class_h); ++k){  // traverse
    if (kh_exist(class_h, k)){            // test if a bucket contains data
      key = kh_key(class_h, k);
      pNode = access_node(key, 0);
      std::vector<uint64_t> vec = get_extensions(key, class_k);
      //optimizes the loop a bit
      while ((counter < 4) && (fiveprime_counter < 2) ) {
        if (access_node(vec[counter], 0) != nullptr){
          fiveprime_counter++;
        }
        counter++;
      }
      while ((counter < 8) && (threeprime_counter < 2) ) {
        if (access_node(vec[counter], 0) != nullptr){
          threeprime_counter++;
        }
        counter++;
      }
      if (fiveprime_counter >= 2) pNode->bit_on(0);
      if (threeprime_counter >= 2) pNode->bit_on(1);
      fiveprime_counter = 0;
      threeprime_counter = 0;
      counter = 0;
    }
  }
  return 0;
}

/* Marks all of the nodes as unvisited
   */
int DBG::mark_all_as_unvisited(){
  khint_t k;
  uint64_t key;
  DBnode * pNode;
  for (k = kh_begin(class_h); k != kh_end(class_h); ++k){  // traverse
    if (kh_exist(class_h, k)){            // test if a bucket contains data
      key = kh_key(class_h, k);
      pNode = access_node(key, 0);
      pNode->bit_off(2);
    }
  }
  return 0;
}

/* This goes through the graph and marks the potentially homozygous regions
   for deletion.
  */
int DBG::mark_non_het_for_deletion(){
  khint_t k;
  uint64_t key;
  DBnode * pNode;
  uint16_t mask = 3;
  uint16_t fptp = 0;
  for (k = kh_begin(class_h); k != kh_end(class_h); ++k){  // traverse
    if (kh_exist(class_h, k)){            // test if a bucket contains data
      key = kh_key(class_h, k);
      pNode = access_node(key, 0);
      if ( pNode->is_flag_on(2) == 0){ //if not yet visited
        fptp = pNode->flag & mask;
        switch (fptp){
          case 3: // both fp and tp branch
            //std::cout << "both branch\n";
            pNode->bit_on(3); //just delete this node
            break;
          case 1: //just fp branch
            //search in the rp direction
            //std::cout << "search in the rp dir\n";
            _mark_nhfd_helper(key, 1);
            break;
          case 2: //just tp branch
            //search in the fp direction
            //std::cout << "search in the fp dir\n";
            _mark_nhfd_helper(key, 0);
            break;
        }
      }
      pNode->bit_on(2); // mark this one as visited
    }
  }
  return 0;
}

/*
  */
int DBG::_mark_nhfd_helper(uint64_t source,
                           uint32_t dir){
  khint_t k;
  uint64_t ext;
  DBnode * pNode;
  DBnode * pNodeT; //just a temp
  pNode = access_node(source, 0);
  pNode->bit_on(2); //mark source as visited
  pNode->bit_on(3); //mark source as delete
  uint32_t done = 0;

  // THE POINT OF THIS WHOLE BLOCK IS TO FIND EXT
  //if direction == 0 it is a 5p extension and we need to find which it is
  // dir == 0 is a 3p extension
  uint32_t start, stop;
  std::vector<uint64_t> dec = get_extensions(source, class_k);
  if (dir == 0){//0 means search in 5p direction
    start = 0;
    stop = 4;
  }
  else if ( dir == 1) {//1 means search in 3p direction
    start = 4;
    stop = 8;
  }
  else {
    throw std::runtime_error("Dir passed to DBG::_mark_nhfd_helper must be 0 or 1.\n");
  }
  uint32_t t_counter = 0;
  uint32_t ind = 8;
  for (uint32_t i = start; i < stop; i++){
    pNodeT = access_node(dec[i], 0);
    if (pNodeT != nullptr){
      t_counter++;
      ind = i;
    }
  }
  if (t_counter == 1){ //we found the extension
    k = kh_get(64, class_h, dec[ind]); // query the hash table
    ext = kh_key(class_h, k);
  }
  else{ //we found the termination of the homozygous region
    throw std::runtime_error("Trying to extend in a direction with two or more options.\n");
  }
  // END OF BLOCK TO FIND EXT

  while (done == 0){
    //std::cout << "i    Source: " << uint64_to_kmer(source, class_k) << std::endl;
    //std::cout << "i Extension: " << uint64_to_kmer(ext, class_k) << std::endl;
    //std::cout << "i    Source: ";
    //print_uint64_t(source);
    //std::cout << "i Extension: ";
    //print_uint64_t(ext);
    //std::cout << '\n';
    pNode = access_node(ext, 0);
    //first mark the node for deletion. it is no good
    pNode->bit_on(2);
    pNode->bit_on(3);
    //second determine if we stop here from branch or terminal
    //determine which direction we came from
    std::vector<uint64_t> vec = get_extensions(ext, class_k);
    uint32_t pos = 8; //only pos0-7 are legal
    for (uint32_t i = 0; i<8; i++){
      if (vec[i] == source){
        pos = i;
        break;
      }
    }
    if (pos >= 8){ //we should have found the source
      std::cout << "   Source: " << source << std::endl;
      std::cout << "Extension: " << ext << std::endl;
      std::cout << "   Source: ";
      print_uint64_t(source);
      std::cout << "Extension: ";
      print_uint64_t(ext);
      throw std::runtime_error("We didn't find the source for some reason.\n");
    }
    else if ( pos < 4 ){ //the source was 5' relative to ext. extend 3'
      //check for existence of the direction
      start = 4;
      stop = 8;
    }
    else{ // the source was 3' relative to ext. extend 5'
      start = 0;
      stop = 3;
    }

    //now figure out if there are zero, one, or more extensions
    uint32_t counter = 0;
    uint32_t new_ext_index = 8;
    for (uint32_t i = start; i < stop; i++){
      pNodeT = access_node(vec[i], 0);
      if (pNodeT != nullptr){ //not a null ptr,
        counter++;
        new_ext_index = i;
      }
    }
    if (counter == 1){ //we found the extension
      source = ext;
      k = kh_get(64, class_h, vec[new_ext_index]); // query the hash table
      ext = kh_key(class_h, k);
    }
    else{ //we found the termination of the homozygous region
      done = 1;
    }
  }
  return 0;
}

/* print out the graph */
int DBG::print_graph(){
  khint_t k, k2;
  uint64_t key, key2;
  for (k = kh_begin(class_h); k != kh_end(class_h); ++k){  // traverse
    if (kh_exist(class_h, k)){            // test if a bucket contains data
      key = kh_key(class_h, k);
      //std::cout << uint64_to_kmer(key, class_k) << "\n";
      std::vector<uint64_t> vec = get_extensions(key, class_k);
      for (const auto & n: vec){
        k2 = kh_get(64, class_h, n); // query the hash table
        if ( k2 != kh_end(class_h)){  // test if it is missing
          key2 = kh_key(class_h, k2);
          if (key < key2){
            std::cout << uint64_to_kmer(key, class_k) << "\t" << uint64_to_kmer(key2, class_k) << '\n';
          }
        }
      }
    }
  }
  return 0;
}

/* generate HKCs from the remaining reads
   */
int DBG::gen_HKCs(std::string ofilename){
  khint_t k;
  uint64_t key;
  std::ofstream myfile;
  myfile.open (ofilename);

  DBnode * pNode;
  for (k = kh_begin(class_h); k != kh_end(class_h); ++k){  // traverse
    if (kh_exist(class_h, k)){            // test if a bucket contains data
      key = kh_key(class_h, k);
      pNode = access_node(key, 0);
      if ( pNode->is_flag_on(2) == 0){
        pNode->bit_on(2);
        std::cout << "\n\n\n\n\n\n\n\n\n\nHEYYY looking at kmer 5p direction" << uint64_to_kmer(key, class_k) << "\n";
        std::string kmer = uint64_to_kmer(key, class_k);
        //extend in the 5p direction
        kmer = _gen_HKC_helper(key, 0, 0) + kmer;
        std::cout << "kmer is now " << kmer << "\n";
        //extend in the 3p direction
        std::cout << "\n now looking at kmer " << uint64_to_kmer(key, class_k) << " in the 3p direction\n";
        kmer = kmer + _gen_HKC_helper(key, 1, 1);
        myfile << kmer << std::endl;
      }
    }
  }
  myfile.close();
  return 0;
}

/* a helper function used to get the extensions in either direction
   */
std::string DBG::_gen_HKC_helper(uint64_t origin, uint32_t char_cat,
                                 uint32_t extend){
  std::cout << "inside recursion looking at kmer " << uint64_to_kmer(origin, class_k) << " extend dir: " << extend << "\n";
  uint32_t start, stop;
  DBnode * pNode;
  uint64_t ori_stem_mask, ext_stem_mask;
  //mark the current node as visited...because we visited it
  pNode = access_node(origin, 0);
  pNode->bit_on(2);
  if (extend == 0){ //if extend is 0, extend in the 5' direction of origin
    start = 0;
    stop = 4;
    ori_stem_mask = ((1ULL << ((class_k -1) * 2)) - 1) << 2;
    ext_stem_mask = (1ULL << ((class_k -1) * 2)) - 1;
  }
  else if (extend == 1){ //if extend is 1, extend in the 3' direction of origin
    start = 4;
    stop = 8;
    ori_stem_mask = (1ULL << ((class_k -1) * 2)) - 1;
    ext_stem_mask = ((1ULL << ((class_k -1) * 2)) - 1) << 2;
  }
  else {
    throw std::runtime_error("extend must be either zero or one\n");
  }
  std::vector<uint64_t> vec = get_extensions(origin, class_k);
  uint32_t ext_pos = 8;
  for (uint32_t i = start; i < stop; i++){
    std::cout << "  - Pos: " << i << " Ext: " << uint64_to_kmer(vec[i], class_k) << "\n";
    if (vec[i] != origin){ //stops the edge case of homopolymers extending into themselves
      //std::string test = "CCTCC";
      pNode = access_node(vec[i], 0);
      //if (vec[i] == kmer_to_uint64(test, class_k)){
      //  std::cout << "found CCTCC and its node exists?: " << (pNode != nullptr) << std::endl;
      //  std::cout << "found CCTCC and its visited flag on?: " << (pNode->is_flag_on(2)) << std::endl;
      //}
      if ((pNode != nullptr) && (pNode->is_flag_on(2) == 0)){//we have found the extension if it hasn't been visited
        ext_pos = i;
        std::cout << "  - found the extension " << uint64_to_kmer(vec[ext_pos], class_k) << "\n";
        break;
      }
    }
  }
  if (ext_pos == 8){ //we have not found an extension and need to exit
    std::cout << "    - found no extension\n";
    return ""; //nothing
  }
  else {
    //we now have the extension. We just need to figure out if the the extension
    // is in the same orientation as the origin, or flipped
    uint32_t same_orientation;
    // the origin and the extension are both in the same orientation
    if (extend == 0) {
      if (((origin & ori_stem_mask)>>2) == (vec[ext_pos] & ext_stem_mask)){
        same_orientation = 1;
      } else {
        same_orientation = 0;
      }
    } else {
      std::cout <<"eyy in else\n";
      if ( (origin & ori_stem_mask) == ( (vec[ext_pos] & ext_stem_mask) >> 2) ){
        same_orientation = 1;
      } else {
        same_orientation = 0;
      }
    }
    std::cout << " * origin " << uint64_to_kmer(origin, class_k) << " same ori as ext " << uint64_to_kmer(vec[ext_pos], class_k) << " ?: " << same_orientation << "\n";
    std::cout << "origin: " << origin << std::endl;
    std::cout << "   ext: " << vec[ext_pos] << std::endl;
    std::cout << "origin: ";
    print_uint64_t(origin);
    std::cout << "orimas: ";
    print_uint64_t(ori_stem_mask);
    std::cout << "   ext: ";
    print_uint64_t(vec[ext_pos]);
    std::cout << "stemas: ";
    print_uint64_t(ext_stem_mask);

    /* there are now four scenarios.
       0 - 000  extending 5p and the ext is flipped w/r/t the ori. source is 5p.
                 base is complement of 3p-most base of ext
       1 - 001  extending 5p and the ext is flipped w/r/t the ori. source is 3p.
                 base is 3p-most base of ext.
       2 - 010  extending 3p and the ext is flipped wrt the ori. source is 5p.
                 base is 5p-most base of ext.
       3 - 011   extending 3p and the ext is flipped wrt the ori. source is 3p.
                  base is complement of 5p-most base of ext.
       4 - 100  extending 5p and the ori/ext are in the same orientation. source is 5p
                 base is complement of 5p-most of the ext
       5 - 101   extending 5p and the ori/ext are in the same orientation. source is 3p
                  base is 5p-most of the ext
       6 - 110  extending 3p and the ori/ext are in the same orientation. source is 5p.
                 base is the RC of 3p-most of the ext
       7 - 111  extending 3p and the ori/ext are in the same orientation. source is 3p.
                 base is the 3p-most of the ext
     */
    uint32_t decision = (same_orientation << 2) + (extend << 1) + char_cat;
    std::string new_base;
    uint32_t run_dir;
    std::cout << "  - case : " << decision << "\n";
    switch (decision){
      case (0):
        new_base = int_base_to_str[3 - (vec[ext_pos] & 3)];
        run_dir = 1;
        break;
      case (1):
        new_base = int_base_to_str[(vec[ext_pos] & 3)];
        run_dir = 1;
        break;
      case (2):
        new_base = int_base_to_str[(vec[ext_pos] >> ((class_k -1) * 2))];
        std::cout <<"  in case 1 and new base is: " << new_base << "\n";
        run_dir = 0;
        break;
      case (3):
        new_base = int_base_to_str[3-(vec[ext_pos] >> ((class_k -1) * 2))];
        run_dir = 0;
        break;
      case (4):
        new_base = int_base_to_str[(vec[ext_pos] >> ((class_k -1) * 2))];
        run_dir = 0;
        break;
      case (5):
        new_base = int_base_to_str[3 - (vec[ext_pos] >> ((class_k -1) * 2))];
        run_dir = 0;
        break;
      case (6):
        new_base = int_base_to_str[3- (vec[ext_pos] & 3)];
        run_dir = 1;
        break;
      case (7):
        new_base = int_base_to_str[(vec[ext_pos] & 3)];
        run_dir = 1;
        break;
    }
    std::cout << "    - new base is : " << new_base << "\n";
    if (char_cat == 0){
      return _gen_HKC_helper(vec[ext_pos], char_cat, run_dir) + new_base;
    } else {
      return new_base + _gen_HKC_helper(vec[ext_pos], char_cat, run_dir) ;
    }
  }
  return "N";
}

template<class Stream, class T>
int parse_kmer_dump(DBG & G, Stream & stream, T k){
  std::string line;
  int skip = 0;
  uint64_t counter = 0;
  DBnode * pNode;
  uint64_t index;
  uint16_t this_count;
  while(std::getline(stream, line)){
    std::stringstream ss(line);
    while(std::getline(ss, line, '\t') && !skip && counter < 2){
      if (counter == 0){
        index = kmer_to_uint64(line, k);
      } else {
        uint32_t temp = std::stoul(line);
        if (temp < std::numeric_limits<uint16_t>::max()){
          this_count = static_cast<uint16_t>(temp);
        }
        else{
          this_count = std::numeric_limits<uint16_t>::max();
        }
      }
      counter++;
    }
    pNode = G.access_node(index, 1);
    pNode->count = this_count;
    counter = 0;
  }
  return 0;
}

#endif
