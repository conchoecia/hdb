#define CATCH_CONFIG_MAIN  // This tells Catch to provide a main() - only do this in one cpp file
#include <iostream>
#include <random>
#include <string>
#include <cstdio>
#include <set>
#include <typeinfo>
#include "catch.hpp"
#include "hdb.hpp"
//#include "gen_HKC.hpp"



TEST_CASE ("tests if the kmer_to_uint64(str, k) works", "[kmer_to_uint64_test]"){
  std::string str = "AAAAAAAAAAAAAAAAAAAAA";
  REQUIRE(kmer_to_uint64(str, 21) == 0);
  str = "TTTTTTTTTTTTTTTTTTTTT";
  REQUIRE(kmer_to_uint64(str, 21) == 0);

  //str and revcomp
  str = "AACTATCAGTATCACAATATC"; //000001110011010010110011010001000011001101
  REQUIRE(kmer_to_uint64(str, 21) == 123795738829);
  REQUIRE(kmer_to_uint64(str, 21) != 123795738827);
  str = "GATATTGTGATACTGATAGTT";
  REQUIRE(kmer_to_uint64(str, 21) == 123795738829);

  //str and revcomp
  str = "GTTACCGCAGTAGCAAATTGA"; //101111000101100100101100100100000011111000
  REQUIRE(kmer_to_uint64(str, 21) == 3235799777528);
  REQUIRE(kmer_to_uint64(str, 21) != 3235799777525);
  REQUIRE(kmer_to_uint64(str, 20) == 0);
  str = "TCAATTTGCTACTGCGGTAAC";
  REQUIRE(kmer_to_uint64(str, 21) == 3235799777528);

  //str and revcomp
  str = "ATAAGTATGTATTACTATTTG"; //001100001011001110110011110001110011111110
  REQUIRE(kmer_to_uint64(str, 21) == 836693335294);
  str = "CAAATAGTAATACATACTTAT";
  REQUIRE(kmer_to_uint64(str, 21) == 836693335294);

  //str and revcomp
  str = "CAGGTGACCATGTTCCATTAC"; //010010101110000101001110111101010011110001
  REQUIRE(kmer_to_uint64(str, 21) == 1286430512369);
  str = "GTAATGGAACATGGTCACCTG"; //101100001110100000010011101011010001011110
  REQUIRE(kmer_to_uint64(str, 21) == 1286430512369);

  int k = 5;
  str = "TACAC"; //1100010001 785
  //equals 748 10 11 10 11 00 which is GTGTA
  REQUIRE(kmer_to_uint64(str, k) == 748);
  str = "ATACA"; // 0011000100 196
  REQUIRE(kmer_to_uint64(str, k) == 196);
}

TEST_CASE( "confirm the node is 16 bytes", "[confirm node size]"){
  DBnode dbn;
  REQUIRE(sizeof(dbn) == 4);
}

TEST_CASE( "tests if the is_flag_on function in the DB node structs works correctly", "[test_DBnode_is_flag_on]"){
  DBnode dbn;
  dbn.flag = 1;
  for (uint16_t i = 0; i < 16; i++){
    for (uint16_t j = 15; (j >= 0) && (j < 16); j--){
      if (i==j){
        REQUIRE( dbn.is_flag_on(j) == 1);
      }
      else{
        REQUIRE( dbn.is_flag_on(j) == 0);
      }
      std::cout << dbn.is_flag_on(j);
    }
    dbn.flag <<= 1;
    std::cout << '\n';
  }
  std::cout << '\n';
}

TEST_CASE( "tests if this turns on bits correctly", "[test_DBnode_bit_on]"){
  DBnode dbn;
  //outer loop turns on flags
  for (uint16_t i = 0; i < 16; i++){
    //inner loop checks statuses
    for (uint16_t j = 15; (j >= 0) && (j < 16); j--){
      if (j >= i){ //haven't turned it on yet
        REQUIRE(dbn.is_flag_on(j) == 0);
      }
      else{
        REQUIRE(dbn.is_flag_on(j) == 1);
      }
      std::cout << dbn.is_flag_on(j);
    }
    //now turn on i for next loop
    dbn.bit_on(i);
    std::cout << '\n';
  }
  for (uint16_t j = 15; (j >= 0) && (j < 16); j--){
    REQUIRE(dbn.is_flag_on(j) == 1);
    std::cout << dbn.is_flag_on(j);
  }
  std::cout << "\n\n";
}

TEST_CASE( "tests if this turns off bits correctly", "[test_DBnode_bit_off]"){
  DBnode dbn;
  dbn.flag = std::numeric_limits<uint16_t>::max();

  //outer loop turns on flags
  for (uint16_t i = 0; i < 16; i++){
    //inner loop checks statuses
    for (uint16_t j = 15; (j >= 0) && (j < 16); j--){
      if (j >= i){ //haven't turned it on yet
        REQUIRE(dbn.is_flag_on(j) == 1);
      }
      else{
        REQUIRE(dbn.is_flag_on(j) == 0);
      }
      std::cout << dbn.is_flag_on(j);
    }
    //now turn on i for next loop
    dbn.bit_off(i);
    std::cout << '\n';
  }
  for (uint16_t j = 15; (j >= 0) && (j < 16); j--){
    REQUIRE(dbn.is_flag_on(j) == 0);
    std::cout << dbn.is_flag_on(j);
  }
  std::cout << "\n\n";
}

TEST_CASE( "tests if this toggles bits correctly", "[test_DBnode_bit_off]"){
  DBnode dbn;
  dbn.flag = std::numeric_limits<uint16_t>::max();

  //outer loop turns on flags
  for (uint16_t i = 0; i < 16; i++){
    REQUIRE(dbn.is_flag_on(i) == 1);
    dbn.bit_toggle(i);
    REQUIRE(dbn.is_flag_on(i) == 0);
    //inner loop prints
    for (uint16_t j = 15; (j >= 0) && (j < 16); j--){
      std::cout << dbn.is_flag_on(j);
    }
    std::cout << '\n';
    dbn.bit_toggle(i);
    REQUIRE(dbn.is_flag_on(i) == 1);
  }
  std::cout << '\n';
}

TEST_CASE( "tests if this sets nth bit to 0/1 correctly", "[test_DBnode_bit_set]"){
  DBnode dbn;
  dbn.flag = 0;

  //outer loop turns on flags
  for (uint16_t i = 0; i < 16; i++){
    REQUIRE(dbn.is_flag_on(i) == 0);
    dbn.bit_set(i, 1);
    REQUIRE(dbn.is_flag_on(i) == 1);
    //inner loop prints
    for (uint16_t j = 15; (j >= 0) && (j < 16); j--){
      std::cout << dbn.is_flag_on(j);
    }
    std::cout << '\n';
    dbn.bit_set(i, 0);
    REQUIRE(dbn.is_flag_on(i) == 0);
  }
  std::cout << '\n';
}

TEST_CASE( "tests if the extension generating function works", "[test_extensions]"){
  uint32_t k = 11;
  uint64_t kmer = 0;
  //just make a bunch of Gs
  for (uint32_t i = 0; i < k; i++){
    kmer = ( (kmer << 2) | 2 );
  }
  std::vector<uint64_t> answer;
  answer.push_back( 699050); //AGGGGGGGGGG
  //answer.push_back(1747626); //CGGGGGGGGGG
  answer.push_back(1398102); //CCCCCCCCCCG canonical of above
  //answer.push_back(2796202); //GGGGGGGGGGG
  answer.push_back(1398101); // CCCCCCCCCCC canonical of above
  //answer.push_back(3844778); //TGGGGGGGGGG
  answer.push_back(1398100); //CCCCCCCCCCA canonical of above
  answer.push_back(2796200); //GGGGGGGGGGA
  //answer.push_back(2796201); //GGGGGGGGGGC
  answer.push_back(2446677); //GCCCCCCCCCC canonical of above
  //answer.push_back(2796202); //GGGGGGGGGGG
  answer.push_back(1398101); // CCCCCCCCCCC canonical of above
  //answer.push_back(2796203); //GGGGGGGGGGT
  answer.push_back(349525); //ACCCCCCCCCC canonical of above

  auto results = get_extensions(kmer, k);
  for (uint16_t i = 0; i < 8; i++){
    REQUIRE(results[i] == answer[i]);
    REQUIRE(results[i] != 0);
  }


  // now do TACAC 1100010001 785
  k = 5;
  answer.clear();
  // ATACA canon ATACA 0011000100 196
  // CTACA canon CTACA 0111000100 452
  // GTACA canon GTACA 1011000100 708
  // TTACA canon TGTAA 1110110000 944
  // ACACA canon ACACA 0001000100 68
  // ACACC canon ACACC 0001000101 69
  // ACACG canon ACACG 0001000110 70
  // ACACT canon ACACT 0001000111 71
  answer.push_back(196);
  answer.push_back(452);
  answer.push_back(708);
  answer.push_back(944);
  answer.push_back(68);
  answer.push_back(69);
  answer.push_back(70);
  answer.push_back(71);
  kmer = 785;
  results = get_extensions(kmer, k);
  for (uint16_t i = 0; i < 8; i++){
    //std::cout << "i: " << i << " results[i]: " << uint64_to_kmer(results[i],5) << " answer[i]: " << uint64_to_kmer(answer[i],5) << '\n';
    REQUIRE(results[i] == answer[i]);
    REQUIRE(results[i] != 0);
  }

  std::string tester = "tacac"; //785
  results = get_extensions(tester, k);
  for (uint16_t i = 0; i < 8; i++){
    std::cout << "i: " << i << " results[i]: " << uint64_to_kmer(results[i],5) << " answer[i]: " << uint64_to_kmer(answer[i],5) << '\n';
    REQUIRE(results[i] == answer[i]);
    REQUIRE(results[i] != 0);
  }

  //now do AAGGT canon AAGGT 0000101011 43
  // AAAGG canon AAAGG 0000001010  10
  // CAAGG canon CAAGG 0100001010 266
  // GAAGG canon CCTTC 0101111101 381
  // TAAGG canon CCTTA 0101111100 380
  // AGGTA canon AGGTA 0010101100 172
  // AGGTC canon AGGTC 0010101101 173
  // AGGTG canon AGGTG 0010101110 174
  // AGGTT canon AACCT 0000010111  23
  answer.clear();
  answer.push_back(  10 );
  answer.push_back( 266 );
  answer.push_back( 381 );
  answer.push_back( 380 );
  answer.push_back( 172 );
  answer.push_back( 173 );
  answer.push_back( 174 );
  answer.push_back(  23 );
  kmer = 43;
  results = get_extensions(kmer, k);
  for (uint16_t i = 0; i < 8; i++){
    //std::cout << "i: " << i << " results[i]: " << uint64_to_kmer(results[i],5) << " answer[i]: " << uint64_to_kmer(answer[i],5) << '\n';
    REQUIRE(results[i] == answer[i]);
    REQUIRE(results[i] != 0);
  }
  tester = "AAGGT"; //785
  results = get_extensions(tester, k);
  for (uint16_t i = 0; i < 8; i++){
    std::cout << "i: " << i << " results[i]: " << uint64_to_kmer(results[i],5) << " answer[i]: " << uint64_to_kmer(answer[i],5) << '\n';
    REQUIRE(results[i] == answer[i]);
    REQUIRE(results[i] != 0);
  }


}

TEST_CASE( "tests if the rc function works properly", "[test_rc_converter]"){
  uint32_t k = 6;
  uint64_t kmer = 0;

  uint64_t rc = get_rc(kmer, k);
  REQUIRE(get_rc(kmer,k) == 4095);

  //this makes the sequence 'gacgtg' -- the rc is cacgtc`
  kmer = 2;
  for (uint64_t i = 0; i < 4; i++){
    kmer = ( (kmer << 2) | i );
  }
  kmer = ( (kmer << 2) | 2 );

  rc = get_rc(kmer, k);
  REQUIRE(rc == 1133);
}

TEST_CASE( "tests if the k return works", "[test_k_return]"){
  DBG G = DBG(29);
  REQUIRE(G.get_k() == 29);

  DBG H = DBG(1);
  REQUIRE(H.get_k() == 1);
}

/*Takes a sequence that is max 32 bp long and converts it to its 2-bit
  encoded canonical form.
  */
template<class T1, class T2>
uint64_t char_array_to_uint64_t(std::string seq, T1 len, T2 k){
  int i, l;
  uint64_t x[2], mask = (1ULL<<k*2) - 1, shift = (k - 1) * 2;
  for (i = l = 0, x[0] = x[1] = 0; i < len; ++i) {
    uint32_t c = (uint8_t)seq[i] < 128? seq_nt4_table[(uint8_t)seq[i]] : 4;
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


TEST_CASE( "tests if the hash table works as intended at all", "[test_hash_table]"){
  khint_t k;
  khash_t(64) *h = kh_init(64);  // allocate a hash table
  uint32_t counter = 0;
  for (k = kh_begin(h); k != kh_end(h); ++k){  // traverse
    if (kh_exist(h, k)) {            // test if a bucket contains data
      std::cout << "counter " << counter << ": " << kh_val(h, k).count << '\n';
    }
  } std::cout << '\n';
  REQUIRE(counter == 0); // there shouldn't be anything in the beginning

  int absent;
  uint64_t index = 0;
  k = kh_put(64, h, index, &absent);  // insert a key to the hash table
  kh_value(h, k).count = 10;             // set the value
  REQUIRE((uint32_t)kh_value(h,k).count == 10);
  k = kh_get(64, h, 10);          // query the hash table
  REQUIRE(1 == (int)(k == kh_end(h))); // 10 as a key should be missing

  counter = 0;
  std::cout << "Should be something here.\n";
  for (k = kh_begin(h); k != kh_end(h); ++k){  // traverse
    if (kh_exist(h, k)) {            // test if a bucket contains data
      std::cout << "counter " << counter++ << ": " << kh_val(h, k).count << '\n';
    }
  } std::cout << '\n';
  REQUIRE(counter ==1);

  index = 7;
  k = kh_put(64, h, index, &absent);  // insert a key to the hash table
  kh_value(h, k).count = 5;             // set the value
  REQUIRE((uint32_t)kh_value(h,k).count == 5);
  counter = 0;
  for (k = kh_begin(h); k != kh_end(h); ++k){  // traverse
    if (kh_exist(h, k)) {            // test if a bucket contains data
      std::cout << "counter " << counter++ << ": " << kh_val(h, k).count << '\n';
    }
  } std::cout << '\n';
  REQUIRE(counter == 2);

  k = kh_get(64, h, 0);          // query the hash table
  kh_del(64, h, k);               // remove a key-value pair
  k = kh_get(64, h, 7);          // query the hash table
  kh_del(64, h, k);               // remove a key-value pair
  counter = 0;
  for (k = kh_begin(h); k != kh_end(h); ++k){  // traverse
    if (kh_exist(h, k)) {            // test if a bucket contains data
      std::cout << "counter " << counter++ << ": " << kh_val(h, k).count << '\n';
    }
  } std::cout << '\n';
  REQUIRE(counter ==0);
  kh_destroy(64, h);              // deallocate the hash table
}

/* We need to test if we are able to correctly get the key from the iterator
   of the hash. If so, we don't have to worry about storing the key in the
   data structure, thereby reducing the memory footprint.
  */
TEST_CASE("tests_get_key", "[test_hash_table_get_key_from_iterator]"){
  khint_t k;
  khash_t(64) *h = kh_init(64);  // allocate a hash table
  int absent = 0;
  //make sure the keys are the same while adding them
  for (uint64_t index = 0; index < 100; index ++){
    k = kh_put(64, h, index, &absent);  // insert a key to the hash table
    REQUIRE(kh_key(h,k) == index);
  }
  //NOTE iterating through the keys does not guarantee that the keys
  // will be lexicographically sorted
  //uint64_t counter = 0;
  //for (k = kh_begin(h); k != kh_end(h); ++k){  // traverse
  //  REQUIRE(kh_key(h,k) == counter++);
  //}
}

TEST_CASE( "tests if the hash table access works. also tests if the DBG.size() func works", "[test_hash_table_class_and_size]"){
  DBG G = DBG(1);
  uint64_t index = 1;
  uint64_t rc = get_rc(index, 1);
  std::cout << "*AN* calculating lookup\n";
  uint64_t lookup = index < rc ? index : rc;
  DBnode* pNode = G.access_node(index, 0) ;
  std::cout << "Accessed the node OK\n";
  REQUIRE(pNode == nullptr);
  REQUIRE(G.size() == 0);
  pNode = G.access_node(index, 1) ;
  REQUIRE(pNode != nullptr);
  REQUIRE(pNode->count == 1);
  REQUIRE(G.size() == 1);

  //now add a second entry
  index+=2;
  rc = get_rc(index, 1);
  lookup = index < rc ? index : rc;
  pNode = G.access_node(index, 1) ;
  REQUIRE(pNode != nullptr);
  REQUIRE(pNode->count == 1);
  REQUIRE(G.size() == 2);
  pNode->count++;
  REQUIRE(pNode->count == 2);
  REQUIRE(G.size() == 2);
}

TEST_CASE("does a thorough test of adding all kmers to the hash table", "[thorough_add]"){
  //generate 10 tests for each odd kmer from 5 up to 32
  for (uint32_t i = 1; i<=5; i+=2){
    //std::cout << "i: " << i << '\n';
    uint32_t counter = 0;
    DBG G = DBG(i);
    std::set<uint64_t> tried;
    DBnode* pNode;
    REQUIRE(G.size() == 0);

    uint64_t max_for_k = 1;
    for (uint64_t m = 0; m<(i*2)-1; m++){
      max_for_k = (max_for_k << 1) | 1ULL;
    }
    for (uint64_t index = 0; index <= max_for_k; index++){
      //std::cout << "INDEX IS " << index << '\n';
      uint64_t rc = get_rc(index, i);
      if ((tried.count(index) == 0) && (tried.count(rc) == 0)){
        //std::cout << "inserting into tried\n";
        tried.insert(index);
        tried.insert(rc);
        //make sure that the node doesn't yet exist
        //std::cout << "requiring nonexistence\n";
        pNode = G.access_node(index, 0);
        REQUIRE(pNode == nullptr);
        //add the node and make sure that its count is 1;
        //std::cout << "adding a node\n";
        pNode = G.access_node(index, 1);
        //std::cout << "checking the node's count\n";
        REQUIRE(pNode->count == 1);
        ////add another value and make sure the count is 2
        //std::cout << "accessing the node again\n";
        pNode = G.access_node(index, 1);
        //std::cout << "checking the node's count again\n";
        REQUIRE(pNode->count == 2);
        ////add another value and make sure the count is 3
        //std::cout << "accessing the node again for the third time\n";
        pNode = G.access_node(index, 1);
        //std::cout << "checking the node's count for the third time\n";
        REQUIRE(pNode->count == 3);
        //std::cout << "incrementing counter\n";
        REQUIRE(G.size() == counter+1);
        counter++;
      }
      else { // just skip this
        ;
        //std::cout << " - already in tried\n";
      }
      //std::cout <<'\n';
    }
    std::cout << "k=" << i << ": added " << counter << " sequences\n";
  }
}

TEST_CASE("tests the canonical function", "[makes sure that the canonical kmer checker is correct]"){
  //87381  = 010101010101010101
  //174762 = 101010101010101010
  uint64_t index = 87381;
  REQUIRE(canon(index, 9) == 87381);
  index = 174762;
  REQUIRE(canon(index, 9) == 87381);

  //2796202 = 1010101010101010101010
  //1398101 = 0101010101010101010101
  index = 2796202;
  REQUIRE(canon(index, 11) == 1398101);
  index = 1398101;
  REQUIRE(canon(index, 11) == 1398101);
}

TEST_CASE( "tests if the remove_node function works", "[test_remove_node]"){
  DBG G = DBG(5);

  //std::cout << "*AN* calculating lookup\n";
  DBnode* pNode;

  //add everything
  for (uint64_t index = 0; index <100; index++){
    if (index == canon(index, 5)){
      pNode = G.access_node(index, 1);
    }
  }

  //go back and delete it all.
  int ret;
  for (uint64_t index = 0; index <100; index++){
    if (index == canon(index, 5)){
      pNode = G.access_node(index, 0);
      REQUIRE(pNode != nullptr);
      ret = G.remove_node(index);
      pNode = G.access_node(index, 0);
      REQUIRE(pNode == nullptr);
      REQUIRE(ret == 1);
      ret = G.remove_node(index);
      REQUIRE(ret == 0);
    }
  }
}


TEST_CASE( "tests if the parser for kmer dump input works", "[test_kmer_dump_parser]"){
  std::string ofilename = "test_dump.tsv";
  std::ofstream myfile;
  myfile.open (ofilename);
  myfile << "AAAAAAAAAAAAAAAAAAAAA\t133274\n";
  myfile << "AACTATCAGTATCACAATATC\t64\n";
  myfile << "GTTACCGCAGTAGCAAATTGA\t50\n";
  myfile << "ATAAGTATGTATTACTATTTG\t82\n";
  myfile << "CAGGTGACCATGTTCCATTAC\t51\n";
  myfile << "GGAATATTAATTCAAAAGATA\t57\n";
  myfile << "CTTTATAGTTTTGGGACCCAA\t58\n";
  myfile << "AAATTCTTTTACCACTTATTC\t79\n";
  myfile << "CGTAAATTTCCTATCGATATA\t51\n";
  myfile << "CAATGAAATCTCTGTTGGCGG\t118\n";
  myfile << "AAAAGTACATGTATTGATTTT\t105\n";
  myfile << "CATGCCAAAGTCAAATTTTGA\t42\n";
  myfile << "AATAAACTAGGGATTCCAACG\t45\n";
  myfile << "AGACACATTTAGAGAGAGGGA\t20\n";
  myfile.close();

  DBG G = DBG(21);
  std::ifstream file(ofilename.c_str());
  parse_kmer_dump(G, file, G.get_k());

  //make sure all the kmers are in the graph
  std::string str = "AAAAAAAAAAAAAAAAAAAAA\t133274\n";
  DBnode * pNode = G.access_node(str, 0);
  REQUIRE(pNode != nullptr);
  REQUIRE(pNode->count != 133274);
  REQUIRE(pNode->count == 65535);
  pNode = nullptr;

  str = "AACTATCAGTATCACAATATC";
  pNode = G.access_node(str, 0);
  REQUIRE(pNode != nullptr);
  REQUIRE(pNode->count == 64);
  pNode = nullptr;

  str = "GTTACCGCAGTAGCAAATTGA";
  pNode = G.access_node(str, 0);
  REQUIRE(pNode != nullptr);
  REQUIRE(pNode->count == 50);
  pNode = nullptr;

  str = "ATAAGTATGTATTACTATTTG";
  pNode = G.access_node(str, 0);
  REQUIRE(pNode != nullptr);
  REQUIRE(pNode->count == 82);
  pNode = nullptr;

  str = "CAGGTGACCATGTTCCATTAC";
  pNode = G.access_node(str, 0);
  REQUIRE(pNode != nullptr);
  REQUIRE(pNode->count == 51);
  pNode = nullptr;

  str = "GGAATATTAATTCAAAAGATA";
  pNode = G.access_node(str, 0);
  REQUIRE(pNode != nullptr);
  REQUIRE(pNode->count == 57);
  pNode = nullptr;

  str = "CTTTATAGTTTTGGGACCCAA";
  pNode = G.access_node(str, 0);
  REQUIRE(pNode != nullptr);
  REQUIRE(pNode->count == 58);
  pNode = nullptr;

  str = "AAATTCTTTTACCACTTATTC";
  pNode = G.access_node(str, 0);
  REQUIRE(pNode != nullptr);
  REQUIRE(pNode->count == 79);
  pNode = nullptr;

  str = "CGTAAATTTCCTATCGATATA";
  pNode = G.access_node(str, 0);
  REQUIRE(pNode != nullptr);
  REQUIRE(pNode->count == 51);
  pNode = nullptr;

  str = "CAATGAAATCTCTGTTGGCGG";
  pNode = G.access_node(str, 0);
  REQUIRE(pNode != nullptr);
  REQUIRE(pNode->count == 118);
  pNode = nullptr;

  str = "AAAAGTACATGTATTGATTTT";
  pNode = G.access_node(str, 0);
  REQUIRE(pNode != nullptr);
  REQUIRE(pNode->count == 105);
  pNode = nullptr;

  str = "CATGCCAAAGTCAAATTTTGA";
  pNode = G.access_node(str, 0);
  REQUIRE(pNode != nullptr);
  REQUIRE(pNode->count == 42);
  pNode = nullptr;

  str = "AATAAACTAGGGATTCCAACG";
  pNode = G.access_node(str, 0);
  REQUIRE(pNode != nullptr);
  REQUIRE(pNode->count == 45);
  pNode = nullptr;

  str = "AGACACATTTAGAGAGAGGGA";
  pNode = G.access_node(str, 0);
  REQUIRE(pNode != nullptr);
  REQUIRE(pNode->count == 20);
  pNode = nullptr;
}

TEST_CASE( "tests if the uint64_to_str works", "[test_uint64_to_str_works]"){
  std::string str;
  uint64_t result;
  result = 123795738829;
  str = "AACTATCAGTATCACAATATC"; //000001110011010010110011010001000011001101
  REQUIRE(uint64_to_kmer(result, 21) == str);

  //str and revcomp
  result = 3235799777528;
  str = "GTTACCGCAGTAGCAAATTGA"; //101111000101100100101100100100000011111000
  REQUIRE(uint64_to_kmer(result, 21) == str);

  //str and revcomp
  result = 836693335294;
  str = "ATAAGTATGTATTACTATTTG"; //001100001011001110110011110001110011111110
  REQUIRE(uint64_to_kmer(result, 21) == str);

  //str and revcomp
  result = 1286430512369;
  str = "CAGGTGACCATGTTCCATTAC"; //010010101110000101001110111101010011110001
  REQUIRE(uint64_to_kmer(result, 21) == str);
}


TEST_CASE( "tests if the mark_branching function works on a double branch",
           "[test_kmer_dump_parser_double]"){
  /*make sure all the kmers are in the graph

      This De Bruijn Graph is shaped like this
       \          /
        \________/
        /        \
       /          \
  */
  DBG G = DBG(5);
  std::string str = "GACAC";
  DBnode * pNode = G.access_node(str, 1);
  REQUIRE(pNode != nullptr);
  REQUIRE(pNode->is_flag_on(0) == 0);
  REQUIRE(pNode->is_flag_on(1) == 0);

  str = "TACAC";
  pNode = G.access_node(str, 1);
  REQUIRE(pNode != nullptr);
  REQUIRE(pNode->is_flag_on(0) == 0);
  REQUIRE(pNode->is_flag_on(1) == 0);

  str =  "ACACC";
  pNode = G.access_node(str, 1);
  REQUIRE(pNode != nullptr);
  REQUIRE(pNode->is_flag_on(0) == 0);
  REQUIRE(pNode->is_flag_on(1) == 0);

  str =   "CACCC";
  pNode = G.access_node(str, 1);
  REQUIRE(pNode != nullptr);
  REQUIRE(pNode->is_flag_on(0) == 0);
  REQUIRE(pNode->is_flag_on(1) == 0);

  str =   "CACCG";
  pNode = G.access_node(str, 1);
  REQUIRE(pNode != nullptr);
  REQUIRE(pNode->is_flag_on(0) == 0);
  REQUIRE(pNode->is_flag_on(1) == 0);

  G.mark_branching();

  str = "GACAC";
  pNode = G.access_node(str, 0);
  REQUIRE(pNode != nullptr);
  REQUIRE(pNode->is_flag_on(0) == 0);
  REQUIRE(pNode->is_flag_on(1) == 0);

  str = "TACAC";
  pNode = G.access_node(str, 0);
  REQUIRE(pNode != nullptr);
  REQUIRE(pNode->is_flag_on(0) == 0);
  REQUIRE(pNode->is_flag_on(1) == 0);

  str =  "ACACC";
  pNode = G.access_node(str, 0);
  REQUIRE(pNode != nullptr);
  REQUIRE(pNode->is_flag_on(0) == 1);
  REQUIRE(pNode->is_flag_on(1) == 1);

  str =   "CACCC";
  pNode = G.access_node(str, 0);
  REQUIRE(pNode != nullptr);
  REQUIRE(pNode->is_flag_on(0) == 0);
  REQUIRE(pNode->is_flag_on(1) == 0);

  str =   "CACCG";
  pNode = G.access_node(str, 0);
  REQUIRE(pNode != nullptr);
  REQUIRE(pNode->is_flag_on(0) == 0);
  REQUIRE(pNode->is_flag_on(1) == 0);
}

TEST_CASE( "tests if the mark_branching function works on multiple branches", \
           "[test_kmer_dump_parser_multiple]"){
  /*make sure all the kmers are in the graph

      This De Bruijn Graph is shaped like this
       \          /      /
        \________/______/
        /               \
       /                 \
  */
  DBG G = DBG(5);
  std::string str = "GACAC";
  DBnode * pNode = G.access_node(str, 1);
  REQUIRE(pNode != nullptr);
  REQUIRE(pNode->is_flag_on(0) == 0);
  REQUIRE(pNode->is_flag_on(1) == 0);

  str = "TACAC";
  pNode = G.access_node(str, 1);
  REQUIRE(pNode != nullptr);
  REQUIRE(pNode->is_flag_on(0) == 0);
  REQUIRE(pNode->is_flag_on(1) == 0);

  str =  "ACACC";
  pNode = G.access_node(str, 1);
  REQUIRE(pNode != nullptr);
  REQUIRE(pNode->is_flag_on(0) == 0);
  REQUIRE(pNode->is_flag_on(1) == 0);

  str =   "CACCC";
  pNode = G.access_node(str, 1);
  REQUIRE(pNode != nullptr);
  REQUIRE(pNode->is_flag_on(0) == 0);
  REQUIRE(pNode->is_flag_on(1) == 0);

  str =   "CACCG";
  pNode = G.access_node(str, 1);
  REQUIRE(pNode != nullptr);
  REQUIRE(pNode->is_flag_on(0) == 0);
  REQUIRE(pNode->is_flag_on(1) == 0);

  str =   "ACCGA";
  pNode = G.access_node(str, 1);
  REQUIRE(pNode != nullptr);
  REQUIRE(pNode->is_flag_on(0) == 0);
  REQUIRE(pNode->is_flag_on(1) == 0);

  str =   "ACCGT";
  pNode = G.access_node(str, 1);
  REQUIRE(pNode != nullptr);
  REQUIRE(pNode->is_flag_on(0) == 0);
  REQUIRE(pNode->is_flag_on(1) == 0);

  G.mark_branching();

  str = "GACAC";
  pNode = G.access_node(str, 0);
  REQUIRE(pNode != nullptr);
  REQUIRE(pNode->is_flag_on(0) == 0);
  REQUIRE(pNode->is_flag_on(1) == 0);

  str = "TACAC";
  pNode = G.access_node(str, 0);
  REQUIRE(pNode != nullptr);
  REQUIRE(pNode->is_flag_on(0) == 0);
  REQUIRE(pNode->is_flag_on(1) == 0);

  str =  "ACACC";
  pNode = G.access_node(str, 0);
  REQUIRE(pNode != nullptr);
  REQUIRE(pNode->is_flag_on(0) == 1);
  REQUIRE(pNode->is_flag_on(1) == 1);

  str =   "CACCC";
  pNode = G.access_node(str, 0);
  REQUIRE(pNode != nullptr);
  REQUIRE(pNode->is_flag_on(0) == 0);
  REQUIRE(pNode->is_flag_on(1) == 0);

  str =   "CACCG";
  pNode = G.access_node(str, 0);
  REQUIRE(pNode != nullptr);
  REQUIRE(pNode->is_flag_on(0) == 0);
  REQUIRE(pNode->is_flag_on(1) == 1);

  str =   "ACCGA";
  pNode = G.access_node(str, 0);
  REQUIRE(pNode != nullptr);
  REQUIRE(pNode->is_flag_on(0) == 0);
  REQUIRE(pNode->is_flag_on(1) == 0);

  str =   "ACCGT";
  pNode = G.access_node(str, 0);
  REQUIRE(pNode != nullptr);
  REQUIRE(pNode->is_flag_on(0) == 0);
  REQUIRE(pNode->is_flag_on(1) == 0);
}


TEST_CASE("tests if the mark_branching works on right branches",\
           "[test_kmer_dump_parser_multiple_right]"){
  /*make sure all the kmers are in the graph

      This De Bruijn Graph is shaped like this
                /      /
       ________/______/
                      \
                       \
  */

  DBG H = DBG(5);
  std::string str;
  DBnode * pNode;
  str =  "ACACC";
  pNode = H.access_node(str, 1);
  REQUIRE(pNode != nullptr);
  REQUIRE(pNode->is_flag_on(0) == 0);
  REQUIRE(pNode->is_flag_on(1) == 0);

  str =   "CACCC";
  pNode = H.access_node(str, 1);
  REQUIRE(pNode != nullptr);
  REQUIRE(pNode->is_flag_on(0) == 0);
  REQUIRE(pNode->is_flag_on(1) == 0);

  str =   "CACCG";
  pNode = H.access_node(str, 1);
  REQUIRE(pNode != nullptr);
  REQUIRE(pNode->is_flag_on(0) == 0);
  REQUIRE(pNode->is_flag_on(1) == 0);

  str =   "ACCGA";
  pNode = H.access_node(str, 1);
  REQUIRE(pNode != nullptr);
  REQUIRE(pNode->is_flag_on(0) == 0);
  REQUIRE(pNode->is_flag_on(1) == 0);

  str =   "ACCGT";
  pNode = H.access_node(str, 1);
  REQUIRE(pNode != nullptr);
  REQUIRE(pNode->is_flag_on(0) == 0);
  REQUIRE(pNode->is_flag_on(1) == 0);

  H.mark_branching();
  str =  "ACACC";
  pNode = H.access_node(str, 0);
  REQUIRE(pNode != nullptr);
  REQUIRE(pNode->is_flag_on(0) == 0);
  REQUIRE(pNode->is_flag_on(1) == 1);

  str =   "CACCC";
  pNode = H.access_node(str, 0);
  REQUIRE(pNode != nullptr);
  REQUIRE(pNode->is_flag_on(0) == 0);
  REQUIRE(pNode->is_flag_on(1) == 0);

  str =   "CACCG";
  pNode = H.access_node(str, 0);
  REQUIRE(pNode != nullptr);
  REQUIRE(pNode->is_flag_on(0) == 0);
  REQUIRE(pNode->is_flag_on(1) == 1);

  str =   "ACCGA";
  pNode = H.access_node(str, 0);
  REQUIRE(pNode != nullptr);
  REQUIRE(pNode->is_flag_on(0) == 0);
  REQUIRE(pNode->is_flag_on(1) == 0);

  str =   "ACCGT";
  pNode = H.access_node(str, 0);
  REQUIRE(pNode != nullptr);
  REQUIRE(pNode->is_flag_on(0) == 0);
  REQUIRE(pNode->is_flag_on(1) == 0);
}

TEST_CASE( "tests if the mark_branching function works", "[test_kmer_dump_parser]"){
  DBG G = DBG(5);
  //make sure all the kmers are in the graph
  //
  //  \  This De Bruijn Graph is shaped like this
  //   \________
  //   /
  //  /
  std::string str;
  DBnode * pNode;
  str = "GACAC";
  pNode = G.access_node(str, 1);
  REQUIRE(pNode != nullptr);
  REQUIRE(pNode->is_flag_on(0) == 0);
  REQUIRE(pNode->is_flag_on(1) == 0);
  str = "TACAC";
  pNode = G.access_node(str, 1);
  REQUIRE(pNode != nullptr);
  REQUIRE(pNode->is_flag_on(0) == 0);
  REQUIRE(pNode->is_flag_on(1) == 0);
  str =  "ACACC";
  pNode = G.access_node(str, 1);
  REQUIRE(pNode != nullptr);
  REQUIRE(pNode->is_flag_on(0) == 0);
  REQUIRE(pNode->is_flag_on(1) == 0);
  str =   "CACCC";
  pNode = G.access_node(str, 1);
  REQUIRE(pNode != nullptr);
  REQUIRE(pNode->is_flag_on(0) == 0);
  REQUIRE(pNode->is_flag_on(1) == 0);
  G.mark_branching();
  str = "GACAC";
  pNode = G.access_node(str, 0);
  REQUIRE(pNode != nullptr);
  REQUIRE(pNode->is_flag_on(0) == 0);
  REQUIRE(pNode->is_flag_on(1) == 0);
  str = "TACAC";
  pNode = G.access_node(str, 0);
  REQUIRE(pNode != nullptr);
  REQUIRE(pNode->is_flag_on(0) == 0);
  REQUIRE(pNode->is_flag_on(1) == 0);
  str =  "ACACC";
  pNode = G.access_node(str, 0);
  REQUIRE(pNode != nullptr);
  REQUIRE(pNode->is_flag_on(0) == 1);
  REQUIRE(pNode->is_flag_on(1) == 0);
  str =   "CACCC";
  pNode = G.access_node(str, 0);
  REQUIRE(pNode != nullptr);
  REQUIRE(pNode->is_flag_on(0) == 0);
  REQUIRE(pNode->is_flag_on(1) == 0);
}

TEST_CASE( "another case of mark_branching. has some non-canonical input.", "[mark_branching_again]"){
  DBG G = DBG(5);
  DBnode * pNode;
  std::string str;

  str = "CGGGG"; // canon CCCCG - should not branch for either
  pNode = G.access_node(str, 1);
  REQUIRE(pNode != nullptr);
  REQUIRE(pNode->is_flag_on(0) == 0);
  REQUIRE(pNode->is_flag_on(1) == 0);

  str = "TGGGG"; // canon CCCCA - should not branch for either
  pNode = G.access_node(str, 1);
  REQUIRE(pNode != nullptr);
  REQUIRE(pNode->is_flag_on(0) == 0);
  REQUIRE(pNode->is_flag_on(1) == 0);

  str =  "GGGGA"; // canon GGGGA should branch 5' not 3'
  pNode = G.access_node(str, 1);
  REQUIRE(pNode != nullptr);
  REQUIRE(pNode->is_flag_on(0) == 0);
  REQUIRE(pNode->is_flag_on(1) == 0);

  str =   "GGGAA"; // GGGAA should not branch
  pNode = G.access_node(str, 1);
  REQUIRE(pNode != nullptr);
  REQUIRE(pNode->is_flag_on(0) == 0);
  REQUIRE(pNode->is_flag_on(1) == 0);

  str =   "GGAAA"; // GGAAA should  branch 3' not 5'
  pNode = G.access_node(str, 1);
  REQUIRE(pNode != nullptr);
  REQUIRE(pNode->is_flag_on(0) == 0);
  REQUIRE(pNode->is_flag_on(1) == 0);

  str =   "GAAAA"; // GAAAA should not branch
  pNode = G.access_node(str, 1);
  REQUIRE(pNode != nullptr);
  REQUIRE(pNode->is_flag_on(0) == 0);
  REQUIRE(pNode->is_flag_on(1) == 0);

  str =   "GAAAT"; // ATTTC should not branch
  pNode = G.access_node(str, 1);
  REQUIRE(pNode != nullptr);
  REQUIRE(pNode->is_flag_on(0) == 0);
  REQUIRE(pNode->is_flag_on(1) == 0);

  G.mark_branching();

  str = "CGGGG"; // canon CCCCG - should not branch for either
  pNode = G.access_node(str, 0);
  REQUIRE(pNode != nullptr);
  REQUIRE(pNode->is_flag_on(0) == 0);
  REQUIRE(pNode->is_flag_on(1) == 0);

  str = "TGGGG"; // canon CCCCA - should not branch for either
  pNode = G.access_node(str, 0);
  REQUIRE(pNode != nullptr);
  REQUIRE(pNode->is_flag_on(0) == 0);
  REQUIRE(pNode->is_flag_on(1) == 0);

  str =  "GGGGA"; // canon GGGGA should branch 5' not 3'
  pNode = G.access_node(str, 0);
  REQUIRE(pNode != nullptr);
  REQUIRE(pNode->is_flag_on(0) == 1);
  REQUIRE(pNode->is_flag_on(1) == 0);

  str =   "GGGAA"; // GGGAA should not branch
  pNode = G.access_node(str, 0);
  REQUIRE(pNode != nullptr);
  REQUIRE(pNode->is_flag_on(0) == 0);
  REQUIRE(pNode->is_flag_on(1) == 0);

  str =   "GGAAA"; // GGAAA should  branch 3' not 5'
  pNode = G.access_node(str, 0);
  REQUIRE(pNode != nullptr);
  REQUIRE(pNode->is_flag_on(0) == 0);
  REQUIRE(pNode->is_flag_on(1) == 1);

  str =   "GAAAA"; // GAAAA should not branch
  pNode = G.access_node(str, 0);
  REQUIRE(pNode != nullptr);
  REQUIRE(pNode->is_flag_on(0) == 0);
  REQUIRE(pNode->is_flag_on(1) == 0);

  str =   "GAAAT"; // ATTTC should not branch
  pNode = G.access_node(str, 0);
  REQUIRE(pNode->is_flag_on(0) == 0);
  REQUIRE(pNode->is_flag_on(1) == 0);

}

TEST_CASE( "tests if the mark_all_as_unvisited function works",
           "[test_mark_all_as_unvisited]"){
  DBG G = DBG(5);
  std::string str;
  DBnode * pNode;
  str = "GACAC";
  pNode = G.access_node(str, 1);
  REQUIRE(pNode != nullptr);
  REQUIRE(pNode->is_flag_on(2) == 0);
  pNode->bit_on(2);
  REQUIRE(pNode->is_flag_on(2) == 1);
  str = "TACAC";
  pNode = G.access_node(str, 1);
  REQUIRE(pNode != nullptr);
  REQUIRE(pNode->is_flag_on(2) == 0);
  pNode->bit_on(2);
  REQUIRE(pNode->is_flag_on(2) == 1);
  str =  "ACACC";
  pNode = G.access_node(str, 1);
  REQUIRE(pNode != nullptr);
  REQUIRE(pNode->is_flag_on(2) == 0);
  pNode->bit_on(2);
  REQUIRE(pNode->is_flag_on(2) == 1);
  str =   "CACCC";
  pNode = G.access_node(str, 1);
  REQUIRE(pNode != nullptr);
  REQUIRE(pNode->is_flag_on(2) == 0);
  pNode->bit_on(2);
  REQUIRE(pNode->is_flag_on(2) == 1);
  //now mark all of the visited flags as off
  G.mark_all_as_unvisited();
  str = "GACAC";
  pNode = G.access_node(str, 0);
  REQUIRE(pNode != nullptr);
  REQUIRE(pNode->is_flag_on(2) == 0);
  str = "TACAC";
  pNode = G.access_node(str, 0);
  REQUIRE(pNode != nullptr);
  REQUIRE(pNode->is_flag_on(2) == 0);
  str =  "ACACC";
  pNode = G.access_node(str, 0);
  REQUIRE(pNode != nullptr);
  REQUIRE(pNode->is_flag_on(2) == 0);
  str =   "CACCC";
  pNode = G.access_node(str, 0);
  REQUIRE(pNode != nullptr);
  REQUIRE(pNode->is_flag_on(2) == 0);
}


/* Tests if source is a 5' extension relative to ext
  */
TEST_CASE( "tests how to match kmer stems",
           "[test_fivep_ext]"){
  uint32_t k = 23;
  std::string str;
  str = "AAATACCGCAGTAGCAAATTGGG";
  uint64_t source = kmer_to_uint64(str, k);
  std::cout << "source: ";

  print_uint64_t(source);
  //3' extension
  str =  "AATACCGCAGTAGCAAATTGGGG";
  uint64_t ext = kmer_to_uint64(str, k);
  std::cout << "   ext: ";
  print_uint64_t(ext);

  uint64_t sMask = (1ULL <<(k-1)*2)-1; // for looking at the kmer stems
  std::cout << "  mask: ";
  print_uint64_t(sMask);

  uint64_t temp = (ext >> 2) & sMask;
  std::cout << "  temp: ";
  print_uint64_t(temp);

  uint64_t temp2 = source & sMask;
  std::cout << "   s&m: ";
  print_uint64_t(temp2);

  REQUIRE((source & sMask) == ((ext >> 2) & sMask ));
}

TEST_CASE( "test if we can correctly identify the potentially homozygous stretches",
           "[test_potentially_homozygous]"){
   /*make sure all the kmers are in the graph

      This De Bruijn Graph is shaped like this
      A\                   /C
        \__1_|__2__|__3__/
        /                \
      B/                  \D

      A
      B
      1
      2
      3
      C
      D
  */
  DBG G = DBG(5);
  DBnode * pNode;
  std::string str;

  str = "CGGGG"; // canon CCCCG
  pNode = G.access_node(str, 1);
  str = "TGGGG"; // canon CCCCA
  pNode = G.access_node(str, 1);
  str =  "GGGGA"; // canon GGGGA
  pNode = G.access_node(str, 1);
  str =   "GGGAA"; // GGGAA
  pNode = G.access_node(str, 1);
  str =   "GGAAA"; // GGAAA
  pNode = G.access_node(str, 1);
  str =   "GAAAA"; // GAAAA
  pNode = G.access_node(str, 1);
  str =   "GAAAT"; // ATTTC
  pNode = G.access_node(str, 1);
  //G.print_graph();
  G.mark_branching();

  std::vector<uint64_t> vec;
  uint32_t match;
  str = "CGGGG"; //canonical CCCCG
  pNode = G.access_node(str, 0);
  REQUIRE(pNode->is_flag_on(0) == 0);
  REQUIRE(pNode->is_flag_on(1) == 0);
  vec = get_extensions(str, 5);
  match = 4;
  std::cout << "Looking at CGGGG\n";
  for (uint32_t i = 0; i <8; i++){
    if (i==match){
      std::cout << i << " " << uint64_to_kmer(vec[i], 5) << " match\n";
      REQUIRE(G.access_node(vec[i], 0) != nullptr);
    }
    else{
      std::cout << i << " " << uint64_to_kmer(vec[i], 5) << " unmatch\n";
      REQUIRE(G.access_node(vec[i], 0) == nullptr);
    }
  }

  str = "TGGGG"; // canon CCCCA
  pNode = G.access_node(str, 0);
  REQUIRE(pNode->is_flag_on(0) == 0);
  REQUIRE(pNode->is_flag_on(1) == 0);
  vec = get_extensions(str, 5);
  match = 4;
  std::cout << "Looking at TGGGG\n";
  for (uint32_t i = 0; i <8; i++){
    if (i==match){
      std::cout << i << " " << uint64_to_kmer(vec[i], 5) << " match\n";
      REQUIRE(G.access_node(vec[i], 0) != nullptr);
    }
    else{
      std::cout << i << " " << uint64_to_kmer(vec[i], 5) << " unmatch\n";
      REQUIRE(G.access_node(vec[i], 0) == nullptr);
    }
  }

  str =  "GGGGA"; // canon GGGGA - should branch 5' not 3'
  pNode = G.access_node(str, 0);
  REQUIRE(pNode->is_flag_on(0) == 1);
  REQUIRE(pNode->is_flag_on(1) == 0);
  vec = get_extensions(str, 5);
  for (uint32_t i = 0; i <8; i++){
    switch(i){
      case(1):
        REQUIRE(G.access_node(vec[i], 0) != nullptr);
        break;
      case(3):
        REQUIRE(G.access_node(vec[i], 0) != nullptr);
        break;
      case(4):
        REQUIRE(G.access_node(vec[i], 0) != nullptr);
        break;
      default:
      REQUIRE(G.access_node(vec[i], 0) == nullptr);
    }
  }

  str =   "GGGAA"; // GGGAA no branch
  pNode = G.access_node(str, 0);
  REQUIRE(pNode->is_flag_on(0) == 0);
  REQUIRE(pNode->is_flag_on(1) == 0);
  vec = get_extensions(str, 5);
  for (uint32_t i = 0; i <8; i++){
    if ((i==2) || (i==4)){
      REQUIRE(G.access_node(vec[i], 0) != nullptr);
    }
    else{
      REQUIRE(G.access_node(vec[i], 0) == nullptr);
    }
  }

  str =   "GGAAA"; // GGAAA should be 3' branching
  pNode = G.access_node(str, 0);
  REQUIRE(pNode->is_flag_on(0) == 0);
  REQUIRE(pNode->is_flag_on(1) == 1);
  vec = get_extensions(str, 5);
  std::cout << "Looking at GGAAA\n";
  for (uint32_t i = 0; i <8; i++){
    if ((i==2) || (i==4) || (i==7)){
      std::cout << i << " " << uint64_to_kmer(vec[i], 5) << " match\n";
      REQUIRE(G.access_node(vec[i], 0) != nullptr);
    }
    else{
      std::cout << i << " " << uint64_to_kmer(vec[i], 5) << " unmatch\n";
      REQUIRE(G.access_node(vec[i], 0) == nullptr);
    }
  }

  str =   "GAAAA"; // GAAAA
  pNode = G.access_node(str, 0);
  REQUIRE(pNode->is_flag_on(0) == 0);
  REQUIRE(pNode->is_flag_on(1) == 0);
  vec = get_extensions(str, 5);
  match = 2;
  for (uint32_t i = 0; i <8; i++){
    if (i==match){
      REQUIRE(G.access_node(vec[i], 0) != nullptr);
    }
    else{
      REQUIRE(G.access_node(vec[i], 0) == nullptr);
    }
  }

  str =   "GAAAT"; //canonical ATTTC
  pNode = G.access_node(str, 0);
  REQUIRE(pNode->is_flag_on(0) == 0);
  REQUIRE(pNode->is_flag_on(1) == 0);
  vec = get_extensions(str, 5);
  match = 2;
  for (uint32_t i = 0; i <8; i++){
    if (i==match){
      REQUIRE(G.access_node(vec[i], 0) != nullptr);
    }
    else{
      REQUIRE(G.access_node(vec[i], 0) == nullptr);
    }
  }

  //now mark the potentially heterozygous sections
  G.mark_non_het_for_deletion();
  str = "CGGGG"; //canonical - 3p ext to T
  pNode = G.access_node(str, 0);
  REQUIRE(pNode->is_flag_on(2) == 1);
  REQUIRE(pNode->is_flag_on(3) == 0);
  str = "TGGGG"; //canonical GTGTA
  pNode = G.access_node(str, 0);
  REQUIRE(pNode->is_flag_on(2) == 1);
  REQUIRE(pNode->is_flag_on(3) == 0);
  str =  "GGGGA"; //canonical - should branch 5' not 3'
  pNode = G.access_node(str, 0);
  REQUIRE(pNode->is_flag_on(2) == 1);
  REQUIRE(pNode->is_flag_on(3) == 1);
  str =   "GGGAA"; //canonical AAGTG - should not branch
  pNode = G.access_node(str, 0);
  REQUIRE(pNode->is_flag_on(2) == 1);
  REQUIRE(pNode->is_flag_on(3) == 1);
  str =   "GGAAA"; //canonical AAAGT - should branch 5' not 3'
  pNode = G.access_node(str, 0);
  REQUIRE(pNode->is_flag_on(2) == 1);
  REQUIRE(pNode->is_flag_on(3) == 1);
  str =   "GAAAA"; //canonical CAAAG - should not branch
  pNode = G.access_node(str, 0);
  REQUIRE(pNode->is_flag_on(2) == 1);
  REQUIRE(pNode->is_flag_on(3) == 0);
  str =   "GAAAT"; //canonical AAAAG - should not branch
  pNode = G.access_node(str, 0);
  REQUIRE(pNode->is_flag_on(2) == 1);
  REQUIRE(pNode->is_flag_on(3) == 0);

  G.delete_flagged();

  str = "CGGGG"; //canonical - 3p ext to T
  pNode = G.access_node(str, 0);
  REQUIRE(pNode->is_flag_on(2) == 1);
  REQUIRE(pNode->is_flag_on(3) == 0);
  str = "TGGGG"; //canonical GTGTA
  pNode = G.access_node(str, 0);
  REQUIRE(pNode->is_flag_on(2) == 1);
  REQUIRE(pNode->is_flag_on(3) == 0);
  str =  "GGGGA"; //canonical - should branch 5' not 3'
  pNode = G.access_node(str, 0);
  REQUIRE( pNode == nullptr);
  str =   "GGGAA"; //canonical AAGTG - should not branch
  pNode = G.access_node(str, 0);
  REQUIRE( pNode == nullptr);
  str =   "GGAAA"; //canonical AAAGT - should branch 5' not 3'
  pNode = G.access_node(str, 0);
  REQUIRE( pNode == nullptr);
  str =   "GAAAA"; //canonical CAAAG - should not branch
  pNode = G.access_node(str, 0);
  REQUIRE(pNode->is_flag_on(2) == 1);
  REQUIRE(pNode->is_flag_on(3) == 0);
  str =   "GAAAT"; //canonical AAAAG - should not branch
  pNode = G.access_node(str, 0);
  REQUIRE(pNode->is_flag_on(2) == 1);
  REQUIRE(pNode->is_flag_on(3) == 0);
}

TEST_CASE("tests if deleting kmers works without crashing everything", "[delete_many_kmers]"){
  std::set<uint64_t> visited;
  uint32_t k_size = 5;
  uint64_t canonical;
  DBnode * pNode;
  DBG G = DBG(5);
  for (uint64_t i = 0; i < 1024; i++){
    canonical = canon(i, k_size);
    if (visited.count(canonical) == 0){
      pNode = G.access_node(canonical, 1);
      REQUIRE(pNode != nullptr);
      REQUIRE(pNode->is_flag_on(3) == 0);
      pNode->bit_on(3);
      REQUIRE(pNode->is_flag_on(3) == 1);
    }
    visited.insert(canonical);
  }

  //std::cout << "Size before deleting is: " << G.size() << "\n";
  REQUIRE(G.size() == 512);
  G.delete_flagged(); //delete everything
  REQUIRE(G.size() == 0);
  //std::cout << "Size after deleting is: " << G.size() << "\n";

  visited.clear();
  for (uint64_t i = 0; i < 1024; i++){
    canonical = canon(i, k_size);
    if (visited.count(canonical) == 0){
      pNode = G.access_node(canonical, 0);
      REQUIRE(pNode == nullptr);
    }
    visited.insert(canonical);
  }
}

TEST_CASE("tests the graph printing function", "[test_graph_print]"){
  DBG G = DBG(5);
  DBnode * pNode;
  std::string str;
  str = "CGGGG"; // canon CCCCG
  pNode = G.access_node(str, 1);
  str = "TGGGG"; // canon CCCCA
  pNode = G.access_node(str, 1);
  str =  "GGGGA"; // canon GGGGA
  pNode = G.access_node(str, 1);
  str =   "GGGAA"; // GGGAA
  pNode = G.access_node(str, 1);
  str =   "GGAAA"; // GGAAA
  pNode = G.access_node(str, 1);
  str =   "GAAAA"; // GAAAA
  pNode = G.access_node(str, 1);
  str =   "GAAAT"; // ATTTC
  pNode = G.access_node(str, 1);
  std::cout << "\n\n\n";
  G.print_graph();
}

TEST_CASE("tests the HKC generating function", "[test_hkc_gen]"){
  /*Here it is important to construct a sequence that will be seeded from the middle
    GGGGGATTTTT
    */
  uint32_t ksize = 5;
  DBG G = DBG(ksize);
  DBnode * pNode;
  std::string str = "GGGGGATTTTT";
  std::string filename = "test_case_HKCs";
  for (uint32_t i = 0; i < str.size()-ksize+1; i++){
    std::string str2 = str.substr(i,ksize);
    pNode = G.access_node(str2, 1);
    REQUIRE(pNode != nullptr);
  }
  G.gen_HKCs(filename);
  std::ifstream file(filename.c_str());
  std::string line;
  while(std::getline(file, line)){
    REQUIRE(line == "AAAAATCCCCC");
    std::cout << "line: " << line << std::endl;
  }
}

TEST_CASE("tests the HKC gen number 2", "[test_hkc_gen2]"){
  /*Here it is important to construct a sequence that will be seeded from the middle
    GGGGGATTTTT
    */
  std::cout << "\n\n\nNow looking at HKC gen number 2.\n";
  uint32_t ksize = 5;
  DBG G = DBG(ksize);
  DBnode * pNode;
  std::string str = "GGAATGCTA";
  std::string filename = "test_case_HKCs";
  for (uint32_t i = 0; i < str.size()-ksize+1; i++){
    std::string str2 = str.substr(i,ksize);
    pNode = G.access_node(str2, 1);
    REQUIRE(pNode != nullptr);
  }
  G.gen_HKCs(filename);
  std::ifstream file(filename.c_str());
  std::string line;
  while(std::getline(file, line)){
    REQUIRE(line == "TAGCATTCC");
    std::cout << "line: " << line << std::endl;
  }
  //G.print_graph();
}

TEST_CASE("tests the HKC gen number 3", "[test_hkc_gen3]"){
  /*Here it is important to construct a sequence that will be seeded from the middle
    GGGGGATTTTT
    */
  std::cout << "\n\n\nNow looking at HKC gen number 2.\n";
  uint32_t ksize = 7;
  DBG G = DBG(ksize);
  DBnode * pNode;
  std::string str = "AAGTAAAATGTGACACACGTTGTCTA";
  std::string filename = "test_case_HKCs";
  for (uint32_t i = 0; i < str.size()-ksize+1; i++){
    std::string str2 = str.substr(i,ksize);
    pNode = G.access_node(str2, 1);
    REQUIRE(pNode != nullptr);
  }
  G.gen_HKCs(filename);
  std::ifstream file(filename.c_str());
  std::string line;
  while(std::getline(file, line)){
    REQUIRE(line == "AAGTAAAATGTGACACACGTTGTCTA");
    std::cout << "line: " << line << std::endl;
  }
  //G.print_graph();
}

TEST_CASE("tests the HKC gen number 4", "[test_hkc_gen4]"){
  /*Here it is important to construct a sequence that will be seeded from the middle
    GGGGGATTTTT
    */
  std::cout << "\n\n\nNow looking at HKC gen number 4.\n";
  uint32_t ksize = 13;
  DBG G = DBG(ksize);
  DBnode * pNode;
  std::string str = "TTCCATACGCCATTGGAGGAGTAACAAAATCGTTCAGTAACAGAGTGAAATGCGTCATAACAAGGGACCATCTGCCGTCAGCTCTACGTTCCAAGTACATGTTGAGTTACAA";
  std::string filename = "test_case_HKCs";
  for (uint32_t i = 0; i < str.size()-ksize+1; i++){
    std::string str2 = str.substr(i,ksize);
    pNode = G.access_node(str2, 1);
    REQUIRE(pNode != nullptr);
  }
  G.gen_HKCs(filename);
  std::ifstream file(filename.c_str());
  std::string line;
  while(std::getline(file, line)){
    REQUIRE(line == "TTCCATACGCCATTGGAGGAGTAACAAAATCGTTCAGTAACAGAGTGAAATGCGTCATAACAAGGGACCATCTGCCGTCAGCTCTACGTTCCAAGTACATGTTGAGTTACAA");
    std::cout << "line: " << line << std::endl;
  }
  REQUIRE(G.count_nulls() == 0);
  //G.print_graph();
}

TEST_CASE("some tests for the get_extensions function", "[get_extensions_str_test]"){
  std::string kmer = "GGACNATAC";
  uint32_t k = 9;
  REQUIRE( get_extensions(kmer, k+1).size() == 0);
  REQUIRE_THROWS_AS( get_extensions(kmer, k), std::runtime_error);
}

TEST_CASE("test for N in kmer", "[N_kmer_to_uint64]"){
  std::string kmer = "GGACNATAC";
  uint32_t k = 9;
  REQUIRE( kmer_to_uint64(kmer, k) == 0);
}

TEST_CASE("test if singleton double branches are identified",
          "[remove_double_branches]"){
  DBG G = DBG(5);
  std::string str = "GACAC";
  DBnode * pNode = G.access_node(str, 1);
  REQUIRE(pNode != nullptr);
  REQUIRE(pNode->is_flag_on(0) == 0);
  REQUIRE(pNode->is_flag_on(1) == 0);

  str = "TACAC";
  pNode = G.access_node(str, 1);
  REQUIRE(pNode != nullptr);
  REQUIRE(pNode->is_flag_on(0) == 0);
  REQUIRE(pNode->is_flag_on(1) == 0);

  str =  "ACACC";
  pNode = G.access_node(str, 1);
  REQUIRE(pNode != nullptr);
  REQUIRE(pNode->is_flag_on(0) == 0);
  REQUIRE(pNode->is_flag_on(1) == 0);

  str =   "CACCC";
  pNode = G.access_node(str, 1);
  REQUIRE(pNode != nullptr);
  REQUIRE(pNode->is_flag_on(0) == 0);
  REQUIRE(pNode->is_flag_on(1) == 0);

  str =   "CACCG";
  pNode = G.access_node(str, 1);
  REQUIRE(pNode != nullptr);
  REQUIRE(pNode->is_flag_on(0) == 0);
  REQUIRE(pNode->is_flag_on(1) == 0);

  G.mark_branching();
  G.mark_non_het_for_deletion();

  str = "GACAC";
  pNode = G.access_node(str, 0);
  REQUIRE(pNode != nullptr);
  REQUIRE(pNode->is_flag_on(0) == 0);
  REQUIRE(pNode->is_flag_on(1) == 0);

  str = "TACAC";
  pNode = G.access_node(str, 0);
  REQUIRE(pNode != nullptr);
  REQUIRE(pNode->is_flag_on(0) == 0);
  REQUIRE(pNode->is_flag_on(1) == 0);

  str =  "ACACC";
  pNode = G.access_node(str, 0);
  REQUIRE(pNode != nullptr);
  REQUIRE(pNode->is_flag_on(0) == 1);
  REQUIRE(pNode->is_flag_on(1) == 1);

  str =   "CACCC";
  pNode = G.access_node(str, 0);
  REQUIRE(pNode != nullptr);
  REQUIRE(pNode->is_flag_on(0) == 0);
  REQUIRE(pNode->is_flag_on(1) == 0);

  str =   "CACCG";
  pNode = G.access_node(str, 0);
  REQUIRE(pNode != nullptr);
  REQUIRE(pNode->is_flag_on(0) == 0);
  REQUIRE(pNode->is_flag_on(1) == 0);
}

TEST_CASE("test mark_non_het_for_deletion case 2",
          "[mnhfd_case2]"){
  uint32_t ksize = 5;
  DBG G = DBG(ksize);
  DBnode * pNode;
  std::string str = "CCCCG";
  std::string br1 = "CCCGA";
  std::string br2 = "CCCGT";
  pNode = G.access_node(str, 1);
  pNode = G.access_node(br1, 1);
  pNode = G.access_node(br2, 1);
  G.print_graph();
  G.mark_branching();
  G.mark_non_het_for_deletion();
  pNode = G.access_node(str, 0);
  REQUIRE(pNode->is_flag_on(0) == 0);
  REQUIRE(pNode->is_flag_on(1) == 1);
}
