#include<string>

typedef struct Vars {
  uint32_t k = 0;
  uint32_t delete_hairs = 0;
  uint32_t min_count = 0;
  uint32_t max_count = 0;
  std::string out_prefix;
} Vars;

Vars process_cl_args(int argc, char **argv);
