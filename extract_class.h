#include "common.h"
#include "peak_detect.h"
#include "kmer_compare.h"

class Extract{
  typedef struct{ //option
    //  string f,f1,f2,o; //file
    string o; //output
    vector<string> files,files1,files2; //files
    double c1; //copy difference
    uint64_t c2; //copy difference
    uint64_t t; //num_threads
    uint64_t n_f; //the number of file
    uint64_t k; //k-mer
    bool d; //data type
    string jf; //jf file
  } option_extract_t;

  void print_extract_usage(void);
  void option_extract_init(option_extract_t *opt);
  void option_extract_destroy(const option_extract_t *opt);
  uint8_t option_extract_parse(option_extract_t *opt, int argc, char**argv);
  void extract_exe(const option_extract_t *opt, uint64_t *peak, uint64_t *peak1, uint64_t *peak2);
}
