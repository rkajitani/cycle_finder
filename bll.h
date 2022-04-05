#include "common.h"
#include "extract.h"
//#include "cycle.h"
#include "cluster.h"
#include "intersperse.h"

//class All:public Extract,public Cycle{
class All{
 public:
  All();
  void all_exe(int argc, char**argv);
  string o,o1,kmer_for_cycle_find,kmer_for_copy_estimate,kmer_for_path_find; //output
  vector<string> files,files1,files2; //files
  double c1; //copy difference
  string c1_string;
  uint64_t c2; //copy difference
  uint64_t t; //num_threads
  bool d; //data type
  int repeat_type;

  //mode
  bool bool_single, bool_compare;
  
  //jellyfish_file
  bool jf_file_exist, jf_file1_exist, jf_file2_exist;
  bool dump_file_exist, dump_file1_exist, dump_file2_exist;
  bool histo_file_exist, histo_file1_exist, histo_file2_exist;
  
  //stopwatch
  bool bool_stopwatch;

  //tandem
  bool tandem_exist;
  int max_l; //max_length
  int max_n; //max_nodes
  int max_d; //max_depth

  //cluster
  int mismatch;

  //intersperse
  int max_L; //max_length
  int max_N; //max_nodes
  int max_D; //max_depth

  //peak
  uint64_t p;
  uint64_t p1;
  uint64_t p2;

  double read_coverage;
  string jf_file, jf_file1, jf_file2;
  string dump_file, dump_file1, dump_file2;
  string histo_file, histo_file1, histo_file2;

  void print_all_usage(void);
  void option_all_init(void);
  int option_all_parse(int argc, char**argv);
  void opt_all_to_extract(Extract *extract);
  void peak_to_opt_all(uint64_t *peak, uint64_t *peak1, uint64_t *peak2);
  void opt_all_to_cycle(Cycle *cycle);
  void opt_all_to_cluster(Cluster *cluster);
  void opt_all_to_intersperse(Intersperse *intersperse);
};
