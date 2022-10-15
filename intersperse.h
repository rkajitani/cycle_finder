#include "common.h"
#include "cycle.h"

class Intersperse:public Cycle{
 public:
  Intersperse();
  string kmer_align_file;
  string kmer_id_file;
  string kmer_for_path_find;
  bool intersperse_exist;
  int length_min_thr;

  void print_intersperse_usage(void);
  int option_intersperse_parse(int argc, char **argv);
  void intersperse_exe(void);
  void extract_intersperse_kmer(void);
  int path_find_parent(void);
  int map_read_intersperse(void);
};
