#include "common.h"
//#include "blastn2fa.h"
//#include "intersperse.h"

#ifndef CYCLE_H
#define CYCLE_H

class Cycle{
 public:
  Cycle();
  void cycle_exe(void);
  int option_cycle_parse(int argc, char**argv);
  void print_cycle_usage(void);
  string fa,fq,o,fq_part,fa_filter,kmer_fasta_id; //file
  int max_l; //max_length
  int max_n; //max_nodes
  int max_d; //max_depth
  int d; //data type
  int p; //peak
  int p1; //peak1
  int p2; //peak2
  int t; //num_threads
  int method;
  int repeat_type;
  double read_coverage; //read coverage for mapping(downsampling)
  string trf_out,trf_out_tandem,blast_read_out,blast_out,repeat_num_out,repeat_num_out_min;
  bool tandem_exist;
  bool bool_single, bool_compare;

  int cycle_find_parent(void);
  void make_output_file(void);
  void trf_filter(void);
  int map_read(void);
  void repeat_num(void);
  void repeat_num_exe(void);
  void remove_cycle(void);
  int blastn2fa(void);
};

#endif
