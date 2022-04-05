#include "common.h"
#ifndef EXTRACT_H
#define EXTRACT_H


class Extract{
 public:
  Extract();
  void extract_exe();
  int option_extract_parse(int argc, char**argv);
  void print_extract_usage();
  string o,o1, kmer_for_cycle_find, kmer_for_copy_estimate,kmer_peak_file; //output
  vector<std::string> files,files1,files2; //files
  double c1; //copy difference
  string c1_string;
  uint64_t c2; //copy difference
  uint64_t t; //num_threads
  uint64_t k; //k-mer
  bool d; //data type
  string dump, dump1, dump2; //dump file(fasta)
  string read_fastq;
  string read_fastq1;
  string read_fastq2;

  string jf_file, jf_file1, jf_file2;
  string dump_file, dump_file1, dump_file2;
  string histo_file, histo_file1, histo_file2;

  string cmd_count, cmd_count1, cmd_count2, cmd_dump, cmd_dump1, cmd_dump2, cmd_histo, cmd_histo1, cmd_histo2;

  bool jf_file_exist, jf_file1_exist, jf_file2_exist;
  bool dump_file_exist, dump_file1_exist, dump_file2_exist;
  bool histo_file_exist, histo_file1_exist, histo_file2_exist;

  bool bool_single, bool_compare;

  uint64_t peak, peak1, peak2;
  void cat_multi_file();
  void make_output_file();
  void jellyfish_count();
  void jellyfish_dump();
  void jellyfish_histo();
  void jellyfish();
  void kmer_compare();
  void peak_detect();
  void write_kmer_fasta(unordered_map<bitset<K_MER*2>, uint64_t> mp1);
  void write_kmer_fasta(unordered_map<bitset<K_MER*2>, uint64_t> mp1, unordered_map<bitset<K_MER*2>, uint64_t> mp2);
  void write_peak_to_file();
};

#endif
