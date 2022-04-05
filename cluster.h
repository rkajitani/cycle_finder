#include "common.h"

class Cluster{
 public:
  Cluster();
  string f1,f2,f3,o,kmer_align_file,blast_target; //file
  int r_t; //repeat_type(1:tandem, 0:intersperse)
  int p; //peak
  int d; //data type
  int p1; //peak1
  int p2; //peak2
  int t; //num_threads
  int mismatch; //the threshold mismatch of kmer alignment
  bool bool_single, bool_compare;

  //function
  void option_cluster_init();
  void print_cluster_usage();
  int option_cluster_parse(int argc, char**argv);
  void cluster_exe();
  void cdhit(string fasta_file, string file_out, double align_cov);
  void kmer_align();
  void blast_to_copynum();
  void tandem(string repeat_num_file, string file_out);
  void tandem_undo(string repeat_num_file, string file_out);
  void cdhit_integrage(string repeat_num_file, string cdhit1_file, string cdhit2_file, string file_out);
  void max_copy_fasta(string repeat_num_file, string cdhit_file, string file_out);
  void cluster_by_blast(string repeat_num_file, string family_file, string blast_self_file, string max_copy_fasta, string file_out,int repeat_type);
  void clst_to_family(string open_file1, string open_file2, string output);
  void cluster_tandem();
  void cluster_intersperse();
  void make_output_file();
  void cluster_cdhit();
  void copy_estimation();
  void cluster_blast();
  void remove_cluster();

  //file_name
  string repeat_num_tandem,repeat_num_tandem_clst,repeat_num_tandem_clst_max,repeat_num_not_tandem_clst,repeat_cluster,kmer_fasta,represent_fa,blast_self,cmd_blast_self,len_thr,ide_thr,filter_blast_self,inherit_tsv,inherit_fa,blst_family,blst_family_clst,blst_family_clst_fa,blst_family_clst_element,final_output_tsv,final_output_fa,final_output_element,blst_self;
};
