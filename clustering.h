#include "common.h"

void tandem(string repeat_num_file, string file_out);
void tandem_undo(string repeat_num_file, string file_out);
void cdhit_integrage(string repeat_num_file, string cdhit1_file, string cdhit2_file, string file_out);
void max_copy_fasta(string repeat_num_file, string cdhit_file, string file_out);
void cluster_by_blast(string repeat_num_file, string family_file, string blast_self_file, string max_copy_fasta, string file_out, int repeat_type);
