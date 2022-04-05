#include "common.h"
#include "extract.h"

#define PEAK_RATIO 0.5 //threshold of k-mer copy number for estimating copy number

namespace std{
  template <>
  class hash<std::bitset<K_MER*2> >
  {
  public:
    size_t operator()(const std::bitset<K_MER*2>& x) const
    {
      size_t value=0;
      for(int i = 0; i < K_MER*2; ++i){
        if(x[i] == 1){
          value = ((value << 1) | 1);
        } else {
          value = value << 1;
        }
      }
      return value;
    }
  };
}

inline void kmer_fasta_store_part(string file, unordered_map<bitset<K_MER*2>, uint64_t> &mp, const int peak, const double thr_num1){ //store the information of limited kmer and frequency
  FILE *fp;
  char str[1024];
  string line;
  fp = fopen(file.c_str(),"r");
  if(fp == NULL){
    cout << "can't open " << file << " in kmer_compare.cpp" << endl;
    return;
  }
  unsigned long long number = 0;
  while((fgets(str,sizeof(str),fp)) != NULL){
      str[strlen(str)-1] = '\0';
      line = str;
      if(line[0] == '>'){
	number = stoll(line.substr(1,line.size()-1));
      } else if(number >= thr_num1 * peak){ //store the k-mer only when the frequency of k-mer is more than thr_num1
	mp.insert(make_pair(str_to_bit(line),number));
      }
  }
  fclose(fp);
}

inline void kmer_fasta_store_all(string file, unordered_map<bitset<K_MER*2>, uint64_t> &mp, const int peak, const double thr_num1){ //store the information of all kmer and frequency
  FILE *fp;
  char str[1024];
  string line;
  fp = fopen(file.c_str(),"r");
  if(fp == NULL){
    cout << "can't open " << file << " in kmer_compare.cpp" << endl;
    return;
  }
  unsigned long long number = 0;
  while((fgets(str,sizeof(str),fp)) != NULL){
      str[strlen(str)-1] = '\0';
      line = str;
      if(line[0] == '>'){
	number = stoll(line.substr(1,line.size()-1));
      } else { //store all the k-mer
	mp.insert(make_pair(str_to_bit(line),number));
      }
  }
  fclose(fp);
}

void Extract::write_kmer_fasta(unordered_map<bitset<K_MER*2>,uint64_t> mp1){ //write two types of kmer for which cycle finding and copy esitimation when single mode
  string output = this->kmer_for_cycle_find;
  string output1 = this->kmer_for_copy_estimate;
  uint64_t peak = this->peak;
  double thr_num1 = this->c1;
  int thr_num2 = this->c2;
  FILE *fp_out,*fp_out1;
  fp_out = fopen(output.c_str(),"w");
  fp_out1 = fopen(output1.c_str(),"w");
  uint64_t counter = 0;
  for(auto it = mp1.begin(); it != mp1.end(); ++it){
    if(it->second >= thr_num1 * peak){ //まずpeak以上かcheck
      ++counter;
      fprintf(fp_out1,"%s%lu%s%lu\n%s\n",
	      ">", it->second, "_", counter, bit_to_str(it->first).c_str());
      if(it->second  >= thr_num2 * peak){ //閾値を超えてるかをcheck
	fprintf(fp_out,"%s%lu%s%lu\n%s\n",
		">", it->second, "_", counter, bit_to_str(it->first).c_str());
      }
    }
  }
  fclose(fp_out);
  fclose(fp_out1);
}

void Extract::write_kmer_fasta(unordered_map<bitset<K_MER*2>,uint64_t> mp1, unordered_map<bitset<K_MER*2>,uint64_t> mp2){ //write two types of kmer for which cycle finding and copy esitimation when compare mode
  string output = this->kmer_for_cycle_find;
  string output1 = this->kmer_for_copy_estimate;
  uint64_t peak1 = this->peak1;
  uint64_t peak2 = this->peak2;
  double thr_num1 = this->c1;
  int thr_num2 = this->c2;
  FILE *fp_out,*fp_out1;
  fp_out = fopen(output.c_str(),"w");
  fp_out1 = fopen(output1.c_str(),"w");
  uint64_t counter = 0;
  for(auto it = mp1.begin(); it != mp1.end(); ++it){
    if(mp2.find(it->first) != mp2.end()){ //mp2に格納されているなら
      if(it->second * peak2 >= mp2[it->first] * peak1 + thr_num1 * peak1 * peak2){ //over thr_num1
	++counter;
	fprintf(fp_out1,"%s%lu%s%lu%s%lu\n%s\n",
		">Testis", it->second, "Somatic", mp2[it->first], "_", counter, bit_to_str(it->first).c_str());
	if(it->second * peak2 >= mp2[it->first] * peak1 + thr_num2 * peak1 * peak2){ //over thr_num2
	fprintf(fp_out,"%s%lu%s%lu%s%lu\n%s\n",
		">Testis", it->second, "Somatic", mp2[it->first], "_", counter, bit_to_str(it->first).c_str());
	}
      }
    } else if(it->second >= thr_num1 * peak1){ //まずpeak以上なら
	++counter;
	fprintf(fp_out1,"%s%lu%s%lu\n%s\n",
		">Testis", it->second, "Somatic0_", counter, bit_to_str(it->first).c_str());
	if(it->second >= thr_num2 * peak1){ //mp2になく,閾値以上ならsomaでは0
	  fprintf(fp_out,"%s%lu%s%lu\n%s\n",
		  ">Testis", it->second, "Somatic0_", counter, bit_to_str(it->first).c_str());
	}
    }
  }
  fclose(fp_out);
  fclose(fp_out1);
}

void Extract::kmer_compare(){
  unordered_map<bitset<K_MER*2>, uint64_t> mp1,mp2;
  if(this->d == 1){
    cout << "\t" << "extracting high frequency k-mer..." << endl;
    kmer_fasta_store_part(this->dump_file, mp1, this->peak, this->c1);
    write_kmer_fasta(mp1);
    mp1.clear();
  } else if(this->d == 0){
    cout << "\t" << "comparing k-mer...(c1:" << this->c1 << " c2:" << this->c2 << ")" << endl;
    kmer_fasta_store_part(this->dump_file1, mp1, this->peak1, this->c1);
    kmer_fasta_store_all(this->dump_file2, mp2, this->peak2, this->c1);
    write_kmer_fasta(mp1, mp2);
    mp1.clear();
    mp2.clear();
  } else {
    cout << "error" << endl;
  }
}
