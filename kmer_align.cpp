#include "cluster.h"

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

struct OUTPUT_ELEMENT{
  string query_name,db_name;
  int snp_num;
};

bool asc( const OUTPUT_ELEMENT& left, const OUTPUT_ELEMENT& right ) { //構造体OUTPUT_ELEMENTの比較関数
  return left.query_name < right.query_name;
}


  //排他的論理和をとったもので何個の変異があるかを出力
inline int kmer_dif(bitset<K_MER*2> &bs1, bitset<K_MER*2> &bs2, int &mismatch){
  bitset<K_MER*2> bs4,bs5,bs_3,bs_0;
  int number = 0;
  bs_3 = 3; //00~11
  bs_0 = 0; //00~00
  for(int i = 0; i < K_MER; ++i){
    bs4 = (bs1 ^ bs2) >> (K_MER - i)*2;
    bs5 = bs4 & bs_3;
    if(bs5 != bs_0){
      ++number;
      if(number > mismatch){
	return (mismatch+1);
      }
    }
  }
  return number;
}

void Cluster::kmer_align(){
  FILE *fp1,*fp2;
  char str[16384]; //1行16384文字までしか読み込めない
  string line,open_file1,open_file2,output;
  open_file1 = this->f1;
  open_file2 = this->blast_target;
  output = this->kmer_align_file;
  int mismatch_thr = this->mismatch;
  cout << "allowance of mismatch is " << mismatch_thr << " or less than " << mismatch_thr << endl;
  int num_threads = this->t;
  fp1 = fopen(open_file1.c_str(),"r");
  fp2 = fopen(open_file2.c_str(),"r");

  if(fp1 == NULL){
    cout << "cant't open " << open_file1 << " in kmer_align.cpp" << endl;
    return;
  } else if(fp2 == NULL){
    cout << "cant't open " << open_file2 << " in kmer_align.cpp" << endl;
    return;
  }

  string seq_name;
  unordered_map<bitset<K_MER*2>,string> mp_kmer;
  cout << "reading " << open_file1 << "..." << endl;
  while((fgets(str,sizeof(str),fp1)) != NULL){
    str[strlen(str)-1] = '\0';
    line = str;
    if(line[0] == '>'){
      seq_name = line.substr(1);
    } else {
      mp_kmer.insert(make_pair(str_to_bit(line),seq_name));
    }
  }
  fclose(fp1);

  unordered_map<string,string> mp_repeat; //配列名、配列
  vector<string> Repeat; //repeatを順番通りに格納
  while((fgets(str,sizeof(str),fp2)) != NULL){
    str[strlen(str)-1] = '\0';
    line = str;
    if(line[0] == '>'){
      seq_name = line.substr(1);
    } else {
      mp_repeat.insert(make_pair(seq_name,line));
      Repeat.push_back(seq_name);
    }
  }
  fclose(fp2);

  bitset<K_MER*2> bs1,bs2,bs3,bs4;
  int dif_num = 0;
  string repeat_str;
  bitset<K_MER*2> repeat_str_kmer,repeat_str_kmer_com;
  OUTPUT_ELEMENT output_element;
  vector<OUTPUT_ELEMENT> v_out,v_out_tmp;
  v_out.reserve(mp_kmer.size());
  v_out_tmp.reserve(mp_kmer.size());
  omp_set_num_threads(num_threads);
#pragma omp parallel
{
#pragma omp for firstprivate(repeat_str_kmer,repeat_str_kmer_com,repeat_str,dif_num,bs1,bs2,bs3,bs4,output_element,mp_kmer,mp_repeat,v_out_tmp)
  for(unsigned i_repeat = 0; i_repeat < Repeat.size(); ++i_repeat){
      repeat_str = mp_repeat[Repeat[i_repeat]]; //タンデムのリピート配列
      repeat_str_kmer = str_to_bit(repeat_str.substr(0,K_MER)); //一番初めのkmer
      repeat_str_kmer >>= 2; //次のfor文では左に2シフトさせるので、初めのk-merも次のfor文で扱うために右に2シフトしておく
      for(unsigned i = 0; i < repeat_str.size() - K_MER + 1; ++i){ //bit演算でkmerを抽出していく
	switch(repeat_str[K_MER + i - 1]){
	case 'A':
	  repeat_str_kmer <<= 2;
	  break;
	case 'G':
	  (repeat_str_kmer <<= 2) |=1;
	  break;
	case 'C':
	  (repeat_str_kmer <<= 2) |=2;
	  break;
	case 'T':
	  (repeat_str_kmer <<= 2) |=3;
	  break;
	default:
	  cout << "including N" << endl;
	  break;
	}
	repeat_str_kmer_com = bit_to_bit_com(repeat_str_kmer);
	for(auto it = mp_kmer.begin(); it != mp_kmer.end(); ++it){
	  bs4 = it->first;
	  dif_num = kmer_dif(repeat_str_kmer, bs4, mismatch_thr);
	  if(dif_num <= mismatch_thr){
	    output_element = {Repeat[i_repeat], it->second, dif_num};
	    v_out_tmp.push_back(output_element);
	  }
	  dif_num = kmer_dif(repeat_str_kmer_com, bs4, mismatch_thr);
	  if(dif_num <= mismatch_thr){
	    output_element = {Repeat[i_repeat], it->second, dif_num};
	    v_out_tmp.push_back(output_element);
	  }
	}
#pragma omp critical (write)
	{
	  v_out.insert(v_out.end(), v_out_tmp.begin(), v_out_tmp.end());
	  v_out_tmp.clear();
	}
      }
    }
  }
 
  //sort
  cout << "\t" << "sorting the result of alignment..." << endl;
  sort(v_out.begin(), v_out.end(), asc);

  //出力
  FILE *fp_out;
  fp_out = fopen(output.c_str(),"w");
  for(unsigned i = 0; i < v_out.size(); ++i){
    output_element = v_out[i];
    dif_num = output_element.snp_num;
    //    fprintf(fp_out, "%s\t%s\t%f\t%d\t%d\t%d\t%s\t%s\t%s\t%s\t%s\t%d\n",
    //	    output_element.query_name.c_str(), output_element.db_name.c_str(), 100*(double)(K_MER-dif_num)/K_MER,K_MER-dif_num,
    //	    0, 0, "column7", "column8", "column9", "column10", "column11", K_MER-dif_num);
    fprintf(fp_out, "%s\t%s\t%d\n",
	    output_element.query_name.c_str(), output_element.db_name.c_str(), K_MER-dif_num);
  }
  fclose(fp_out);
}
