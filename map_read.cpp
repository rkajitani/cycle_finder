#include "cycle.h"
#include "intersperse.h"
//#define COVERAGE 0.5
#define READ_LEN 100
#define EXTRACT_MIN_LEN 20

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

int Cycle::map_read(){
  string trf_out_tandem = this->trf_out_tandem;
  string blast_read_out = this->blast_read_out;
  string fastq = this->fq;
  int peak;
  if(this->d == 1){
    peak = this->p;
  } else if(this->d == 0){
    peak = this->p1;
  } else {
    cerr <<  "error : data type is incorrect" << endl;
    return -1;
  }
  int num_threads = this->t;
  //  string fq_part = this->o + "_read_part";
  string fq_part = this->fq_part;
  FILE *fp1,*out_file;
  char str[16384];
  string line;
  fp1 = fopen(fastq.c_str(),"r");
  if(fp1 == NULL){
    cerr << "ERROR : COULDN'T OPEN FASTQ FILE FOR MAPPING READS TO REPEATS" << endl;
    return -1;
  }

  //fastqの行数を得る  
  uint64_t line_num = 0, extract_num = 0;
  while((fgets(str,sizeof(str),fp1)) != NULL){
    ++line_num;
  }
  fseek(fp1,0,SEEK_SET);

  int depth = (peak * READ_LEN)/(READ_LEN - K_MER + 1);
  //  extract_num = (line_num/4)*COVERAGE/depth;
  cout << "read_coverage:" << this->read_coverage << endl;
  extract_num = (line_num/4)*(this->read_coverage)/depth;
  cout << "\t" << "estimated depth:" << depth << "\t" << "using read number:" << extract_num << "\t" << "k-mer peak:" << peak << endl;
  uint64_t extract_per = line_num/(extract_num*4); //extract_per本に1本配列を抽出
  uint64_t counter = 0, counter_extract = 0; //counter:読み込んだ行数、　counter_extract:抽出した本数
  unordered_map<string,string> mp_f,mp_r;
  vector<string> a,Seq;
  out_file = fopen(fq_part.c_str(),"w");
  while((fgets(str,sizeof(str),fp1)) != NULL){
    if(counter >= 4 * extract_per * counter_extract){
      if(counter % 4 == 1 && strlen(str) > EXTRACT_MIN_LEN){
	++counter_extract;
	fprintf(out_file,">%ld\n%s", counter_extract ,str);
	if(counter_extract > extract_num){
	  break;
	}
      }
    }
    ++counter;
  }
  cout << "\t" << "extracted read number:" << counter_extract << endl;
  fclose(fp1);
  fclose(out_file);

  //makeblastdb
  makeblastdb(fq_part,this->o+".log");

  //blast
  string cmd_blast = "blastn -db " + fq_part + " -query " + trf_out_tandem 
    + " -outfmt 6 -task megablast -out " + blast_read_out + " -dust no -num_threads " + to_string(static_cast<long long>(num_threads));
  cerr << "CMD: " << cmd_blast << endl;
  system(cmd_blast.c_str());

  //rm fq_part.*
  string cmd_rm_index = "rm " + fq_part + ".*";
  cerr << "CMD: " << cmd_rm_index << endl;
  system(cmd_rm_index.c_str());

  if(file_exist(blast_read_out) == 0){ //readにblastした結果のファイルがあるか
    return -1;
  } else {
    //blast結果->k-mer fasta
    if(blastn2fa() == -1){
      return -1;
    }
  }
  return 0;
}

int Intersperse::map_read_intersperse(){
  string trf_out = this->trf_out;
  string blast_read_out = this->blast_read_out;
  string fastq = this->fq;
  int peak;
  if(this->d == 1){
    peak = this->p;
  } else if(this->d == 0){
    peak = this->p1;
  } else {
    cerr <<  "error : data type is incorrect" << endl;
    return -1;
  }
  int num_threads = this->t;
  //  string fq_part = this->o + "_read_part";
  string fq_part = this->fq_part;
  FILE *fp1,*out_file;
  char str[16384];
  string line;
  fp1 = fopen(fastq.c_str(),"r");
  if(fp1 == NULL){
    cerr << "ERROR : COULDN'T OPEN FASTQ FILE FOR MAPPING READS TO REPEATS" << endl;
    return -1;
  }

  //fastqの行数を得る  
  uint64_t line_num = 0, extract_num = 0;
  while((fgets(str,sizeof(str),fp1)) != NULL){
    ++line_num;
  }
  fseek(fp1,0,SEEK_SET);

  int depth = (peak * READ_LEN)/(READ_LEN - K_MER + 1);
  //  extract_num = (line_num/4)*COVERAGE/depth;
  cout << "read_coverage:" << this->read_coverage << endl;
  extract_num = (line_num/4)*(this->read_coverage)/depth;
  cout << "\t" << "estimated depth:" << depth << "\t" << "using read number:" << extract_num << "\t" << "k-mer peak:" << peak << endl;
  uint64_t extract_per = line_num/(extract_num*4); //extract_per本に1本配列を抽出
  uint64_t counter = 0, counter_extract = 0; //counter:読み込んだ行数、　counter_extract:抽出した本数
  unordered_map<string,string> mp_f,mp_r;
  vector<string> a,Seq;
  out_file = fopen(fq_part.c_str(),"w");
  while((fgets(str,sizeof(str),fp1)) != NULL){
    if(counter >= 4 * extract_per * counter_extract){
      if(counter % 4 == 1 && strlen(str) > EXTRACT_MIN_LEN){
	++counter_extract;
	fprintf(out_file,">%ld\n%s", counter_extract ,str);
	if(counter_extract > extract_num){
	  break;
	}
      }
    }
    ++counter;
  }
  cout << "\t" << "extracted read number:" << counter_extract << endl;
  fclose(fp1);
  fclose(out_file);

  //makeblastdb
  makeblastdb(fq_part,this->o+".log");

  ifstream query_ifs(trf_out);
  unsigned word_size = 28;
  while (getline(query_ifs, line)) {
    if (line.size() > 0 && line[0] != '>' && line.size() < word_size) {
      word_size = line.size();
    }
  }
  query_ifs.close();

  //blast
  string cmd_blast = "blastn -db " + fq_part + " -query " + trf_out
    + " -outfmt 6 -task megablast -word_size " + to_string(word_size) + " -out " + blast_read_out + " -dust no -num_threads " + to_string(static_cast<long long>(num_threads));
  cerr << "CMD: " << cmd_blast << endl;
  system(cmd_blast.c_str());

  //rm fq_part.*
  string cmd_rm_index = "rm " + fq_part + ".*";
  cerr << "CMD: " << cmd_rm_index << endl;
  system(cmd_rm_index.c_str());

  if(file_exist(blast_read_out) == 0){ //readにblastした結果のファイルがあるか
    return -1;
  } else {
    //blast結果->k-mer fasta
    if(blastn2fa() == -1){
      return -1;
    }
  }
  return 0;
}
