#include "intersperse.h"
#include <tuple>
#include <utility>

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

Intersperse::Intersperse(){
  this->repeat_type = 0;
  this->max_l = 5000;
  this->max_n = 10000;
  this->max_d = 5;
  this->d = 0;
  this->t = 1;
  this->o = "out_I";
  this->intersperse_exist = false;
  this->read_coverage = 0.5;
  this->length_min_thr = 100;
}

void Intersperse::print_intersperse_usage(void){
  cerr << "intersperse usage: cycle_finder intersperse [option] " << endl;
  cerr << "------------------------------------------------" << endl;
  cerr << "-f : k-mer id file [string] " << endl;
  cerr << "-b : blast file that kmer to tandem repeats [string] " << endl;
  cerr << "-r : fastq file [string]" << endl;
  cerr << "-c : single->1, compare->0 [bool] (defualt:" << this->d << ")" << endl;
  cerr << "-L : length max threshold [int] (defualt:" << this->max_l << ")" << endl;
  cerr << "-N : number of node threshold [int] (defualt:" << this->max_n << ")" << endl;
  cerr << "-D : depth threshold [int] (defualt:" << this->max_d << ")" << endl;
  cerr << "-l : length min threshold [int] (defualt:" << this->length_min_thr << ")" << endl;
  cerr << "-p : peak  [int] " << endl;
  cerr << "-p1: peak1 [int] " << endl;
  cerr << "-p2: peak2 [int] " << endl;
  cerr << "-o : filename of output [string] (defualt:" << this->o << ")" << endl;
  cerr << "-t : num_threads[int](default:" << this->t << ")" << endl;
  cerr << "-rc: read coverage [double] (default:" << this->read_coverage << ")" << endl;
  cerr << "-------------------------------------------------" << endl;
}

int Intersperse::option_intersperse_parse(int argc, char**argv){
  cout << PROGRAM_NAME << " intersperse";
  for(int i = 2; i < argc;){
    if(!strcmp(argv[i],"-c")){
      this->d = stoi(argv[i+1]);
      cout << " -c " << this->d;
      i = i + 2;
    } else if(!strcmp(argv[i],"-L")){
      this->max_l = stoi(argv[i+1]);
      cout << " -L " << this->max_l;
      i = i + 2;
    } else if(!strcmp(argv[i],"-N")){
      this->max_n = stoi(argv[i+1]);
      cout << " -N " << this->max_n;
      i = i + 2;
    } else if(!strcmp(argv[i],"-D")){
      this->max_d = stoi(argv[i+1]);
      cout << " -D " << this->max_d;
      i = i + 2;
    } else if(!strcmp(argv[i],"-l")){
      this->length_min_thr = stoi(argv[i+1]);
      cout << " -l " << this->length_min_thr;
      i = i + 2;
    } else if(!strcmp(argv[i],"-p")){
      this->p = stoi(argv[i+1]);
      cout << " -p " << this->p;
      this->d = 1;
      i = i + 2;
      this->bool_single = true;
    } else if(!strcmp(argv[i],"-p1")){
      this->p1 = stoi(argv[i+1]);
      cout << " -p1 " << this->p1;
      this->d = 0;
      i = i + 2;
      this->bool_compare = true;
    } else if(!strcmp(argv[i],"-p2")){
      this->p2 = stoi(argv[i+1]);
      cout << " -p2 " << this->p2;
      this->d = 0;
      i = i + 2;
      this->bool_compare = true;
    } else if(!strcmp(argv[i], "-f")){
      this->kmer_id_file = argv[i+1];
      cout << " -f " << this->kmer_id_file;
      i = i + 2;
    } else if(!strcmp(argv[i], "-b")){
      this->kmer_align_file = argv[i+1];
      cout << " -b " << this->kmer_align_file;
      i = i + 2;
    } else if(!strcmp(argv[i], "-r")){
      this->fq = argv[i+1];
      cout << " -r " << this->fq;
      i = i + 2;
    } else if(!strcmp(argv[i],"-rc")){
      this->read_coverage = stod(argv[i+1]);
      cout << " -rc " << this->read_coverage;    
      i = i + 2;
    } else if(!strcmp(argv[i], "-t")){
      this->t = stoi(argv[i+1]);
      cout << " -t " << this->t;
      i = i + 2;
    } else if(!strcmp(argv[i],"-o")){
      this->o = strcat(argv[i+1], "_I");
      cout << " -o " << this->o;    
      i = i + 4;
    } else {
      cerr << "error : wrong option :" << argv[i] << endl;
      return -1;
    }
  }
  cout << endl;
  if((bool_single == true && bool_compare == true) || (bool_single == false && bool_compare == false)){
    cerr << "error : input file correctly" << endl;
    return -1;
  }
  return 0;
}

void Intersperse::extract_intersperse_kmer(void){
  FILE *fp1,*fp2,*fp_out;
  char str[16384];
  string open_file1 = this->kmer_align_file;
  string open_file2 = this->kmer_id_file;
  fp1 = fopen(open_file1.c_str(),"r");
  fp2 = fopen(open_file2.c_str(),"r");
  //  this->kmer_no_cycle = this->kmer_id_file + "_no_cycle";
  this->kmer_for_path_find = this->kmer_id_file + "_no_cycle";
  this->fa = this->kmer_for_path_find;
  //  string output = this->kmer_no_cycle;
  string output = this->kmer_for_path_find;
  fp_out = fopen(output.c_str(),"w");
  if(fp1 == NULL || fp2 == NULL || fp_out == NULL){
    if(fp1 == NULL){
      cerr << "there are no k-mer alignment in tandem repeat" << endl;
    }
    if(fp2 == NULL){
      cerr << "error in extract_intersperse_kmer of intersperse.cpp fp2:" << open_file2 << endl;
      return;
    }
    if(fp_out == NULL){
      cerr << "error in extract_intersperse_kmer of intersperse.cpp fp_out:" << output << endl;
      return;
    }
  }
  string line;
  vector<string> a;
  unordered_set<string> st;
  if(fp1 != NULL){
  //read kmer alignment file
  while((fgets(str,sizeof(str),fp1)) != NULL){
    str[strlen(str)-1] = '\0';
    line = str;
    a = split(line,'\t');
    st.insert(a[1]);
  }
  fclose(fp1);
  }
  //read kmer_id fasta file
  unordered_map<string,string> mp;
  string kmer_id;
  while((fgets(str,sizeof(str),fp2)) != NULL){
    str[strlen(str)-1] = '\0';
    line = str;
    if(line[0] == '>'){
      kmer_id = line.substr(1);
    } else {
      if(st.find(kmer_id) == st.end()){
	mp.insert(make_pair(kmer_id, line));
      }
    }
  }
  fclose(fp2);

  //output
  for(auto it = mp.begin(); it != mp.end(); ++it){
    fprintf(fp_out, "%s%s\n%s\n",
	    ">", it->first.c_str(), it->second.c_str());
  }
  fclose(fp_out);
  return;
}


void show_intersperse(){
  cout << endl;
  cout << "  +--------------------------------------------------+" << endl;
  cout << "  | Detect intersperse repeats from de bruijn graph  |" << endl;
  cout << "  +--------------------------------------------------+" << endl;
}

void Intersperse::intersperse_exe(){
  show_intersperse();
  make_output_file();
  extract_intersperse_kmer();
  try{
    if(path_find_parent() == -1){
      throw "No intersperse repeats found";
    } else {
      trf_filter();
      if(map_read() == -1){
	throw "all interspersed repeats were filtered";
      }
      repeat_num();
      intersperse_exist = true;
    }
  } catch(const char *error_message) {
    cerr << error_message << endl;
  }
  remove_cycle();
}
