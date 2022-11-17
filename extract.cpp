#include "extract.h"
#include "common.h"

//constructor
Extract::Extract(){
  this->c1 = 100;
  this->c1_string = "100";
  this->c2 = 10;
  this->t = 1;
  this->k = K_MER;
  this->o = "out";
  this->d = 1;
  this->jf_file = "";
  this->jf_file1 = "";
  this->jf_file2 = "";
  this->dump_file = "";
  this->dump_file1 = "";
  this->dump_file2 = "";
  this->histo_file = "";
  this->histo_file1 = "";
  this->histo_file2 = "";
  this->peak = 0;
  this->peak1 = 0;
  this->peak2 = 0;
  this->jf_file_exist = false;
  this->jf_file1_exist = false;
  this->jf_file2_exist = false;
  this->dump_file_exist = false;
  this->dump_file1_exist = false;
  this->dump_file2_exist = false;
  this->histo_file_exist = false;
  this->histo_file1_exist = false;
  this->histo_file2_exist = false;

  this->bool_single = false;
  this->bool_compare = false;
}

//variable argument
inline string str_sum(){
  return string();
}
template<typename First, typename... Rest>
string str_sum(const First& first, const Rest&... rest){
  return first + " " + str_sum(rest...) + " ";
}

void Extract::print_extract_usage(){
  cerr << "extract usage: cycle_finder extract [option]" << endl;
  cerr << "------------------------------------------------" << endl;
  cerr << "-f  : input fastq file(for single mode)[string]" << endl;
  cerr << "-f1 : input main fastq file1(for double mode)[string]" << endl;
  cerr << "-f2 : input comparing fastq file2(for double mode)[string]" << endl;
  cerr << "-c1 : copy difference(default:" << this->c1 << ")[double]" << endl;
  cerr << "-c2 : copy difference(default:" << this->c2 << ")[int] (it must be less than c1 or c1)" << endl;
  cerr << "-t  : num_threads(default:" << this->t << ")[int]                  " << endl;
  cerr << "-o  : filename of output(defualt:" << this->o << ")[string]" << endl;
  cerr << "-jf : jf file for single mode" << endl;
  cerr << "-jf1: jf file correspond to file1 for compare mode" << endl;
  cerr << "-jf2: jf file correspond to file2 for compare mode" << endl;
  cerr << "-dm : dump file for single mode" << endl;
  cerr << "-dm1: dump file correspond to file1 for compare mode" << endl;
  cerr << "-dm2: dump file correspond to file2 for compare mode" << endl;
  cerr << "-hs : histo file for single mode" << endl;
  cerr << "-hs1: histo file correspond to file1 for compare mode" << endl;
  cerr << "-hs2: histo file correspond to file2 for compare mode" << endl;
  cerr << "-------------------------------------------------" << endl;
}

int Extract::option_extract_parse(int argc, char**argv){
  cout << PROGRAM_NAME;
  for(int i = 2; i < argc;){
    if(!strcmp(argv[i],"-c1")){
      c1 = stod(argv[i+1]);
      c1_string = argv[i+1];
      cout << " -c1 " << c1;
      i = i + 2;
    } else if(!strcmp(argv[i],"-c2")){
      c2 = stoi(argv[i+1]);
      cout << " -c2 " << c2;
      i = i + 2;
    } else if(!strcmp(argv[i],"-t")){
      this->t = stoi(argv[i+1]);
      cout << " -t " << this->t;
      i = i + 2;
    } else if(!strcmp(argv[i], "-f")){
      this->files = option_multi_file(argc,argv,&i); //複数ファイルをvectorに格納
      cout << " -f";
      for(uint64_t j = 0; j < this->files.size(); ++j){
         cout << " " << this->files[j];
      }
      this->bool_single = true;
      this->d = 1; //input file;1
    } else if(!strcmp(argv[i],"-f1")){
      this->files1 = option_multi_file(argc,argv,&i); //複数ファイルをvectorに格納
      cout << " -f1";
      for(uint64_t j = 0; j < this->files1.size(); ++j){
        cout << " " << this->files1[j];
      }
      this->bool_compare = true;
      this->d = 0; //compare mode
    } else if(!strcmp(argv[i], "-f2")){
      this->files2 = option_multi_file(argc,argv,&i); //複数ファイルをvectorに格納
      cout << " -f2";
      for(uint64_t j = 0; j < this->files2.size(); ++j){
        cout << " " << this->files2[j];
      }
      this->bool_compare = true;
      this->d = 0; //compare mode
    } else if(!strcmp(argv[i],"-o")){
      this->o = argv[i+1];
      cout << " -o " << this->o;
      i = i + 2;
    } else if(!strcmp(argv[i],"-jf")){
      this->jf_file = argv[i+1];
      cout << " -jf " << this->jf_file;
      jf_file_exist = true;
      i = i + 2;
      this->bool_single = true;
    } else if(!strcmp(argv[i],"-jf1")){
      this->jf_file1 = argv[i+1];
      jf_file1_exist = true;
      cout << " -jf1 " << this->jf_file1;
      i = i + 2;
      this->bool_compare = true;
    } else if(!strcmp(argv[i],"-jf2")){
      this->jf_file2 = argv[i+1];
      jf_file2_exist = true;
      cout << " -jf2 " << this->jf_file2;
      i = i + 2;
      this->bool_compare = true;
    } else if(!strcmp(argv[i],"-dm")){
      this->dump_file = argv[i+1];
      dump_file_exist = true;
      cout << " -dm " << this->dump_file;
      i = i + 2;
      this->bool_single = true;
    } else if(!strcmp(argv[i],"-dm1")){
      this->dump_file1 = argv[i+1];
      dump_file1_exist = true;
      cout << " -dm1 " << this->dump_file1;
      i = i + 2;
      this->bool_compare = true;
    } else if(!strcmp(argv[i],"-dm2")){
      this->dump_file2 = argv[i+1];
      dump_file2_exist = true;
      cout << " -dm2 " << this->dump_file2;
      i = i + 2;
      this->bool_compare = true;
    } else if(!strcmp(argv[i],"-hs")){
      this->histo_file = argv[i+1];
      histo_file_exist = true;
      cout << " -hs " << this->histo_file;
      i = i + 2;
      this->bool_single = true;
    } else if(!strcmp(argv[i],"-hs1")){
      this->histo_file1 = argv[i+1];
      histo_file1_exist = true;
      cout << " -hs1 " << this->histo_file1;
      i = i + 2;
      this->bool_compare = true;
    } else if(!strcmp(argv[i],"-hs2")){
      this->histo_file2 = argv[i+1];
      histo_file2_exist = true;
      cout << " -hs2 " << this->histo_file2;
      i = i + 2;
      this->bool_compare = true;
    } else {
      cerr << "error : wrong option" << endl;
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

inline void skip_message(string command, string *file){
  cerr << "skip jellyfish " << command << " and use " << *file << endl;
}

inline void either_command(string command, bool file_exist, string file, string cmd){
  if(file_exist){
    skip_message(command, &file);
  } else {
    cerr << "CMD: " <<  cmd << endl;
    system(cmd.c_str());
  }
}

void Extract::jellyfish_count(){
  if(this->d == 1){ //single mode
    either_command("count", this->jf_file_exist, this->jf_file, this->cmd_count);
  } else if(this->d == 0){ //compare mode
    either_command("count", this->jf_file1_exist, this->jf_file1, this->cmd_count1);    
    either_command("count", this->jf_file2_exist, this->jf_file2, this->cmd_count2);    
  } else {
    cout << "error in jellyfish_count()" << endl;
  }
}

void Extract::jellyfish_dump(){
  if(this->d == 1){ //single mode
    either_command("dump", this->dump_file_exist, this->dump_file, this->cmd_dump);
  } else if(this->d == 0){ //compare mode
    either_command("dump", this->dump_file1_exist, this->dump_file1, this->cmd_dump1);    
    either_command("dump", this->dump_file2_exist, this->dump_file2, this->cmd_dump2);    
  } else {
    cout << "error in jellyfish_dump()" << endl;
  }
}

void Extract::jellyfish_histo(){ //jellyfish histoの実行
  if(this->d == 1){ //single mode
    either_command("histo", this->histo_file_exist, this->histo_file, this->cmd_histo);
  } else if(this->d == 0){ //compare mode
    either_command("histo", this->histo_file1_exist, this->histo_file1, this->cmd_histo1);    
    either_command("histo", this->histo_file2_exist, this->histo_file2, this->cmd_histo2);    
  } else {
    cout << "error in jellyfish_histo()" << endl;
  }
}

void Extract::cat_multi_file(){
  if(this->d == 1){
    for(uint64_t i = 0; i < this->files.size(); ++i){
      this->read_fastq += this->files[i] + " ";
    }
  } else if(this->d == 0){
    for(uint64_t i = 0; i < this->files1.size(); ++i){
      this->read_fastq1 += this->files1[i] + " ";
    }
    for(uint64_t i = 0; i < this->files2.size(); ++i){
      this->read_fastq2 += this->files2[i] + " ";
    }
  }
}

void Extract::jellyfish(){
  jellyfish_count();
  jellyfish_dump();
  jellyfish_histo();
}

void Extract::make_output_file(){ //name output file
  //jf file
  if(!jf_file_exist) this->jf_file = this->o + ".jf";
  if(!jf_file1_exist) this->jf_file1 = this->o + ".jf1";
  if(!jf_file2_exist) this->jf_file2 = this->o + ".jf2";

  //dump file
  //  if(!dump_file_exist) this->dump_file = this->o + "_" + to_string(static_cast<long long>(k)) + ".fa";
  //  if(!dump_file1_exist) this->dump_file1 = this->o + "_" + to_string(static_cast<long long>(k)) + ".fa1";
  //  if(!dump_file2_exist) this->dump_file2 = this->o + "_" + to_string(static_cast<long long>(k)) + ".fa2";
  if(!dump_file_exist) this->dump_file = this->o + ".dm";
  if(!dump_file1_exist) this->dump_file1 = this->o + ".dm1";
  if(!dump_file2_exist) this->dump_file2 = this->o + ".dm2";

  //histo file
  //  if(!histo_file_exist) this->histo_file = this->o + "_histo";
  //  if(!histo_file1_exist) this->histo_file1 = this->o + "_histo1";
  //  if(!histo_file2_exist) this->histo_file2 = this->o + "_histo2";
  if(!histo_file_exist) this->histo_file = this->o + ".hs";
  if(!histo_file1_exist) this->histo_file1 = this->o + ".hs1";
  if(!histo_file2_exist) this->histo_file2 = this->o + ".hs2";
  this->kmer_peak_file = this->o + ".kmer_peak";

  //kmer file
  //    this->kmer_for_cycle_find = this->o + "_" + to_string((long long)k) + "mer_" + to_string((long long)this->c2)+"_for_detect";
  //    this->kmer_for_copy_estimate = this->o + "_" + to_string((long long)k) + "mer_" + this->c1_string + "_for_estimate";
  this->kmer_for_cycle_find = this->o + "_for_detect" + to_string((long long)this->c2) + ".fa";
  this->kmer_for_copy_estimate = this->o + "_for_estimate" + this->c1_string +".fa";

  //cat fastq file
  cat_multi_file();

  //jellyfish count command
  this->cmd_count = "jellyfish count -C -m " + to_string(static_cast<long long>(this->k)) + " -o " + this->jf_file + " -t " + to_string(static_cast<long long>(this->t)) + " -s 100M " + this->read_fastq;
  this->cmd_count1 = "jellyfish count -C -m " + to_string(static_cast<long long>(this->k)) + " -o " + this->jf_file1 + " -t " + to_string(static_cast<long long>(this->t)) + " -s 100M " + this->read_fastq1;
  this->cmd_count2 = "jellyfish count -C -m " + to_string(static_cast<long long>(this->k)) + " -o " + this->jf_file2 + " -t " + to_string(static_cast<long long>(this->t)) + " -s 100M " + this->read_fastq2;

  //jellyfish dump command
  this->cmd_dump = "jellyfish dump " + this->jf_file + " > " +this->dump_file;
  this->cmd_dump1 = "jellyfish dump " + this->jf_file1 + " > " +this->dump_file1;
  this->cmd_dump2 = "jellyfish dump " + this->jf_file2 + " > " +this->dump_file2;

  //jellyfish histo command
  this->cmd_histo = "jellyfish histo " + this->jf_file + " -l 1 -o " + histo_file;
  this->cmd_histo1 = "jellyfish histo " + this->jf_file1 + " -l 1 -o " + histo_file1;
  this->cmd_histo2 = "jellyfish histo " + this->jf_file2 + " -l 1 -o " + histo_file2;
}

void show_extract(){
  cout << endl;
  cout << "  +--------------------------------+" << endl;
  cout << "  |  Extract high frequency k-mer  |" << endl;
  cout << "  +--------------------------------+" << endl;
}

void Extract::extract_exe(){
  show_extract();
  make_output_file();
  jellyfish();
  peak_detect();
  kmer_compare();
}
