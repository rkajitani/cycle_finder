#include "cycle.h"
#include "common.h"

//constructor
Cycle::Cycle(){
  this->max_l = 1000;
  this->max_n = 10000;
  this->max_d = 10;
  this->d = 0;
  this->t = 1;
  this->o = "out_T";
  this->read_coverage = 0.5;
  this->trf_out = "";
  this->trf_out_tandem = "";
  this->blast_read_out = "";
  this->blast_out = "";
  this->repeat_num_out = "";
  this->repeat_num_out_min = "";
  this->repeat_type = 1;
  this->tandem_exist = false;

  this->bool_single = false;
  this->bool_compare = false;
}

void Cycle::print_cycle_usage(){
  cerr << "cycle usage: cycle_finder cycle [option] " << endl;
  cerr << "------------------------------------------------" << endl;
  cerr << "-f : k-mer file [string] " << endl;
  cerr << "-r : fastq file [string]" << endl;
  cerr << "-c : single->1, compare->0 [bool] (defualt:" << this->d << ")" << endl;
  cerr << "-l : length threshold [int] (defualt:" << this->max_l << ")" << endl;
  cerr << "-n : number of node threshold [int] (defualt:" << this->max_n << ")" << endl;
  cerr << "-d : depth threshold [int] (defualt:" << this->max_d << ")" << endl;
  cerr << "-p : peak  [int] " << endl;
  cerr << "-p1: peak1 [int] " << endl;
  cerr << "-p2: peak2 [int] " << endl;
  cerr << "-o : filename of output [string] (defualt:" << this->o << ")" << endl;
  cerr << "-t : num_threads [int] (default:" << this->t << ")" << endl;
  cerr << "-rc: read coverage [double] (default:" << this->read_coverage << ")" << endl;
  cerr << "-------------------------------------------------" << endl;
}

int Cycle::option_cycle_parse(int argc, char**argv){
  cout << PROGRAM_NAME;
  for(int i = 2; i < argc;){
    if(!strcmp(argv[i],"-c")){
      this->d = stoi(argv[i+1]);
      cout << " -c " << this->d;
      i = i + 2;
    } else if(!strcmp(argv[i],"-l")){
      this->max_l = stoi(argv[i+1]);
      cout << " -l " << this->max_l;
      i = i + 2;
    } else if(!strcmp(argv[i],"-n")){
      this->max_n = stoi(argv[i+1]);
      cout << " -n " << this->max_n;
      i = i + 2;
    } else if(!strcmp(argv[i],"-d")){
      this->max_d = stoi(argv[i+1]);
      cout << " -d " << this->max_d;
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
      this->fa = argv[i+1];
      cout << " -f " << this->fa;
      i = i + 2;
    } else if(!strcmp(argv[i], "-r")){
      this->fq = argv[i+1];
      cout << " -r " << this->fq;
      i = i + 2;
    } else if(!strcmp(argv[i], "-t")){
      this->t = stoi(argv[i+1]);
      cout << " -t " << this->t;
      i = i + 2;
    } else if(!strcmp(argv[i], "-o")){
      this->o = strcat(argv[i+1], "_T");
      cout << " -o " << this->o;
      i = i + 4;
    } else if(!strcmp(argv[i], "-rc")){
      this->read_coverage = stod(argv[i+1]);
      cout << " -rc " << this->read_coverage;
      i = i + 2;
    } else {
      cerr << "error : wrong option : " << argv[i] << endl;
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

void Cycle::make_output_file(){
  this->trf_out = this->o + "_trf_out";
  this->trf_out_tandem = this->o + "_trf_out_tandem";
  this->blast_read_out = this->o + "_read_out"; //blast結果
  this->fq_part = this->o + "_read_part"; //readのダウンサンプリング後のfasta
  this->fa_filter = this->o + "_filter"; //readをマップしてフィルターしたあとのk-mer
  this->repeat_num_out = this->o + "_repeat_num";
  this->repeat_num_out_min = this->o + "_repeat_num_min";
}

void Cycle::remove_cycle(){
  remove((this->trf_out).c_str());
  remove((this->trf_out_tandem).c_str());
  remove((this->fq_part).c_str());
  remove((this->blast_read_out).c_str());
  remove((this->o).c_str());
  remove((this->fa_filter).c_str());
  remove((this->o + ".log").c_str());
}

void show_cycle(){
    cout << endl;
    cout << "  +----------------------------------------------+" << endl;
    cout << "  |  Detect tandem repeats from de bruijn graph  |" << endl;
    cout << "  +----------------------------------------------+" << endl;
}

void Cycle::cycle_exe(){
  show_cycle();
  make_output_file();
  try{
    if(cycle_find_parent() == -1){
      throw "no cycle";
    }
    trf_filter();
    if(map_read() == -1){
      throw "filtered by map_read()";
    }
    repeat_num();    
    tandem_exist = true;
  }catch(const char *error_message){
    cerr << error_message << endl;
  }
//  remove_cycle();
}
