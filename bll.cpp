#include "bll.h"

//constructor
All::All(){
  this->c1 = 10;
  this->c1_string = "10";
  this->c2 = 10;
  this->t = 1;
  this->o = "out";
  this->d = 1;
  this->mismatch = 3;

  this->max_l = 1000;
  this->max_n = 10000;
  this->max_d = 10;
  this->max_L = 5000;
  this->max_N = 10000;
  this->max_D = 5;

  this->read_coverage = 0.5;
  this->repeat_type = -1;

  this->jf_file = "";
  this->jf_file1 = "";
  this->jf_file2 = "";

  this->dump_file = "";
  this->dump_file1 = "";
  this->dump_file2 = "";

  this->histo_file = "";
  this->histo_file1 = "";
  this->histo_file2 = "";

  this->jf_file_exist = false;
  this->jf_file1_exist = false;
  this->jf_file2_exist = false;

  this->dump_file_exist = false;
  this->dump_file1_exist = false;
  this->dump_file2_exist = false;

  this->histo_file_exist = false;
  this->histo_file1_exist = false;
  this->histo_file2_exist = false;

  this->tandem_exist = false;

  this->bool_stopwatch = false;

  this->bool_single = false;
  this->bool_compare = false;
}

void All::print_all_usage(){
  cerr << "all usage: cycle_finder all [option]" << endl;
  cerr << "------------------------------------------------" << endl;
  cerr << "-f   : input fastq file(for single mode)[string]" << endl;
  cerr << "-f1  : input main fastq file1(for double mode)[string]" << endl;
  cerr << "-f2  : input comparing fastq file2(for double mode)[string]" << endl;
  cerr << "-c1  : copy difference (default:" << this->c1 << ") [double]" << endl;
  cerr << "-c2  : copy difference (default:" << this->c2 << ") [int]" << endl;
  cerr << "-m   : the threshold of mismatch in kmer alignment[int] (default:" << this->mismatch << ") " << endl;
  cerr << "-t   : num_threads (default:" << this->t << ") [int]" << endl;
  cerr << "-o   : filename of output (defualt:" << this->o << ") [string]" << endl;
  cerr << "-l   : length threshold [int] (defualt:" << this->max_l << ")" << endl;
  cerr << "-n   : number of node threshold [int] (defualt:" << this->max_n << ")" << endl;
  cerr << "-d   : depth threshold [int] (defualt:" << this->max_d << ")" << endl;
  cerr << "-L   : length threshold [int] (defualt:" << this->max_L << ")" << endl;
  cerr << "-N   : number of node threshold [int] (defualt:" << this->max_N << ")" << endl;
  cerr << "-D   : depth threshold [int] (defualt:" << this->max_D << ")" << endl;
  cerr << "-rc  : read coverage(downsampling)  [double] (defualt:" << this->read_coverage << ")" << endl;
  cerr << "-jf  : jf file for single mode [string]" << endl;
  cerr << "-jf1 : jf file correspond to file1 for compare mode [string]" << endl;
  cerr << "-jf2 : jf file correspond to file2 for compare mode [string]" << endl;
  cerr << "-dm  : dump file for single mode" << endl;
  cerr << "-dm1 : dump file correspond to file1 for compare mode" << endl;
  cerr << "-dm2 : dump file correspond to file2 for compare mode" << endl;
  cerr << "-hs  : histo file for single mode" << endl;
  cerr << "-hs1 : histo file correspond to file1 for compare mode" << endl;
  cerr << "-hs2 : histo file correspond to file2 for compare mode" << endl;
  cerr << "-time: use stopwatch" << endl;
  cerr << "-------------------------------------------------" << endl;
}

int All::option_all_parse(int argc, char**argv){
  cout << PROGRAM_NAME << " all";
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
      t = stoi(argv[i+1]);
      cout << " -t " << t;
      i = i + 2;
    } else if(!strcmp(argv[i], "-f")){
      files = option_multi_file(argc,argv,&i); //複数ファイルをvectorに格納
      if(files.size() == 0) cerr << "error:input more than 1 file" << endl;
      cout << " -f";
      for(unsigned j = 0; j < files.size(); ++j){
	cout << " " << files[j];
      }
      bool_single = true;
      d = 1; //input file;1
    } else if(!strcmp(argv[i],"-f1")){
      files1 = option_multi_file(argc,argv,&i); //複数ファイルをvectorに格納
      if(files1.size() == 0) cerr << "error:input more than 1 file" << endl;
      cout << " -f1";
      for(unsigned j = 0; j < files1.size(); ++j){
	cout << " " << files1[j];
      }
      bool_compare = true;
      d = 0; //compare mode
    } else if(!strcmp(argv[i], "-f2")){
      files2 = option_multi_file(argc,argv,&i); //複数ファイルをvectorに格納
      if(files2.size() == 0) cerr << "error:input more than 1 file" << endl;
      cout << " -f2";
      for(unsigned j = 0; j < files2.size(); ++j){
	cout << " " << files2[j];
      }
      bool_compare = true;
      d = 0; //compare mode
    } else if(!strcmp(argv[i],"-o")){
      o = argv[i+1];
      cout << " -o " << o;
      i = i + 2;
    } else if(!strcmp(argv[i],"-l")){
      max_l = stoi(argv[i+1]);
      cout << " -l " << max_l;
      i = i + 2;
    } else if(!strcmp(argv[i],"-n")){
      max_n = stoi(argv[i+1]);
      cout << " -n " << max_n;
      i = i + 2;
    } else if(!strcmp(argv[i],"-d")){
      max_d = stoi(argv[i+1]);
      cout << " -d " << max_d;
      i = i + 2;
    } else if(!strcmp(argv[i],"-L")){
      max_L = stoi(argv[i+1]);
      cout << " -L " << max_L;
      i = i + 2;
    } else if(!strcmp(argv[i],"-N")){
      max_N = stoi(argv[i+1]);
      cout << " -N " << max_N;
      i = i + 2;
    } else if(!strcmp(argv[i],"-D")){
      max_D = stoi(argv[i+1]);
      cout << " -D " << max_D;
      i = i + 2;
    } else if(!strcmp(argv[i],"-m")){
      mismatch = stoi(argv[i+1]);
      cout << " -m " << mismatch;
      i = i + 2;
    } else if(!strcmp(argv[i],"-jf")){
      jf_file = argv[i+1];
      jf_file_exist = true;
      cout << " -jf " << jf_file;
      i = i + 2;
      bool_single = true;
    } else if(!strcmp(argv[i],"-jf1")){
      jf_file1 = argv[i+1];
      jf_file1_exist = true;
      cout << " -jf1 " << jf_file1;
      i = i + 2;
      bool_compare = true;
    } else if(!strcmp(argv[i],"-jf2")){
      jf_file2 = argv[i+1];
      jf_file2_exist = true;
      cout << " -jf2 " << jf_file2;
      i = i + 2;
      bool_compare = true;
    } else if(!strcmp(argv[i],"-dm")){
      dump_file = argv[i+1];
      dump_file_exist = true;
      cout << " -dm " << dump_file;
      i = i + 2;
      bool_single = true;
    } else if(!strcmp(argv[i],"-dm1")){
      dump_file1 = argv[i+1];
      dump_file1_exist = true;
      cout << " -dm1 " << dump_file1;
      i = i + 2;
      bool_compare = true;
    } else if(!strcmp(argv[i],"-dm2")){
      dump_file2 = argv[i+1];
      dump_file2_exist = true;
      cout << " -dm2 " << dump_file2;
      i = i + 2;
      bool_compare = true;
    } else if(!strcmp(argv[i],"-hs")){
      histo_file = argv[i+1];
      histo_file_exist = true;
      cout << " -hs " << histo_file;
      i = i + 2;
      bool_single = true;
    } else if(!strcmp(argv[i],"-hs1")){
      histo_file1 = argv[i+1];
      histo_file1_exist = true;
      cout << " -hs1 " << histo_file1;
      i = i + 2;
      bool_compare = true;
    } else if(!strcmp(argv[i],"-hs2")){
      histo_file2 = argv[i+1];
      histo_file2_exist = true;
      cout << " -hs2 " << histo_file2;
      i = i + 2;
      bool_compare = true;
    } else if(!strcmp(argv[i],"-rc")){
      read_coverage = stod(argv[i+1]);
      cout << " -rc " << read_coverage;
      i = i + 2;
    } else if(!strcmp(argv[i],"-time")){
      bool_stopwatch = true;
      cout << " -time ";
      i = i + 1;
    } else {
      cerr << "error : wrong option" << endl;
      cerr << argv[i] << "\t" << i << endl;
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

      void show_all(){
  cout << endl;
  cout << "  ****************************" << endl;
  cout << "  *  Run the whole pipeline  *" << endl;
  cout << "  ****************************" << endl;
}


void All::all_exe(int argc, char**argv){
  auto stopwatch_start = std::chrono::system_clock::now(); //stopwatch start
  show_all();
  if(option_all_parse(argc, argv) == -1){
    print_all_usage();
    return;
  }

  Extract extract;
  opt_all_to_extract(&extract);
  extract.extract_exe();
  kmer_for_cycle_find = extract.kmer_for_cycle_find;
  peak_to_opt_all(&extract.peak, &extract.peak1, &extract.peak2);  
  kmer_for_copy_estimate = extract.kmer_for_copy_estimate;
  stopwatch(stopwatch_start, bool_stopwatch);

  Cycle cycle;
  opt_all_to_cycle(&cycle);
  cycle.d = extract.d;
  cycle.cycle_exe();
  stopwatch(stopwatch_start, bool_stopwatch);
  
  Cluster cluster_t;
  repeat_type = 1;
  opt_all_to_cluster(&cluster_t);
  if(cycle.tandem_exist){
    cluster_t.cluster_exe();
  }
  stopwatch(stopwatch_start, bool_stopwatch);

  Intersperse intersperse;
  opt_all_to_intersperse(&intersperse);
  intersperse.kmer_id_file = cycle.fa;
  intersperse.kmer_align_file = cluster_t.kmer_align_file;
  intersperse.intersperse_exe();
  stopwatch(stopwatch_start, bool_stopwatch);

  Cluster cluster_i;
  repeat_type = 0;
  opt_all_to_cluster(&cluster_i);
  if(intersperse.intersperse_exist){
    cluster_i.cluster_exe();
  }
  stopwatch(stopwatch_start, bool_stopwatch);
}

void All::opt_all_to_extract(Extract *extract){ //bllのoptionをextractのoptionに移す
  if(d == 1){
    extract->files = files;
    extract->d = d;
  } else if(d == 0){
    extract->files1 = files1;
    extract->files2 = files2;
    extract->d = d;
  }
  extract->c1 = c1;
  extract->c1_string = c1_string;
  extract->c2 = c2;
  extract->t = t;
  extract->o = o;

  extract->jf_file_exist = jf_file_exist;
  extract->jf_file1_exist = jf_file1_exist;
  extract->jf_file2_exist = jf_file2_exist;
  extract->jf_file = jf_file;
  extract->jf_file1 = jf_file1;
  extract->jf_file2 = jf_file2;

  extract->dump_file_exist = dump_file_exist;
  extract->dump_file1_exist = dump_file1_exist;
  extract->dump_file2_exist = dump_file2_exist;
  extract->dump_file = dump_file;
  extract->dump_file1 = dump_file1;
  extract->dump_file2 = dump_file2;  

  extract->histo_file_exist = histo_file_exist;
  extract->histo_file1_exist = histo_file1_exist;
  extract->histo_file2_exist = histo_file2_exist;
  extract->histo_file = histo_file;
  extract->histo_file1 = histo_file1;
  extract->histo_file2 = histo_file2;  
}

void All::peak_to_opt_all(uint64_t *peak, uint64_t *peak1, uint64_t *peak2){
  p = *peak;
  p1 = *peak1;
  p2 = *peak2;
}

void All::opt_all_to_cycle(Cycle *cycle){
  cycle->p = p;
  cycle->p1 = p1;
  cycle->p2 = p2;
  cycle->fa = this->kmer_for_cycle_find;
  cycle->d = d;
  if(cycle->d == 1){
    cycle->fq = files[0]; //ダウンサンプリングをするのでfileの一つめだけを渡す
  } else if(cycle->d == 0){
    cycle->fq = files1[0]; //ダウンサンプリングをするのでfileの一つめだけを渡す
  }
  cycle->t = t;
  cycle->o = o + "_T";
  cycle->max_l = max_l;
  cycle->max_n = max_n;
  cycle->max_d = max_d;
  cycle->read_coverage = read_coverage;
}


void All::opt_all_to_cluster(Cluster *cluster){
  cluster->f1 = kmer_for_copy_estimate;
  string output;
  if(repeat_type == 1){
    output = this->o + "_T";
  } else if(repeat_type == 0){
    output = this->o + "_I";
  } else {
    cerr << "error in bll.cpp" << endl;
  }
  
  cluster->mismatch = this->mismatch;

  cluster->f2 = output + "_repeat_num";
  cluster->f3 = output + "_repeat_num_min";
  cluster->r_t = repeat_type;
  cluster->o = output;

  cluster->p = this->p;
  cluster->p1 = this->p1;
  cluster->p2 = this->p2;
  cluster->t = this->t;  
  cluster->d = this->d;
}

void All::opt_all_to_intersperse(Intersperse *intersperse){
  intersperse->p = this->p;
  intersperse->p1 = this->p1;
  intersperse->p2 = this->p2;
  //  intersperse->fa = this->o + "_" + to_string((long long)this->k) + "mer_" + to_string((long long)this->c2) + ".fa";
  //  intersperse->fa = this->kmer_for_path_find;
  //  intersperse->b = this->o + "_blst.blastn";
  intersperse->d = this->d;
  if(intersperse->d == 1){
    intersperse->fq = this->files[0]; //ダウンサンプリングをするのでfileの一つめだけを渡す
  } else if(intersperse->d == 0){
    intersperse->fq = this->files1[0]; //ダウンサンプリングをするのでfileの一つめだけを渡す
  }
  intersperse->t = this->t;
  intersperse->o = this->o + "_I";
  intersperse->max_l = this->max_L;
  intersperse->max_n = this->max_N;
  intersperse->max_d = this->max_D;
  intersperse->read_coverage = this->read_coverage;
}
