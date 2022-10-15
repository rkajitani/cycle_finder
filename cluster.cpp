#include "cluster.h"
//#include "clustering.h"
//#include "kmer_align.h"
//#include "blast_to_copynum.h"
//#include "clst_to_family.h"
#include "common.h"

Cluster::Cluster(){
  this->o = "out";
  this->t = 1;
  this->mismatch = 3;
  this->r_t = 1;
  this->d = 0;  

  this->bool_single = false;
  this->bool_compare = false;
}

void Cluster::print_cluster_usage(void){
  cerr << "cluster usage: cycle_finder cluster [option] " << endl;
  cerr << "------------------------------------------------" << endl;
  cerr << "-f1: k-mer file " << endl;
  cerr << "-f2: tsv file " << endl;
  cerr << "-f3: tsv file(min) " << endl;
  cerr << "-c : single->1, compare->0 [bool] (defualt:" << this->d << ")" << endl;
  cerr << "-i : tandem->1, intersperse->0 [bool] (defualt:" << this->r_t << ")" << endl;
  cerr << "-p : peak  [int] " << endl;
  cerr << "-p1: peak1 [int] " << endl;
  cerr << "-p2: peak2 [int] " << endl;
  cerr << "-t : num_threads[int](default:" << t << ")" << endl;
  cerr << "-m : threshold mismatch of kmer alignment for copy estimation;less than or equal to 3 is recommended(default:" << mismatch << ")" << endl;
  cerr << "-o : filename of output [int] (defualt:" << this->o << ")" << endl;
  cerr << "-------------------------------------------------" << endl;
}

int Cluster::option_cluster_parse(int argc, char**argv){
  cout << PROGRAM_NAME << " cluster";
  for(int i = 2; i < argc;){
    if(!strcmp(argv[i], "-f1")){
      this->f1 = argv[i+1];
      cout << " -f1 " << this->f1;
      i = i + 2;
    } else if(!strcmp(argv[i], "-f2")){
      this->f2 = argv[i+1];
      cout << " -f2 " << this->f2;
      i = i + 2;
    } else if(!strcmp(argv[i], "-f3")){
      this->f3 = argv[i+1];
      cout << " -f3 " << this->f3;
      i = i + 2;
    } else if(!strcmp(argv[i], "-c")){
      this->d = stoi(argv[i+1]);
      cout << " -c " << this->d;
      i = i + 2;
    } else if(!strcmp(argv[i], "-i")){
      this->r_t = stoi(argv[i+1]); //repeat_type
      cout << " -i " << this->r_t;
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
    } else if(!strcmp(argv[i], "-t")){
      this->t = stoi(argv[i+1]);
      cout << " -t " << this->t;
      i = i + 2;
    } else if(!strcmp(argv[i], "-m")){
      this->mismatch = stoi(argv[i+1]);
      cout << " -m " << this->mismatch;
      i = i + 2;
    } else if(!strcmp(argv[i],"-o")){
      this->o = argv[i+1];
      cout << " -o " << this->o;
      i = i + 2;
    } else {
      cerr << "error : wrong option" << endl;
      cerr << argv[i] << endl;
      return -1;
    }
  }
  cout << endl;
  if((bool_single == true && bool_compare == true) || (bool_single == false && bool_compare == false)){
    cerr << "error : input file correctly" << endl;
    return -1;
  }

  if(this->r_t == 1){ //tandemなら
    this->o += "_T";
  } else if(this->r_t == 0){ //intersperseなら
    this->o += "_I";
  }
  return 0;
}

void Cluster::cluster_tandem(){
  string repeat_num_tandem = this->o + "_repeat_num_tandem.fa";
  tandem(this->f2, repeat_num_tandem); //tsvファイルからtandemのfastaをつくる

  string repeat_num_tandem_clst = this->o + "_repeat_num_tandem_clst";
  cdhit(repeat_num_tandem, repeat_num_tandem_clst, CDHIT_COVERAGE_THR_TANDEM); //fasta, num_threads, output, coverage


  string repeat_num_tandem_clst_max = this->o + "_repeat_num_tandem_clst_max";
  max_copy_fasta(this->f2, repeat_num_tandem_clst + ".clstr", repeat_num_tandem_clst_max); //repeat_num_file, clustering1, output

  string repeat_num_not_tandem_clst = this->o + "_repeat_num_not_tandem_clst";
  cdhit(repeat_num_tandem_clst_max, repeat_num_not_tandem_clst , CDHIT_COVERAGE_THR_NOT_TANDEM);


  string repeat_cluster = this->blast_target;
  cdhit_integrage(this->f2, repeat_num_tandem_clst_max + ".clstr", repeat_num_not_tandem_clst + ".clstr", repeat_cluster);

  system(("mv " + repeat_cluster + "_1 " + this->o + "_clst.fa").c_str());
}

void Cluster::cluster_intersperse(){
  system( ("cat " + this->f2 + " | awk \'{print \">\"NR\"_\"length($1)\"\\n\"$1;}\' > " + this->f2 + ".fa").c_str() );
  //    cd-hit-est -i ${repeat_fasta}_result.fa -n 3 -c $cd_prc -aL $cd_prs2 -o ${repeat_fasta}_result_cdhit -T 0 -r 1 -g 1 -G 0 -M 2000;
  cdhit(this->f2 + ".fa", this->f2 + "_cdhit" , CDHIT_COVERAGE_THR_NOT_TANDEM);
  clst_to_family(this->f2 + ".fa", this->f2 + "_cdhit.clstr", this->o + "_blst"); //tsvファイル, cdhitの結果.clstrファイル, 出力
  system( ("cp " + this->o + "_blst " +  this->o + "_clst.fa" ).c_str() );					     
}

void Cluster::cluster_cdhit(){
  if(this->r_t == 1){ //tandem repeat
    cluster_tandem();
  } else if(this->r_t == 0){ //interspersed repeat
    cluster_intersperse();
  }
}

void Cluster::make_output_file(void){
  this->repeat_num_tandem = this->o + "_repeat_num_tandem.fa";
  this->repeat_num_tandem_clst = this->o + "_repeat_num_tandem_clst";
  this->repeat_num_tandem_clst_max = this->o + "_repeat_num_tandem_clst_max";
  this->repeat_num_not_tandem_clst = this->o + "_repeat_num_not_tandem_clst";
  this->repeat_cluster = this->o + "_blst";
  this->kmer_fasta = this->f1;
  this->blast_target = this->o + "_blst";
  this->kmer_align_file = this->o + "_blst.blastn";
  this->represent_fa = this->o + "_blst_tmp";
  this->blast_self = this->o + "_blst.blastn_self";
  this->blst_self = this->o + "_blst_self";
  this->blst_family = this->o + "_blst_family";
  this->blst_family_clst = this->o + "_blst_family_clst";
  this->blst_family_clst_fa = this->o + "_blst_family_clst.fa";
  this->blst_family_clst_element = this->o + "_blst_family_clst.element";
  this->final_output_tsv = this->o + ".tsv";
  this->final_output_fa = this->o + ".fa";
  this->final_output_element = this->o + ".element";

  this->cmd_blast_self = "blastn -db " + this->represent_fa + " -query " + this->represent_fa + " -outfmt 6 -out " + this->blast_self 
    + " -task blastn-short -dust no -num_threads " + to_string(static_cast<long long>(this->t));
  this->len_thr = "0.6";
  this->ide_thr = "80";

  this->filter_blast_self = "cat " + this->o + "_blst.blastn_self |awk \'$1!=$2{split($1,a,\"_\");split($2,b,\"_\");if(a[2]>b[2]&&(a[3]>b[3]&&$3*$4/100>a[3]*" 
    + len_thr + "&&$3>" + ide_thr + "||a[3]<=b[3] && $3*$4/100>b[3]*" + len_thr + "&&$3>" + ide_thr 
    + "))print;}\' |sort -k1,2 -u |awk \'{split($1,a,\"_\");print a[2]\"\\t\"$0;}\' |sort -k1nr > " + this->blst_self;
  this->inherit_tsv = "cat " + this->o + "_blst_family | awk  \'split($1,a,\"_\"){print \"Family_\"a[2]\"_\"a[3]\"\\t\"$2\"\\t\"$3\"\\t\"$4\"\\t\"$5;}\' > " + this->o + "_blst_family_clst";

  this->inherit_fa = "cat " + this->o + "_blst_tmp"
      + " |grep \"*\" -A1|grep -v \"-\"|awk \'{if(substr($1,1,1)==\">\"){split($1,a,\"_\");print \">Family_\"a[2]\"_\"a[3]\"\\t\"$2\"\\t\"$3\"\\t\"$4\"\\t\"$5;} else{print;}}' > " 
    + this->o + "_blst_family_clst.fa";  
}

void Cluster::copy_estimation(){
  cout << "\t" << "alignmenting k-mer..." << endl;
  kmer_align();
  blast_to_copynum();
}

void Cluster::cluster_blast(){
  system(("grep '*' -A1 " + this->o + "_blst " +  "|grep -v '-' > " + represent_fa).c_str());
  makeblastdb(represent_fa, this->o + ".log");
  cerr << "CMD: " << cmd_blast_self << endl;
  system(cmd_blast_self.c_str());
  system(("rm " + this->o + "_blst_tmp.n*").c_str());
  cerr << "CMD: " << filter_blast_self << endl;
  system(filter_blast_self.c_str());
  file_write(this->o + ".log", filter_blast_self);
  //check_file_size
  if(file_exist(this->blst_self) == 0){
    cout << "no family clustering" << endl;
    cerr << "CMD: " << inherit_tsv << endl;
    system(inherit_tsv.c_str());
    file_write(this->o + ".log", inherit_tsv);
    cerr << "CMD: " << inherit_fa << endl;
    system(inherit_fa.c_str());
    file_write(this->o + ".log", inherit_fa);
  } else {
    cluster_by_blast(this->f2, this->blst_family, this->blst_self, this->o + "_clst.fa", this->blst_family_clst,this->r_t);
  }
  system(("cp " + this->blst_family_clst + " " + this->final_output_tsv).c_str());
  system(("cp " + this->blst_family_clst_fa + " " + this->final_output_fa).c_str());
  if(file_exist(this->blst_self) != 0){ //when family clustering is available
    system(("cp " + this->blst_family_clst_element + " " + this->final_output_element).c_str());
  }
}

void Cluster::remove_cluster(){
  remove(repeat_num_tandem.c_str());
  remove(repeat_num_tandem_clst.c_str());
  remove((repeat_num_tandem_clst + ".clstr").c_str());
  remove(repeat_num_tandem_clst_max.c_str());
  remove((repeat_num_tandem_clst_max + ".clstr").c_str());
  remove(repeat_num_not_tandem_clst.c_str());
  remove((repeat_num_not_tandem_clst + ".clstr").c_str());
  remove((this->f2 + "_cdhit").c_str());
  remove((this->f2 + "_cdhit.clstr").c_str());
  remove((this->represent_fa).c_str());
  remove((this->o + ".log").c_str());
  remove((this->repeat_cluster).c_str());
  remove((this->blast_target).c_str());
  remove((this->blast_self).c_str());
  remove((this->blst_family).c_str());
  remove((this->blst_family_clst).c_str());
  remove((this->blst_family_clst_fa).c_str());
  remove((this->blst_family_clst_element).c_str());
  remove((this->blst_self).c_str());
  
  if(this->r_t == 0) remove((this->kmer_align_file).c_str()); //remove intersperse kmer align file
}

void show_cluster(int r_t){
  cout << endl;
  if(r_t == 1){
    cout << "  +--------------------------+" << endl;
    cout << "  |  Cluster tandem repeats  |" << endl;
    cout << "  +--------------------------+" << endl;
  } else if(r_t == 0){
    cout << "  +--------------------------------+" << endl;
    cout << "  |  Cluster interspersed repeats  |" << endl;
    cout << "  +--------------------------------+" << endl;
  }
}

void Cluster::cluster_exe(void){
  show_cluster(this->r_t);
  make_output_file();
  cluster_cdhit();
  copy_estimation();
  cluster_blast();  
  remove_cluster();
}

void Cluster::cdhit(string fasta_file, string file_out, double align_cov){
  int word_length = 3;
  int accurate = 1;
  int global = 0;
  int memory_lim = 0;
  string root_path = ROOT_PATH;

  string cmd = root_path + "/cd-hit-est -i " + fasta_file  //クラスタリングする配列
    + " -n " + to_string(static_cast<long long>(word_length))  //word_length default:10
    + " -c " + to_string(static_cast<long double>(CDHIT_IDENTITY_THR))  //sequence identity
    + " -aL " + to_string(static_cast<long double>(align_cov))  //alignment coverage for the longer sequence
    + " -o " + file_out  //output
    + " -T " + to_string(static_cast<long long>(this->t))  //num_threads
    + " -g " + to_string(static_cast<long long>(accurate))  //accurate but slow mode
    + " -G " + to_string(static_cast<long long>(global))  //use global identity
    + " -M " + to_string(static_cast<long long>(memory_lim)) //memory limit
    + " >> " + this->o + ".log";
  system(("echo " + cmd).c_str());
  system(cmd.c_str());
}

//kmer_alignment using blast
//kmer_align_blast(){
  //makeblastdb k-mer
  //  makeblastdb(kmer_fasta_id,this->o+".log");

  //blast k-mer
  //  string cmd_blast_kmer = "blastn -db " + this->f1 + "_id -query " + this->o + "_blst -outfmt 6 -out " 
  //    + this->o + "_blst.blastn -max_hsps 1 -dust no -num_threads " + to_string(static_cast<long long>(this->t)) 
  //    + " -evalue 1000 -task blastn-short -max_target_seqs 100000";
  //  system(cmd_blast_kmer.c_str());
  //  file_write(this->o + ".log", cmd_blast_kmer);
  //filtering
  //string cmd_filter = "cat " + this->o + "_blst.blastn | awk \'$4>K_MER*0.75{print;}\' > " + this->o + "_blst.blastn_filter";
  //  string cmd_filter = "cat " + this->o + "_blst.blastn | sort -k1 > " + this->o + "_blst.blastn_filter";
  //  system(cmd_filter.c_str());
  //  file_write(this->o + ".log", cmd_filter);

  //index削除
  //  system(("rm " + this->f1 + "_id.*").c_str());


//}
